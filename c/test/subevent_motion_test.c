/*
 * subevent_motion_alg 的主机回归测试：读取 cs_agc/export_vectors.py 导出的测试向量，
 * 按工程的数据组织方式（channel_select_t + 按信道号排列的 IQ，每次只放一路 IQ）调用模块，
 * 与 Python 参考实现逐点对比。
 *
 * 用法：subevent_motion_test vectors.txt
 * 通过条件（每个用例）：
 *   返回状态、速度可信标志一致；分数相对误差 < 1e-3；
 *   |速度差| < 0.01 m/s；subevent 相位差 < 0.01 rad；补偿后 IQ 相对误差 < 1e-3
 */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "common_util.h"
#include "securec.h"
#include "subevent_motion_alg.h"

#define TOL_SPEED       0.01
#define TOL_PHASE       0.01
#define TOL_IQ_REL      1e-3
#define MAX_STEP_NUM    ALG_CHANNEL_NUM
#define STATUS_NO_SIGNAL 3      // 与 export_vectors.py 中的状态编号一致

typedef struct {
    char name[64];
    float max_speed, speed_guard, speed_resolution, min_score;
    unsigned orders, per_path;
    unsigned path_num, step_num, subevent_num;
    uint8_t channel[MAX_STEP_NUM];
    uint8_t subevent[MAX_STEP_NUM];
    double time_s[MAX_STEP_NUM];
    complex iq[SUBEVENT_MOTION_MAX_PATH_NUM][MAX_STEP_NUM];       // 跳频顺序
    int exp_status, exp_valid;
    double exp_score, exp_speed_est, exp_speed;
    double exp_phase[SUBEVENT_MOTION_MAX_PATH_NUM][SUBEVENT_MOTION_MAX_SUBEVENT_NUM];
    double exp_iq[SUBEVENT_MOTION_MAX_PATH_NUM][MAX_STEP_NUM][2];
} TestCase;

static int expect_tag(FILE* fp, const char* tag)
{
    char buf[32];
    return (fscanf(fp, "%31s", buf) == 1) && (strcmp(buf, tag) == 0);
}

// 读取一个用例，返回 1 成功，0 文件结束，-1 格式错误
static int read_case(FILE* fp, TestCase* tc)
{
    unsigned power_iter, refine, u;
    if (fscanf(fp, " case %63s", tc->name) != 1) {
        return 0;
    }
    if (!expect_tag(fp, "cfg") ||
        fscanf(fp, "%f %f %f %f %u %u %u %u", &tc->max_speed, &tc->speed_guard, &tc->speed_resolution,
               &tc->min_score, &tc->orders, &power_iter, &refine, &tc->per_path) != 8 ||
        !expect_tag(fp, "dims") || fscanf(fp, "%u %u %u", &tc->path_num, &tc->step_num, &tc->subevent_num) != 3 ||
        tc->path_num > SUBEVENT_MOTION_MAX_PATH_NUM || tc->step_num > MAX_STEP_NUM ||
        tc->subevent_num > SUBEVENT_MOTION_MAX_SUBEVENT_NUM) {
        return -1;
    }
    if (!expect_tag(fp, "chan")) return -1;
    for (unsigned n = 0; n < tc->step_num; n++) { if (fscanf(fp, "%u", &u) != 1) return -1; tc->channel[n] = (uint8_t)u; }
    if (!expect_tag(fp, "se")) return -1;
    for (unsigned n = 0; n < tc->step_num; n++) { if (fscanf(fp, "%u", &u) != 1) return -1; tc->subevent[n] = (uint8_t)u; }
    if (!expect_tag(fp, "time")) return -1;
    for (unsigned n = 0; n < tc->step_num; n++) { if (fscanf(fp, "%lf", &tc->time_s[n]) != 1) return -1; }
    if (!expect_tag(fp, "iq")) return -1;
    for (unsigned p = 0; p < tc->path_num; p++)
        for (unsigned n = 0; n < tc->step_num; n++)
            if (fscanf(fp, "%f %f", &tc->iq[p][n].r, &tc->iq[p][n].i) != 2) return -1;
    if (!expect_tag(fp, "exp_status") ||
        fscanf(fp, "%d %d %lf %lf", &tc->exp_status, &tc->exp_valid, &tc->exp_score, &tc->exp_speed_est) != 4)
        return -1;
    if (!expect_tag(fp, "exp_v") || fscanf(fp, "%lf", &tc->exp_speed) != 1) return -1;
    if (!expect_tag(fp, "exp_psi")) return -1;
    for (unsigned p = 0; p < tc->path_num; p++)
        for (unsigned k = 0; k < tc->subevent_num; k++)
            if (fscanf(fp, "%lf", &tc->exp_phase[p][k]) != 1) return -1;
    if (!expect_tag(fp, "exp_iq")) return -1;
    for (unsigned p = 0; p < tc->path_num; p++)
        for (unsigned n = 0; n < tc->step_num; n++)
            if (fscanf(fp, "%lf %lf", &tc->exp_iq[p][n][0], &tc->exp_iq[p][n][1]) != 2) return -1;
    return 1;
}

// 由每一步的时间戳还原工程的时序描述：time_per_channel（us）与 subevent 间隔 t_mes（us）
static int build_channel_select(const TestCase* tc, channel_select_t* cfg, uint8_t* hop, uint8_t* per_subevent,
                                uint16_t* time_per_channel)
{
    long t_us[MAX_STEP_NUM];
    for (unsigned n = 0; n < tc->step_num; n++) {
        t_us[n] = lround((tc->time_s[n] - tc->time_s[0]) * 1e6);
    }
    memset(per_subevent, 0, SUBEVENT_MOTION_MAX_SUBEVENT_NUM);
    long step_us = 0;
    for (unsigned n = 0; n + 1 < tc->step_num; n++) {
        if (tc->subevent[n + 1] == tc->subevent[n]) {
            step_us = t_us[n + 1] - t_us[n];
            break;
        }
    }
    long t_mes = -1;
    hop[0] = tc->channel[0];                                // mode0 占位，算法不使用
    for (unsigned n = 0; n < tc->step_num; n++) {
        hop[n + 1] = tc->channel[n];
        per_subevent[tc->subevent[n]]++;
        bool last_in_subevent = (n + 1 == tc->step_num) || (tc->subevent[n + 1] != tc->subevent[n]);
        time_per_channel[n] = (uint16_t)(last_in_subevent ? step_us : t_us[n + 1] - t_us[n]);
        if (last_in_subevent && (n + 1 < tc->step_num)) {
            long gap = t_us[n + 1] - t_us[n] - step_us;
            if ((t_mes >= 0) && (gap != t_mes)) {
                return -1;                                  // 各 subevent 间隔不相等，无法用单个 t_mes 表示
            }
            t_mes = gap;
        }
    }
    cfg->ch_num = (uint8_t)(tc->step_num + 1);
    cfg->subevent_num = (uint8_t)tc->subevent_num;
    cfg->t_mes = (uint16_t)(t_mes < 0 ? 0 : t_mes);
    cfg->ch_hop_orders = hop;
    cfg->antenna_path_permutation_index = NULL;
    cfg->ch_num_per_subevent = per_subevent;
    cfg->time_per_channel = time_per_channel;
    return 0;
}

static double wrap_pi(double x)
{
    return x - 2.0 * PI * floor((x + PI) / (2.0 * PI));
}

// 按工程方式运行：每次只把一路 IQ（按信道号排列）放进 iq_buf
static int run_case(TestCase* tc, SubeventMotionCtx* ctx, SubeventMotionResult* res,
                    complex out[SUBEVENT_MOTION_MAX_PATH_NUM][ALG_CHANNEL_NUM])
{
    static complex iq_buf[ALG_CHANNEL_NUM];
    uint8_t hop[MAX_STEP_NUM + 1];
    uint8_t per_subevent[SUBEVENT_MOTION_MAX_SUBEVENT_NUM];
    uint16_t time_per_channel[MAX_STEP_NUM];
    channel_select_t cfg;

    // 外部存储中的 IQ：按信道号排列，未测量的信道为 0
    for (unsigned p = 0; p < tc->path_num; p++) {
        memset(out[p], 0, sizeof(out[p]));
        for (unsigned n = 0; n < tc->step_num; n++) {
            out[p][tc->channel[n]] = tc->iq[p][n];
        }
    }
    if (build_channel_select(tc, &cfg, hop, per_subevent, time_per_channel) != 0) {
        return -1;
    }
    errcode_t err = subevent_motion_init(ctx, &cfg, (uint8_t)tc->path_num);
    if (err != ERRCODE_RANGING_ALG_SUCCESS) return -1;
    for (unsigned p = 0; p < tc->path_num; p++) {
        memcpy(iq_buf, out[p], sizeof(iq_buf));
        if (subevent_motion_add_speed_path(ctx, iq_buf) != ERRCODE_RANGING_ALG_SUCCESS) return -1;
    }
    err = subevent_motion_solve_speed(ctx);
    if (err == ERRCODE_RANGING_ALG_NOT_ENOUGH_IQ) return STATUS_NO_SIGNAL;
    if (err != ERRCODE_RANGING_ALG_SUCCESS) return -1;
    for (unsigned p = 0; p < tc->path_num; p++) {
        memcpy(iq_buf, out[p], sizeof(iq_buf));
        if (subevent_motion_add_phase_path(ctx, (uint8_t)p, iq_buf) != ERRCODE_RANGING_ALG_SUCCESS) return -1;
    }
    if (subevent_motion_solve_phase(ctx, res) != ERRCODE_RANGING_ALG_SUCCESS) return -1;
    for (unsigned p = 0; p < tc->path_num; p++) {
        memcpy(iq_buf, out[p], sizeof(iq_buf));
        if (subevent_motion_iq_compensation(ctx, iq_buf) != ERRCODE_RANGING_ALG_SUCCESS) return -1;
        memcpy(out[p], iq_buf, sizeof(iq_buf));
    }
    return 0;
}

int main(int argc, char** argv)
{
    static TestCase tc;
    static SubeventMotionCtx ctx;
    static complex out[SUBEVENT_MOTION_MAX_PATH_NUM][ALG_CHANNEL_NUM];
    SubeventMotionResult res;
    int case_num = 0, fail_num = 0;
    double total_us = 0.0;

    FILE* fp = fopen(argc > 1 ? argv[1] : "vectors.txt", "r");
    if (fp == NULL) {
        perror("open vectors");
        return 2;
    }
    printf("PAIR_ORDER %d: sizeof(SubeventMotionCtx) = %zu B (MAX_PATH %d, MAX_SUBEVENT %d), IQ buffer = 1 path\n\n",
           SUBEVENT_MOTION_PAIR_ORDER, sizeof(SubeventMotionCtx), SUBEVENT_MOTION_MAX_PATH_NUM,
           SUBEVENT_MOTION_MAX_SUBEVENT_NUM);
    printf("%-26s %9s %9s %6s %6s %5s %10s %10s  %s\n", "case", "v_c", "v_py", "score", "valid", "st",
           "phase_err", "iq_relerr", "result");
    int rc;
    while ((rc = read_case(fp, &tc)) == 1) {
        if (tc.orders != SUBEVENT_MOTION_PAIR_ORDER || tc.per_path != 0) {
            continue;                                       // 用例的配对阶数与本次编译不一致
        }
        memset(&res, 0, sizeof(res));
        clock_t t0 = clock();
        int status = run_case(&tc, &ctx, &res, out);
        total_us += 1e6 * (double)(clock() - t0) / CLOCKS_PER_SEC;

        double phase_err = 0.0, err2 = 0.0, ref2 = 0.0;
        for (unsigned k = 0; k < tc.subevent_num; k++) {
            double e = fabs(wrap_pi(res.subevent_phase[k] - tc.exp_phase[0][k]));
            phase_err = e > phase_err ? e : phase_err;
        }
        for (unsigned p = 0; p < tc.path_num; p++) {
            for (unsigned n = 0; n < tc.step_num; n++) {
                complex c = out[p][tc.channel[n]];
                double dr = c.r - tc.exp_iq[p][n][0];
                double di = c.i - tc.exp_iq[p][n][1];
                err2 += dr * dr + di * di;
                ref2 += tc.exp_iq[p][n][0] * tc.exp_iq[p][n][0] + tc.exp_iq[p][n][1] * tc.exp_iq[p][n][1];
            }
        }
        double iq_rel = sqrt(err2 / (ref2 > 0 ? ref2 : 1.0));
        double score_tol = 1e-3 * (fabs(tc.exp_score) > 1.0 ? fabs(tc.exp_score) : 1.0);
        int ok = (status == tc.exp_status) && (res.speed_valid == tc.exp_valid) &&
                 (fabs(res.speed_score - tc.exp_score) < score_tol) && (fabs(res.speed - tc.exp_speed) < TOL_SPEED) &&
                 (phase_err < TOL_PHASE) && (iq_rel < TOL_IQ_REL);
        printf("%-26s %+9.4f %+9.4f %6.1f %6d %5d %10.2e %10.2e  %s\n", tc.name, res.speed, tc.exp_speed,
               res.speed_score, res.speed_valid, status, phase_err, iq_rel, ok ? "PASS" : "FAIL");
        case_num++;
        fail_num += !ok;
    }
    fclose(fp);
    if (rc < 0) {
        printf("vector file parse error after %d cases\n", case_num);
        return 2;
    }
    printf("\n%d cases, %d failed; mean host time %.1f us/case (x86, not MCU)\n", case_num, fail_num,
           case_num ? total_us / case_num : 0.0);
    return fail_num ? 1 : 0;
}
