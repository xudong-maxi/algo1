/**
 * @file    jdps_test.c
 * @brief   Host test: runs jdps_process() on vectors exported by
 *          cs_agc/export_vectors.py and compares with the Python reference.
 *
 * Usage: jdps_test vectors.txt
 * Pass criteria per case:
 *   status and v_valid identical, |score_c - score_py| < 1e-3 * max(1, |score|)
 *   |v_c - v_py| (applied velocity) < 0.01 m/s
 *   max wrapped |psi_c - psi_py|    < 0.01 rad
 *   ||iq_c - iq_py|| / ||iq_py||    < 1e-3
 */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "../jdps.h"

#define TOL_V_MPS     0.01
#define TOL_PSI_RAD   0.01
#define TOL_IQ_REL    1e-3

typedef struct {
    char             name[64];
    jdps_cfg_t    cfg;
    uint8_t          num_ant, num_steps, num_seg;
    uint8_t          chan[JDPS_MAX_STEPS];
    uint8_t          seg[JDPS_MAX_STEPS];
    float            time_s[JDPS_MAX_STEPS];
    complex        iq[JDPS_MAX_ANT][JDPS_MAX_STEPS];
    int              exp_status;
    int              exp_valid;
    double           exp_score;
    double           exp_v_est;
    double           exp_v;
    double           exp_psi[JDPS_MAX_ANT][JDPS_MAX_SEG];
    double           exp_iq[JDPS_MAX_ANT][JDPS_MAX_STEPS][2];
} test_case_t;

static int expect_tag(FILE *fp, const char *tag)
{
    char buf[32];
    return fscanf(fp, "%31s", buf) == 1 && strcmp(buf, tag) == 0;
}

/** Reads one case; returns 0 at end of file. */
static int read_case(FILE *fp, test_case_t *tc)
{
    float vmax, vguard, vstep, vscore;
    unsigned orders, pit, rit, perant, A, N, K, u;

    if (fscanf(fp, " case %63s", tc->name) != 1) {
        return 0;
    }
    if (!expect_tag(fp, "cfg") ||
        fscanf(fp, "%f %f %f %f %u %u %u %u", &vmax, &vguard, &vstep, &vscore, &orders, &pit, &rit,
               &perant) != 8 ||
        !expect_tag(fp, "dims") || fscanf(fp, "%u %u %u", &A, &N, &K) != 3) {
        return -1;
    }
    tc->cfg = jdps_default_cfg();      /* frequency plan: BLE CS defaults */
    tc->cfg.v_max_mps = vmax;
    tc->cfg.v_guard_mps = vguard;
    tc->cfg.v_step_mps = vstep;
    tc->cfg.min_v_score = vscore;
    tc->cfg.num_orders = (uint8_t)orders;
    tc->cfg.power_iter_num = (uint8_t)pit;
    tc->cfg.refine_iter_num = (uint8_t)rit;
    tc->cfg.per_ant_phase = (uint8_t)perant;
    tc->num_ant = (uint8_t)A;
    tc->num_steps = (uint8_t)N;
    tc->num_seg = (uint8_t)K;

    if (!expect_tag(fp, "chan")) return -1;
    for (unsigned n = 0; n < N; n++) { if (fscanf(fp, "%u", &u) != 1) return -1; tc->chan[n] = (uint8_t)u; }
    if (!expect_tag(fp, "se")) return -1;
    for (unsigned n = 0; n < N; n++) { if (fscanf(fp, "%u", &u) != 1) return -1; tc->seg[n] = (uint8_t)u; }
    if (!expect_tag(fp, "time")) return -1;
    for (unsigned n = 0; n < N; n++) { if (fscanf(fp, "%f", &tc->time_s[n]) != 1) return -1; }
    if (!expect_tag(fp, "iq")) return -1;
    for (unsigned a = 0; a < A; a++)
        for (unsigned n = 0; n < N; n++)
            if (fscanf(fp, "%f %f", &tc->iq[a][n].r, &tc->iq[a][n].i) != 2) return -1;
    if (!expect_tag(fp, "exp_status") ||
        fscanf(fp, "%d %d %lf %lf", &tc->exp_status, &tc->exp_valid, &tc->exp_score, &tc->exp_v_est) != 4)
        return -1;
    if (!expect_tag(fp, "exp_v") || fscanf(fp, "%lf", &tc->exp_v) != 1) return -1;
    if (!expect_tag(fp, "exp_psi")) return -1;
    for (unsigned a = 0; a < A; a++)
        for (unsigned k = 0; k < K; k++)
            if (fscanf(fp, "%lf", &tc->exp_psi[a][k]) != 1) return -1;
    if (!expect_tag(fp, "exp_iq")) return -1;
    for (unsigned a = 0; a < A; a++)
        for (unsigned n = 0; n < N; n++)
            if (fscanf(fp, "%lf %lf", &tc->exp_iq[a][n][0], &tc->exp_iq[a][n][1]) != 2) return -1;
    return 1;
}

static double wrap_pi(double x)
{
    return x - 2.0 * M_PI * floor((x + M_PI) / (2.0 * M_PI));
}

int main(int argc, char **argv)
{
    static test_case_t    tc;
    static jdps_work_t work;
    jdps_result_t      res;
    int num_cases = 0, num_fail = 0;
    double total_us = 0.0;

    FILE *fp = fopen(argc > 1 ? argv[1] : "vectors.txt", "r");
    if (fp == NULL) {
        perror("open vectors");
        return 2;
    }
    printf("%-26s %9s %9s %6s %6s %5s %10s %10s  %s\n", "case", "v_c", "v_py", "score", "valid", "st",
           "psi_err", "iq_relerr", "result");
    int rc;
    while ((rc = read_case(fp, &tc)) == 1) {
        jdps_meas_t meas = {
            .num_ant = tc.num_ant, .num_steps = tc.num_steps, .num_seg = tc.num_seg,
            .chan_idx = tc.chan, .seg_idx = tc.seg, .step_time_s = tc.time_s, .iq = tc.iq,
        };
        clock_t t0 = clock();
        memset(&res, 0, sizeof(res));
        jdps_status_t st = jdps_process(&tc.cfg, &meas, &work, &res);
        total_us += 1e6 * (double)(clock() - t0) / CLOCKS_PER_SEC;

        double psi_err = 0.0, err2 = 0.0, ref2 = 0.0;
        for (unsigned a = 0; a < tc.num_ant; a++) {
            for (unsigned k = 0; k < tc.num_seg; k++) {
                double e = fabs(wrap_pi(res.psi_rad[a][k] - tc.exp_psi[a][k]));
                psi_err = e > psi_err ? e : psi_err;
            }
            for (unsigned n = 0; n < tc.num_steps; n++) {
                double dr = tc.iq[a][n].r - tc.exp_iq[a][n][0];
                double di = tc.iq[a][n].i - tc.exp_iq[a][n][1];
                err2 += dr * dr + di * di;
                ref2 += tc.exp_iq[a][n][0] * tc.exp_iq[a][n][0] + tc.exp_iq[a][n][1] * tc.exp_iq[a][n][1];
            }
        }
        double iq_rel = sqrt(err2 / (ref2 > 0 ? ref2 : 1.0));
        double score_tol = 1e-3 * (fabs(tc.exp_score) > 1.0 ? fabs(tc.exp_score) : 1.0);
        int ok = (int)st == tc.exp_status && res.v_valid == tc.exp_valid &&
                 fabs(res.v_score - tc.exp_score) < score_tol &&
                 fabs(res.v_mps - tc.exp_v) < TOL_V_MPS &&
                 psi_err < TOL_PSI_RAD && iq_rel < TOL_IQ_REL;
        printf("%-26s %+9.4f %+9.4f %6.1f %6d %5d %10.2e %10.2e  %s\n", tc.name, res.v_mps, tc.exp_v,
               res.v_score, res.v_valid, (int)st, psi_err, iq_rel, ok ? "PASS" : "FAIL");
        num_cases++;
        num_fail += !ok;
    }
    fclose(fp);
    if (rc < 0) {
        printf("vector file parse error after %d cases\n", num_cases);
        return 2;
    }
    printf("\n%d cases, %d failed; mean host time %.1f us/case (4 ant, x86, not MCU)\n",
           num_cases, num_fail, num_cases ? total_us / num_cases : 0.0);
    return num_fail ? 1 : 0;
}
