/*
 * 用固件打印的 [jdps_dump] 输入（格式见 cs_agc/replay_dump.py）在主机上运行 subevent_motion_alg，
 * 打印与固件相同格式的 [jdps] 结果行，用于和固件结果、Python 参考实现逐项比对。
 *
 * 用法：replay_dump firmware.log
 * 可以把工程中移植后的 subevent_motion_alg.c 替换本目录的版本一起编译，检查移植是否改变了结果。
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common_util.h"
#include "securec.h"
#include "subevent_motion_alg.h"

#define MAX_STEP_NUM    (ALG_CHANNEL_NUM + SUBEVENT_MOTION_MAX_SUBEVENT_NUM)
#define LINE_LEN        16384

typedef struct {
    channel_select_t cfg;
    uint8_t hop[MAX_STEP_NUM];
    uint8_t per_subevent[SUBEVENT_MOTION_MAX_SUBEVENT_NUM];
    uint16_t time_per_channel[MAX_STEP_NUM];
    complex iq[SUBEVENT_MOTION_MAX_PATH_NUM][ALG_CHANNEL_NUM];
    uint8_t path_num;
} DumpProc;

// 解析一行中 "[jdps_dump]" 之后的内容，返回 1 表示一次测量结束
static int parse_line(char* line, DumpProc* d)
{
    char* p = strstr(line, "[jdps_dump]");
    if (p == NULL) {
        return 0;
    }
    char* tok = strtok(p + strlen("[jdps_dump]"), " \t\r\n");
    if (tok == NULL) {
        return 0;
    }
    if (strcmp(tok, "begin") == 0) {
        memset(d, 0, sizeof(*d));
    } else if (strcmp(tok, "cfg") == 0) {
        d->cfg.ch_num = (uint8_t)atoi(strtok(NULL, " "));
        d->cfg.subevent_num = (uint8_t)atoi(strtok(NULL, " "));
        d->cfg.t_mes = (uint16_t)atoi(strtok(NULL, " "));
    } else if (strcmp(tok, "hop") == 0) {
        for (int n = 0; (tok = strtok(NULL, " \r\n")) != NULL && n < MAX_STEP_NUM; n++) {
            d->hop[n] = (uint8_t)atoi(tok);
        }
    } else if (strcmp(tok, "per_subevent") == 0) {
        for (int n = 0; (tok = strtok(NULL, " \r\n")) != NULL && n < SUBEVENT_MOTION_MAX_SUBEVENT_NUM; n++) {
            d->per_subevent[n] = (uint8_t)atoi(tok);
        }
    } else if (strcmp(tok, "time") == 0) {
        for (int n = 0; (tok = strtok(NULL, " \r\n")) != NULL && n < MAX_STEP_NUM; n++) {
            d->time_per_channel[n] = (uint16_t)atoi(tok);
        }
    } else if (strcmp(tok, "path") == 0) {
        int path = atoi(strtok(NULL, " "));
        if (path < 0 || path >= SUBEVENT_MOTION_MAX_PATH_NUM) {
            return 0;
        }
        for (int c = 0; c < ALG_CHANNEL_NUM; c++) {
            char* re = strtok(NULL, " \r\n");
            char* im = strtok(NULL, " \r\n");
            if (re == NULL || im == NULL) {
                break;
            }
            d->iq[path][c].r = strtof(re, NULL);
            d->iq[path][c].i = strtof(im, NULL);
        }
        if (path + 1 > d->path_num) {
            d->path_num = (uint8_t)(path + 1);
        }
    } else if (strcmp(tok, "end") == 0) {
        return 1;
    }
    return 0;
}

static void run_proc(DumpProc* d, int index)
{
    static SubeventMotionCtx ctx;
    static complex iq_buf[ALG_CHANNEL_NUM];
    SubeventMotionResult res;

    d->cfg.ch_hop_orders = d->hop;
    d->cfg.ch_num_per_subevent = d->per_subevent;
    d->cfg.time_per_channel = d->time_per_channel;
    errcode_t err = subevent_motion_init(&ctx, &d->cfg, d->path_num);
    for (uint8_t p = 0; (err == ERRCODE_RANGING_ALG_SUCCESS) && (p < d->path_num); p++) {
        memcpy(iq_buf, d->iq[p], sizeof(iq_buf));
        err = subevent_motion_add_speed_path(&ctx, iq_buf);
    }
    if (err == ERRCODE_RANGING_ALG_SUCCESS) {
        err = subevent_motion_solve_speed(&ctx);
    }
    for (uint8_t p = 0; (err == ERRCODE_RANGING_ALG_SUCCESS) && (p < d->path_num); p++) {
        memcpy(iq_buf, d->iq[p], sizeof(iq_buf));
        err = subevent_motion_add_phase_path(&ctx, p, iq_buf);
    }
    if (err == ERRCODE_RANGING_ALG_SUCCESS) {
        err = subevent_motion_solve_phase(&ctx, &res);
    }
    printf("proc #%d: ", index);
    if (err != ERRCODE_RANGING_ALG_SUCCESS) {
        printf("error %u\n", (unsigned)err);
        return;
    }
    printf("[jdps] v=%+.3f m/s (est %+.3f, score %.1f, valid %u), psi[rad]=", (double)res.speed,
           (double)res.speed_est, (double)res.speed_score, (unsigned)res.speed_valid);
    for (int k = 0; k < d->cfg.subevent_num; k++) {
        printf("%s%+.3f", (k == 0) ? "" : " ", (double)res.subevent_phase[k]);
    }
    printf("\n");
}

int main(int argc, char** argv)
{
    static DumpProc d;
    static char line[LINE_LEN];
    int index = 0;
    FILE* fp = fopen(argc > 1 ? argv[1] : "firmware.log", "r");
    if (fp == NULL) {
        perror("open dump");
        return 2;
    }
    while (fgets(line, sizeof(line), fp) != NULL) {
        if (parse_line(line, &d)) {
            run_proc(&d, index++);
        }
    }
    fclose(fp);
    return 0;
}
