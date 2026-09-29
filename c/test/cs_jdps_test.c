/**
 * @file    cs_jdps_test.c
 * @brief   Host test: runs cs_jdps_process() on vectors exported by
 *          cs_agc/export_vectors.py and compares with the Python reference.
 *
 * Usage: cs_jdps_test vectors.txt
 * Pass criteria per case:
 *   |v_c - v_py|                    < 0.01 m/s
 *   max wrapped |psi_c - psi_py|    < 0.01 rad
 *   ||iq_c - iq_py|| / ||iq_py||    < 1e-3
 */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "../cs_jdps.h"

#define TOL_V_MPS     0.01
#define TOL_PSI_RAD   0.01
#define TOL_IQ_REL    1e-3

typedef struct {
    char             name[64];
    cs_jdps_cfg_t    cfg;
    uint8_t          num_ant, num_steps, num_se;
    uint8_t          chan[CS_JDPS_MAX_STEPS];
    uint8_t          se[CS_JDPS_MAX_STEPS];
    float            time_s[CS_JDPS_MAX_STEPS];
    cs_cplx_t        iq[CS_JDPS_MAX_ANT][CS_JDPS_MAX_STEPS];
    double           exp_v;
    double           exp_psi[CS_JDPS_MAX_ANT][CS_JDPS_MAX_SE];
    double           exp_iq[CS_JDPS_MAX_ANT][CS_JDPS_MAX_STEPS][2];
} test_case_t;

static int expect_tag(FILE *fp, const char *tag)
{
    char buf[32];
    return fscanf(fp, "%31s", buf) == 1 && strcmp(buf, tag) == 0;
}

/** Reads one case; returns 0 at end of file. */
static int read_case(FILE *fp, test_case_t *tc)
{
    float vmax, vstep;
    unsigned orders, pit, rit, perant, A, N, K, u;

    if (fscanf(fp, " case %63s", tc->name) != 1) {
        return 0;
    }
    if (!expect_tag(fp, "cfg") ||
        fscanf(fp, "%f %f %u %u %u %u", &vmax, &vstep, &orders, &pit, &rit, &perant) != 6 ||
        !expect_tag(fp, "dims") || fscanf(fp, "%u %u %u", &A, &N, &K) != 3) {
        return -1;
    }
    tc->cfg.v_max_mps = vmax;
    tc->cfg.v_step_mps = vstep;
    tc->cfg.num_orders = (uint8_t)orders;
    tc->cfg.power_iter_num = (uint8_t)pit;
    tc->cfg.refine_iter_num = (uint8_t)rit;
    tc->cfg.per_ant_phase = (uint8_t)perant;
    tc->num_ant = (uint8_t)A;
    tc->num_steps = (uint8_t)N;
    tc->num_se = (uint8_t)K;

    if (!expect_tag(fp, "chan")) return -1;
    for (unsigned n = 0; n < N; n++) { if (fscanf(fp, "%u", &u) != 1) return -1; tc->chan[n] = (uint8_t)u; }
    if (!expect_tag(fp, "se")) return -1;
    for (unsigned n = 0; n < N; n++) { if (fscanf(fp, "%u", &u) != 1) return -1; tc->se[n] = (uint8_t)u; }
    if (!expect_tag(fp, "time")) return -1;
    for (unsigned n = 0; n < N; n++) { if (fscanf(fp, "%f", &tc->time_s[n]) != 1) return -1; }
    if (!expect_tag(fp, "iq")) return -1;
    for (unsigned a = 0; a < A; a++)
        for (unsigned n = 0; n < N; n++)
            if (fscanf(fp, "%f %f", &tc->iq[a][n].re, &tc->iq[a][n].im) != 2) return -1;
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
    static cs_jdps_work_t work;
    cs_jdps_result_t      res;
    int num_cases = 0, num_fail = 0;
    double total_us = 0.0;

    FILE *fp = fopen(argc > 1 ? argv[1] : "vectors.txt", "r");
    if (fp == NULL) {
        perror("open vectors");
        return 2;
    }
    printf("%-26s %9s %9s %10s %10s  %s\n", "case", "v_c", "v_py", "psi_err", "iq_relerr", "result");
    int rc;
    while ((rc = read_case(fp, &tc)) == 1) {
        cs_jdps_meas_t meas = {
            .num_ant = tc.num_ant, .num_steps = tc.num_steps, .num_se = tc.num_se,
            .chan_idx = tc.chan, .se_idx = tc.se, .step_time_s = tc.time_s, .iq = tc.iq,
        };
        clock_t t0 = clock();
        cs_jdps_status_t st = cs_jdps_process(&tc.cfg, &meas, &work, &res);
        total_us += 1e6 * (double)(clock() - t0) / CLOCKS_PER_SEC;

        double psi_err = 0.0, err2 = 0.0, ref2 = 0.0;
        for (unsigned a = 0; a < tc.num_ant; a++) {
            for (unsigned k = 0; k < tc.num_se; k++) {
                double e = fabs(wrap_pi(res.psi_rad[a][k] - tc.exp_psi[a][k]));
                psi_err = e > psi_err ? e : psi_err;
            }
            for (unsigned n = 0; n < tc.num_steps; n++) {
                double dr = tc.iq[a][n].re - tc.exp_iq[a][n][0];
                double di = tc.iq[a][n].im - tc.exp_iq[a][n][1];
                err2 += dr * dr + di * di;
                ref2 += tc.exp_iq[a][n][0] * tc.exp_iq[a][n][0] + tc.exp_iq[a][n][1] * tc.exp_iq[a][n][1];
            }
        }
        double iq_rel = sqrt(err2 / (ref2 > 0 ? ref2 : 1.0));
        int ok = st == CS_JDPS_OK && fabs(res.v_mps - tc.exp_v) < TOL_V_MPS &&
                 psi_err < TOL_PSI_RAD && iq_rel < TOL_IQ_REL;
        printf("%-26s %+9.4f %+9.4f %10.2e %10.2e  %s\n", tc.name, res.v_mps, tc.exp_v,
               psi_err, iq_rel, ok ? "PASS" : "FAIL");
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
