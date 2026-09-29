/**
 * @file    cs_jdps.c
 * @brief   JDPS: joint Doppler + AGC phase compensation for multi-subevent CS.
 *          See cs_jdps.h for the interface and docs/multi_subevent_agc_design.md
 *          for the derivation.
 */
#include "cs_jdps.h"

#include <math.h>
#include <string.h>

/* ------------------------------------------------------------------------- */
/* Constants                                                                 */
/* ------------------------------------------------------------------------- */
#define JDPS_PI               3.14159265358979f
#define JDPS_TWO_PI           6.28318530717959f
#define JDPS_LIGHT_SPEED      299792458.0f
#define JDPS_CH0_FREQ_HZ      2402.0e6f          /* CS channel 0              */
#define JDPS_CH_SPACING_HZ    1.0e6f
#define JDPS_CENTER_FREQ_HZ   2440.0e6f          /* only used for numerics    */
/* round-trip phase per (Hz * m): 4*pi/c                                    */
#define JDPS_K_ROUND_TRIP     (4.0f * JDPS_PI / JDPS_LIGHT_SPEED)
/* re-normalise recursive rotators every N grid points (float drift)       */
#define JDPS_ROT_RENORM_EVERY 16u

/* ------------------------------------------------------------------------- */
/* Complex helpers                                                           */
/* ------------------------------------------------------------------------- */
static inline cs_cplx_t cplx(float re, float im)
{
    cs_cplx_t z = { re, im };
    return z;
}

static inline cs_cplx_t cplx_mul(cs_cplx_t a, cs_cplx_t b)
{
    return cplx(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}

/** a * conj(b) */
static inline cs_cplx_t cplx_mul_conj(cs_cplx_t a, cs_cplx_t b)
{
    return cplx(a.re * b.re + a.im * b.im, a.im * b.re - a.re * b.im);
}

static inline cs_cplx_t cplx_conj(cs_cplx_t a)
{
    return cplx(a.re, -a.im);
}

static inline float cplx_abs(cs_cplx_t a)
{
    return sqrtf(a.re * a.re + a.im * a.im);
}

/** a / |a|; returns 1+0j for a == 0 */
static inline cs_cplx_t cplx_unit(cs_cplx_t a)
{
    float mag = cplx_abs(a);
    return (mag > 0.0f) ? cplx(a.re / mag, a.im / mag) : cplx(1.0f, 0.0f);
}

/** exp(j*phase), phase wrapped to [-pi, pi] first for accuracy of sinf/cosf */
static inline cs_cplx_t cplx_expj(float phase)
{
    phase -= JDPS_TWO_PI * rintf(phase / JDPS_TWO_PI);
    return cplx(cosf(phase), sinf(phase));
}

static inline float chan_freq_hz(uint8_t chan_idx)
{
    return JDPS_CH0_FREQ_HZ + (float)chan_idx * JDPS_CH_SPACING_HZ;
}

/* ------------------------------------------------------------------------- */
/* Public: default configuration                                             */
/* ------------------------------------------------------------------------- */
cs_jdps_cfg_t cs_jdps_default_cfg(void)
{
    cs_jdps_cfg_t cfg;
    cfg.v_max_mps       = 10.0f;
    cfg.v_step_mps      = 0.25f;
    cfg.num_orders      = 2u;
    cfg.power_iter_num  = 8u;
    cfg.refine_iter_num = 3u;
    cfg.per_ant_phase   = 0u;
    return cfg;
}

/* ------------------------------------------------------------------------- */
/* Step 0: time reference                                                    */
/* ------------------------------------------------------------------------- */
/** Reference time t_ref = mean step time; the output distance refers to it. */
static float calc_ref_time(const cs_jdps_meas_t *meas)
{
    float sum = 0.0f;
    for (uint8_t n = 0; n < meas->num_steps; n++) {
        sum += meas->step_time_s[n];
    }
    return sum / (float)meas->num_steps;
}

/** Mean of (t - t_ref) per subevent; removed from the search phase to keep it small. */
static void calc_se_mean_time(const cs_jdps_meas_t *meas, float t_ref, float *se_mean_tau)
{
    uint8_t count[CS_JDPS_MAX_SE] = { 0 };
    for (uint8_t k = 0; k < meas->num_se; k++) {
        se_mean_tau[k] = 0.0f;
    }
    for (uint8_t n = 0; n < meas->num_steps; n++) {
        se_mean_tau[meas->se_idx[n]] += meas->step_time_s[n] - t_ref;
        count[meas->se_idx[n]]++;
    }
    for (uint8_t k = 0; k < meas->num_se; k++) {
        if (count[k] != 0u) {
            se_mean_tau[k] /= (float)count[k];
        }
    }
}

/* ------------------------------------------------------------------------- */
/* Step 1: adjacent-channel pairs                                            */
/* ------------------------------------------------------------------------- */
/**
 * Pair (lo, hi) with chan[hi] = chan[lo] + order, order = 1..num_orders.
 * Group id = ((order-1)*K + se[lo])*K + se[hi].
 * Pair rate alpha = 4*pi/c * (f_hi*tau_hi - f_lo*tau_lo) minus its group-common
 * part 4*pi/c * f_c * (mean_tau[se_hi] - mean_tau[se_lo]) (does not change |group sum|,
 * keeps float32 values small).
 */
static void build_pairs(const cs_jdps_cfg_t *cfg, const cs_jdps_meas_t *meas,
                        float t_ref, cs_jdps_work_t *work)
{
    int16_t step_of_chan[CS_JDPS_NUM_CHANNELS];
    float   se_mean_tau[CS_JDPS_MAX_SE];
    const uint8_t num_se = meas->num_se;

    memset(step_of_chan, 0xFF, sizeof(step_of_chan));           /* -1 = unused */
    for (uint8_t n = 0; n < meas->num_steps; n++) {
        step_of_chan[meas->chan_idx[n]] = (int16_t)n;
    }
    calc_se_mean_time(meas, t_ref, se_mean_tau);

    work->num_pairs = 0u;
    for (uint8_t order = 1u; order <= cfg->num_orders; order++) {
        for (uint8_t ch = 0u; ch + order < CS_JDPS_NUM_CHANNELS; ch++) {
            int16_t lo = step_of_chan[ch];
            int16_t hi = step_of_chan[ch + order];
            if (lo < 0 || hi < 0) {
                continue;
            }
            uint16_t p     = work->num_pairs++;
            uint8_t  se_lo = meas->se_idx[lo];
            uint8_t  se_hi = meas->se_idx[hi];
            /* offset frequencies from f_c: keeps f*tau products in float32 range */
            float df_lo  = chan_freq_hz(meas->chan_idx[lo]) - JDPS_CENTER_FREQ_HZ;
            float df_hi  = chan_freq_hz(meas->chan_idx[hi]) - JDPS_CENTER_FREQ_HZ;
            float tau_lo = meas->step_time_s[lo] - t_ref;
            float tau_hi = meas->step_time_s[hi] - t_ref;
            float local_dt = (tau_hi - se_mean_tau[se_hi]) - (tau_lo - se_mean_tau[se_lo]);

            work->pair_lo[p]    = lo;
            work->pair_hi[p]    = hi;
            work->pair_group[p] = (uint8_t)(((order - 1u) * num_se + se_lo) * num_se + se_hi);
            work->pair_rate[p]  = JDPS_K_ROUND_TRIP *
                                  (JDPS_CENTER_FREQ_HZ * local_dt + df_hi * tau_hi - df_lo * tau_lo);
        }
    }
}

/** pair_prod[a][p] = y[a][hi] * conj(y[a][lo]) (amplitude-weighted phase difference). */
static void calc_pair_products(const cs_jdps_meas_t *meas, cs_jdps_work_t *work)
{
    for (uint8_t a = 0; a < meas->num_ant; a++) {
        for (uint16_t p = 0; p < work->num_pairs; p++) {
            work->pair_prod[a][p] = cplx_mul_conj(meas->iq[a][work->pair_hi[p]],
                                                  meas->iq[a][work->pair_lo[p]]);
        }
    }
}

/** group_sum[a][g] = sum of pair_prod over the pairs of group g. */
static void calc_group_sums(const cs_jdps_meas_t *meas, cs_jdps_work_t *work, uint16_t num_groups)
{
    for (uint8_t a = 0; a < meas->num_ant; a++) {
        memset(work->group_sum[a], 0, num_groups * sizeof(cs_cplx_t));
        for (uint16_t p = 0; p < work->num_pairs; p++) {
            cs_cplx_t *acc = &work->group_sum[a][work->pair_group[p]];
            acc->re += work->pair_prod[a][p].re;
            acc->im += work->pair_prod[a][p].im;
        }
    }
}

/* ------------------------------------------------------------------------- */
/* Step 2: velocity search                                                   */
/* ------------------------------------------------------------------------- */
/**
 * metric(v) = sum_a sum_g | sum_{p in g} pair_prod[a][p] * exp(j*alpha_p*v) |
 * evaluated on v = -v_max .. +v_max; the rotator of every pair is advanced by
 * one complex multiply per grid point (no trigonometry in the inner loop).
 * Returns the parabolically interpolated arg-max.
 */
static float search_velocity(const cs_jdps_cfg_t *cfg, const cs_jdps_meas_t *meas,
                             cs_jdps_work_t *work, uint16_t num_groups)
{
    const int16_t  half_points = (int16_t)lrintf(cfg->v_max_mps / cfg->v_step_mps);
    const uint16_t num_points  = (uint16_t)(2 * half_points + 1);
    const float    v_first     = -(float)half_points * cfg->v_step_mps;

    for (uint16_t p = 0; p < work->num_pairs; p++) {
        work->pair_rot[p]      = cplx_expj(work->pair_rate[p] * v_first);
        work->pair_rot_step[p] = cplx_expj(work->pair_rate[p] * cfg->v_step_mps);
    }

    uint16_t best = 0u;
    for (uint16_t iv = 0; iv < num_points; iv++) {
        float metric = 0.0f;
        for (uint8_t a = 0; a < meas->num_ant; a++) {
            cs_cplx_t *acc = work->group_sum[a];
            memset(acc, 0, num_groups * sizeof(cs_cplx_t));
            for (uint16_t p = 0; p < work->num_pairs; p++) {
                cs_cplx_t t = cplx_mul(work->pair_prod[a][p], work->pair_rot[p]);
                acc[work->pair_group[p]].re += t.re;
                acc[work->pair_group[p]].im += t.im;
            }
            for (uint16_t g = 0; g < num_groups; g++) {
                metric += cplx_abs(acc[g]);
            }
        }
        work->v_metric[iv] = metric;
        if (metric > work->v_metric[best]) {
            best = iv;
        }
        /* advance every rotator to the next grid point */
        for (uint16_t p = 0; p < work->num_pairs; p++) {
            work->pair_rot[p] = cplx_mul(work->pair_rot[p], work->pair_rot_step[p]);
            if ((iv % JDPS_ROT_RENORM_EVERY) == (JDPS_ROT_RENORM_EVERY - 1u)) {
                work->pair_rot[p] = cplx_unit(work->pair_rot[p]);
            }
        }
    }

    /* parabolic interpolation around the peak */
    float v_hat = v_first + (float)best * cfg->v_step_mps;
    if (best > 0u && best + 1u < num_points) {
        float m_l = work->v_metric[best - 1u];
        float m_c = work->v_metric[best];
        float m_r = work->v_metric[best + 1u];
        float den = m_l - 2.0f * m_c + m_r;
        if (den < 0.0f) {
            v_hat += 0.5f * (m_l - m_r) / den * cfg->v_step_mps;
        }
    }
    return v_hat;
}

/* ------------------------------------------------------------------------- */
/* Step 3: Doppler + range-migration compensation                            */
/* ------------------------------------------------------------------------- */
/** iq[a][n] *= exp(j*4*pi/c * f_n * v * (t_n - t_ref)) */
static void compensate_doppler(cs_jdps_meas_t *meas, float v_mps, float t_ref)
{
    for (uint8_t n = 0; n < meas->num_steps; n++) {
        float phase = JDPS_K_ROUND_TRIP * chan_freq_hz(meas->chan_idx[n]) * v_mps *
                      (meas->step_time_s[n] - t_ref);
        cs_cplx_t rot = cplx_expj(phase);
        for (uint8_t a = 0; a < meas->num_ant; a++) {
            meas->iq[a][n] = cplx_mul(meas->iq[a][n], rot);
        }
    }
}

/* ------------------------------------------------------------------------- */
/* Step 4: subevent phase synchronisation                                    */
/* ------------------------------------------------------------------------- */
/**
 * Model: group_sum[a][o][k][j] ~ |.| * exp(j*(c[a][o] + psi_j - psi_k)).
 * c   : intercept (local group delay) per antenna and order,
 * psi : phase of subevent k (psi_0 = 0).
 * Phases are carried as unit phasors: se_phasor[k] = exp(j*psi_k),
 * icpt_conj[a][o] = exp(-j*c[a][o]).
 *
 * Each refinement round:
 *   Q[k][j]  = sum_{a,o} group_sum * icpt_conj                  (~ exp(j(psi_j-psi_k)))
 *   H        = Q + Q^H, diagonal set to max row sum             (shift: K = 2 would
 *                                                                oscillate otherwise)
 *   x        = principal eigenvector of H by power iteration    (~ exp(-j*psi))
 *   c        = arg sum_{k,j} group_sum * exp(-j(psi_j - psi_k)) (intercept update)
 *
 * @param ant_first, ant_num  antennas that share the estimated phases
 */
static void sync_subevent_phases(const cs_jdps_cfg_t *cfg, const cs_jdps_work_t *work,
                                 uint8_t num_se, uint8_t ant_first, uint8_t ant_num,
                                 cs_cplx_t *se_phasor)
{
    const uint8_t K = num_se;
    cs_cplx_t icpt_conj[CS_JDPS_MAX_ANT][CS_JDPS_MAX_ORDER];
    cs_cplx_t hmat[CS_JDPS_MAX_SE][CS_JDPS_MAX_SE];
    cs_cplx_t x[CS_JDPS_MAX_SE];
    cs_cplx_t x_next[CS_JDPS_MAX_SE];

#define GROUP_SUM(a, o, k, j) (work->group_sum[(a)][((o) * K + (k)) * K + (j)])

    for (uint8_t k = 0; k < K; k++) {
        se_phasor[k] = cplx(1.0f, 0.0f);
    }
    /* initial intercept from intra-subevent pairs only (independent of psi) */
    for (uint8_t a = ant_first; a < ant_first + ant_num; a++) {
        for (uint8_t o = 0; o < cfg->num_orders; o++) {
            cs_cplx_t diag = cplx(0.0f, 0.0f);
            for (uint8_t k = 0; k < K; k++) {
                diag.re += GROUP_SUM(a, o, k, k).re;
                diag.im += GROUP_SUM(a, o, k, k).im;
            }
            icpt_conj[a][o] = cplx_conj(cplx_unit(diag));
        }
    }
    if (K == 1u) {
        return;
    }

    for (uint8_t round = 0; round < cfg->refine_iter_num; round++) {
        /* Q[k][j] (stored in hmat) */
        for (uint8_t k = 0; k < K; k++) {
            for (uint8_t j = 0; j < K; j++) {
                cs_cplx_t q = cplx(0.0f, 0.0f);
                for (uint8_t a = ant_first; a < ant_first + ant_num; a++) {
                    for (uint8_t o = 0; o < cfg->num_orders; o++) {
                        cs_cplx_t t = cplx_mul(GROUP_SUM(a, o, k, j), icpt_conj[a][o]);
                        q.re += t.re;
                        q.im += t.im;
                    }
                }
                hmat[k][j] = q;
            }
        }
        /* H = Q + Q^H, then diagonal := max_k sum_j |H[k][j]| */
        float diag_shift = 0.0f;
        for (uint8_t k = 0; k < K; k++) {
            hmat[k][k] = cplx(2.0f * hmat[k][k].re, 0.0f);
            for (uint8_t j = k + 1u; j < K; j++) {
                cs_cplx_t h = cplx(hmat[k][j].re + hmat[j][k].re, hmat[k][j].im - hmat[j][k].im);
                hmat[k][j] = h;
                hmat[j][k] = cplx_conj(h);
            }
        }
        for (uint8_t k = 0; k < K; k++) {
            float row_sum = 0.0f;
            for (uint8_t j = 0; j < K; j++) {
                row_sum += cplx_abs(hmat[k][j]);
            }
            diag_shift = (row_sum > diag_shift) ? row_sum : diag_shift;
        }
        for (uint8_t k = 0; k < K; k++) {
            hmat[k][k] = cplx(diag_shift, 0.0f);
        }
        /* power iteration, start from current estimate x = exp(-j*psi) */
        for (uint8_t k = 0; k < K; k++) {
            x[k] = cplx_conj(se_phasor[k]);
        }
        for (uint8_t it = 0; it < cfg->power_iter_num; it++) {
            for (uint8_t k = 0; k < K; k++) {
                cs_cplx_t s = cplx(0.0f, 0.0f);
                for (uint8_t j = 0; j < K; j++) {
                    cs_cplx_t t = cplx_mul(hmat[k][j], x[j]);
                    s.re += t.re;
                    s.im += t.im;
                }
                x_next[k] = cplx_unit(s);
            }
            memcpy(x, x_next, K * sizeof(cs_cplx_t));
        }
        /* exp(j*psi_k) = conj(x_k), normalised so that psi_0 = 0 */
        for (uint8_t k = 0; k < K; k++) {
            se_phasor[k] = cplx_mul(cplx_conj(x[k]), x[0]);
        }
        /* intercept update with all groups */
        for (uint8_t a = ant_first; a < ant_first + ant_num; a++) {
            for (uint8_t o = 0; o < cfg->num_orders; o++) {
                cs_cplx_t s = cplx(0.0f, 0.0f);
                for (uint8_t k = 0; k < K; k++) {
                    for (uint8_t j = 0; j < K; j++) {
                        /* * exp(-j(psi_j - psi_k)) = * conj(u_j) * u_k */
                        cs_cplx_t rot = cplx_mul_conj(se_phasor[k], se_phasor[j]);
                        cs_cplx_t t   = cplx_mul(GROUP_SUM(a, o, k, j), rot);
                        s.re += t.re;
                        s.im += t.im;
                    }
                }
                icpt_conj[a][o] = cplx_conj(cplx_unit(s));
            }
        }
    }
#undef GROUP_SUM
}

/** iq[a][n] *= exp(-j*psi[a][se_n]) */
static void remove_subevent_phases(cs_jdps_meas_t *meas,
                                   cs_cplx_t se_phasor[CS_JDPS_MAX_ANT][CS_JDPS_MAX_SE])
{
    for (uint8_t a = 0; a < meas->num_ant; a++) {
        for (uint8_t n = 0; n < meas->num_steps; n++) {
            meas->iq[a][n] = cplx_mul_conj(meas->iq[a][n], se_phasor[a][meas->se_idx[n]]);
        }
    }
}

/* ------------------------------------------------------------------------- */
/* Public: top level                                                         */
/* ------------------------------------------------------------------------- */
static int check_params(const cs_jdps_cfg_t *cfg, const cs_jdps_meas_t *meas,
                        const cs_jdps_work_t *work)
{
    if (cfg == NULL || meas == NULL || work == NULL || meas->iq == NULL ||
        meas->chan_idx == NULL || meas->se_idx == NULL || meas->step_time_s == NULL) {
        return 0;
    }
    if (meas->num_ant == 0u || meas->num_ant > CS_JDPS_MAX_ANT ||
        meas->num_steps == 0u || meas->num_steps > CS_JDPS_MAX_STEPS ||
        meas->num_se == 0u || meas->num_se > CS_JDPS_MAX_SE ||
        cfg->num_orders == 0u || cfg->num_orders > CS_JDPS_MAX_ORDER ||
        !(cfg->v_step_mps > 0.0f) || cfg->v_max_mps < 0.0f ||
        2L * lrintf(cfg->v_max_mps / cfg->v_step_mps) + 1L > (long)CS_JDPS_MAX_V_POINTS) {
        return 0;
    }
    for (uint8_t n = 0; n < meas->num_steps; n++) {
        if (meas->chan_idx[n] >= CS_JDPS_NUM_CHANNELS || meas->se_idx[n] >= meas->num_se) {
            return 0;
        }
    }
    return 1;
}

cs_jdps_status_t cs_jdps_process(const cs_jdps_cfg_t *cfg, cs_jdps_meas_t *meas,
                                 cs_jdps_work_t *work, cs_jdps_result_t *res)
{
    cs_cplx_t se_phasor[CS_JDPS_MAX_ANT][CS_JDPS_MAX_SE];

    if (!check_params(cfg, meas, work)) {
        return CS_JDPS_ERR_PARAM;
    }
    const uint16_t num_groups = (uint16_t)(cfg->num_orders * meas->num_se * meas->num_se);
    const float    t_ref      = calc_ref_time(meas);

    /* 1. pairs + velocity */
    build_pairs(cfg, meas, t_ref, work);
    if (work->num_pairs == 0u) {
        return CS_JDPS_ERR_NO_PAIRS;
    }
    calc_pair_products(meas, work);
    const float v_hat = search_velocity(cfg, meas, work, num_groups);

    /* 2. Doppler + range migration */
    compensate_doppler(meas, v_hat, t_ref);

    /* 3. subevent phases from the compensated data */
    calc_pair_products(meas, work);
    calc_group_sums(meas, work, num_groups);
    if (cfg->per_ant_phase) {
        for (uint8_t a = 0; a < meas->num_ant; a++) {
            sync_subevent_phases(cfg, work, meas->num_se, a, 1u, se_phasor[a]);
        }
    } else {
        sync_subevent_phases(cfg, work, meas->num_se, 0u, meas->num_ant, se_phasor[0]);
        for (uint8_t a = 1; a < meas->num_ant; a++) {
            memcpy(se_phasor[a], se_phasor[0], meas->num_se * sizeof(cs_cplx_t));
        }
    }
    remove_subevent_phases(meas, se_phasor);

    if (res != NULL) {
        memset(res, 0, sizeof(*res));
        res->v_mps = v_hat;
        for (uint8_t a = 0; a < meas->num_ant; a++) {
            for (uint8_t k = 0; k < meas->num_se; k++) {
                res->psi_rad[a][k] = atan2f(se_phasor[a][k].im, se_phasor[a][k].re);
            }
        }
    }
    return CS_JDPS_OK;
}
