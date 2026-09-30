/**
 * @file    jdps.c
 * @brief   JDPS: joint Doppler + AGC phase compensation for multi-segment carrier-phase ranging.
 *          See jdps.h for the interface and docs/multi_subevent_agc_design.md
 *          for the derivation.
 */
#include "jdps.h"

#include <math.h>
#include <string.h>

/* ------------------------------------------------------------------------- */
/* Constants                                                                 */
/* ------------------------------------------------------------------------- */
#define JDPS_PI               3.14159265358979f
#define JDPS_TWO_PI           6.28318530717959f
#define JDPS_LIGHT_SPEED      299792458.0f
/* round-trip phase per (Hz * m): 4*pi/c                                    */
#define JDPS_K_ROUND_TRIP     (4.0f * JDPS_PI / JDPS_LIGHT_SPEED)
/* re-normalise recursive rotators every N grid points (float drift)       */
#define JDPS_ROT_RENORM_EVERY 16u

/* ------------------------------------------------------------------------- */
/* Complex helpers                                                           */
/* ------------------------------------------------------------------------- */
/* complex literals are written { r, i }: relies on the project type's member order */
static inline complex cplx_mul(complex a, complex b)
{
    return (complex){ a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r };
}

/** a * conj(b) */
static inline complex cplx_mul_conj(complex a, complex b)
{
    return (complex){ a.r * b.r + a.i * b.i, a.i * b.r - a.r * b.i };
}

static inline complex cplx_conj(complex a)
{
    return (complex){ a.r, -a.i };
}

static inline float cplx_abs(complex a)
{
    return sqrtf(a.r * a.r + a.i * a.i);
}

/** a / |a|; returns 1+0j for a == 0 */
static inline complex cplx_unit(complex a)
{
    float mag = cplx_abs(a);
    return (mag > 0.0f) ? (complex){ a.r / mag, a.i / mag } : (complex){ 1.0f, 0.0f };
}

/** exp(j*phase), phase wrapped to [-pi, pi] first for accuracy of sinf/cosf */
static inline complex cplx_expj(float phase)
{
    phase -= JDPS_TWO_PI * rintf(phase / JDPS_TWO_PI);
    return (complex){ cosf(phase), sinf(phase) };
}

static inline float chan_freq_hz(const jdps_cfg_t *cfg, uint8_t chan_idx)
{
    return cfg->chan0_freq_hz + (float)chan_idx * cfg->chan_spacing_hz;
}

/* ------------------------------------------------------------------------- */
/* Public: default configuration                                             */
/* ------------------------------------------------------------------------- */
jdps_cfg_t jdps_default_cfg(void)
{
    jdps_cfg_t cfg;
    cfg.chan0_freq_hz   = 2402.0e6f;       /* BLE CS channel 0 */
    cfg.chan_spacing_hz = 1.0e6f;
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
static float calc_ref_time(const jdps_meas_t *meas)
{
    float sum = 0.0f;
    for (uint8_t n = 0; n < meas->num_steps; n++) {
        sum += meas->step_time_s[n];
    }
    return sum / (float)meas->num_steps;
}

/** Mean of (t - t_ref) per segment; removed from the search phase to keep it small. */
static void calc_seg_mean_time(const jdps_meas_t *meas, float t_ref, float *seg_mean_tau)
{
    uint8_t count[JDPS_MAX_SEG] = { 0 };
    for (uint8_t k = 0; k < meas->num_seg; k++) {
        seg_mean_tau[k] = 0.0f;
    }
    for (uint8_t n = 0; n < meas->num_steps; n++) {
        seg_mean_tau[meas->seg_idx[n]] += meas->step_time_s[n] - t_ref;
        count[meas->seg_idx[n]]++;
    }
    for (uint8_t k = 0; k < meas->num_seg; k++) {
        if (count[k] != 0u) {
            seg_mean_tau[k] /= (float)count[k];
        }
    }
}

/* ------------------------------------------------------------------------- */
/* Step 1: adjacent-channel pairs                                            */
/* ------------------------------------------------------------------------- */
/**
 * Pair (lo, hi) with chan[hi] = chan[lo] + order, order = 1..num_orders.
 * Group id = ((order-1)*K + seg[lo])*K + seg[hi].
 * Pair rate alpha = 4*pi/c * (f_hi*tau_hi - f_lo*tau_lo) minus its group-common
 * part 4*pi/c * f_ref * (mean_tau[seg_hi] - mean_tau[seg_lo]) (does not change
 * |group sum|, keeps float32 values small). f_ref = mean measured frequency.
 */
static void build_pairs(const jdps_cfg_t *cfg, const jdps_meas_t *meas,
                        float t_ref, jdps_work_t *work)
{
    int16_t step_of_chan[JDPS_MAX_CHANNELS];
    float   seg_mean_tau[JDPS_MAX_SEG];
    const uint8_t num_seg = meas->num_seg;

    memset(step_of_chan, 0xFF, sizeof(step_of_chan));           /* -1 = unused */
    for (uint8_t n = 0; n < meas->num_steps; n++) {
        step_of_chan[meas->chan_idx[n]] = (int16_t)n;
    }
    calc_seg_mean_time(meas, t_ref, seg_mean_tau);

    /* reference frequency f_ref = frequency of the mean channel index */
    float chan_ref = 0.0f;
    for (uint8_t n = 0; n < meas->num_steps; n++) {
        chan_ref += (float)meas->chan_idx[n];
    }
    chan_ref /= (float)meas->num_steps;
    const float freq_ref = cfg->chan0_freq_hz + chan_ref * cfg->chan_spacing_hz;

    work->num_pairs = 0u;
    for (uint8_t order = 1u; order <= cfg->num_orders; order++) {
        for (uint8_t ch = 0u; ch + order < JDPS_MAX_CHANNELS; ch++) {
            int16_t lo = step_of_chan[ch];
            int16_t hi = step_of_chan[ch + order];
            if (lo < 0 || hi < 0) {
                continue;
            }
            uint16_t p     = work->num_pairs++;
            uint8_t  seg_lo = meas->seg_idx[lo];
            uint8_t  seg_hi = meas->seg_idx[hi];
            /* frequencies as offsets from f_ref: keeps f*tau products accurate in float32 */
            float df_lo  = ((float)meas->chan_idx[lo] - chan_ref) * cfg->chan_spacing_hz;
            float df_hi  = ((float)meas->chan_idx[hi] - chan_ref) * cfg->chan_spacing_hz;
            float tau_lo = meas->step_time_s[lo] - t_ref;
            float tau_hi = meas->step_time_s[hi] - t_ref;
            float local_dt = (tau_hi - seg_mean_tau[seg_hi]) - (tau_lo - seg_mean_tau[seg_lo]);

            work->pair_lo[p]    = lo;
            work->pair_hi[p]    = hi;
            work->pair_group[p] = (uint8_t)(((order - 1u) * num_seg + seg_lo) * num_seg + seg_hi);
            work->pair_rate[p]  = JDPS_K_ROUND_TRIP *
                                  (freq_ref * local_dt + df_hi * tau_hi - df_lo * tau_lo);
        }
    }
}

/** pair_prod[a][p] = y[a][hi] * conj(y[a][lo]) (amplitude-weighted phase difference). */
static void calc_pair_products(const jdps_meas_t *meas, jdps_work_t *work)
{
    for (uint8_t a = 0; a < meas->num_ant; a++) {
        for (uint16_t p = 0; p < work->num_pairs; p++) {
            work->pair_prod[a][p] = cplx_mul_conj(meas->iq[a][work->pair_hi[p]],
                                                  meas->iq[a][work->pair_lo[p]]);
        }
    }
}

/** group_sum[a][g] = sum of pair_prod over the pairs of group g. */
static void calc_group_sums(const jdps_meas_t *meas, jdps_work_t *work, uint16_t num_groups)
{
    for (uint8_t a = 0; a < meas->num_ant; a++) {
        memset(work->group_sum[a], 0, num_groups * sizeof(complex));
        for (uint16_t p = 0; p < work->num_pairs; p++) {
            complex *acc = &work->group_sum[a][work->pair_group[p]];
            acc->r += work->pair_prod[a][p].r;
            acc->i += work->pair_prod[a][p].i;
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
static float search_velocity(const jdps_cfg_t *cfg, const jdps_meas_t *meas,
                             jdps_work_t *work, uint16_t num_groups)
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
            complex *acc = work->group_sum[a];
            memset(acc, 0, num_groups * sizeof(complex));
            for (uint16_t p = 0; p < work->num_pairs; p++) {
                complex t = cplx_mul(work->pair_prod[a][p], work->pair_rot[p]);
                acc[work->pair_group[p]].r += t.r;
                acc[work->pair_group[p]].i += t.i;
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
static void compensate_doppler(const jdps_cfg_t *cfg, jdps_meas_t *meas, float v_mps, float t_ref)
{
    for (uint8_t n = 0; n < meas->num_steps; n++) {
        float phase = JDPS_K_ROUND_TRIP * chan_freq_hz(cfg, meas->chan_idx[n]) * v_mps *
                      (meas->step_time_s[n] - t_ref);
        complex rot = cplx_expj(phase);
        for (uint8_t a = 0; a < meas->num_ant; a++) {
            meas->iq[a][n] = cplx_mul(meas->iq[a][n], rot);
        }
    }
}

/* ------------------------------------------------------------------------- */
/* Step 4: segment phase synchronisation                                    */
/* ------------------------------------------------------------------------- */
/**
 * Model: group_sum[a][o][k][j] ~ |.| * exp(j*(c[a][o] + psi_j - psi_k)).
 * c   : intercept (local group delay) per antenna and order,
 * psi : phase of segment k (psi_0 = 0).
 * Phases are carried as unit phasors: seg_phasor[k] = exp(j*psi_k),
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
static void sync_segment_phases(const jdps_cfg_t *cfg, const jdps_work_t *work,
                                 uint8_t num_seg, uint8_t ant_first, uint8_t ant_num,
                                 complex *seg_phasor)
{
    const uint8_t K = num_seg;
    complex icpt_conj[JDPS_MAX_ANT][JDPS_MAX_ORDER];
    complex hmat[JDPS_MAX_SEG][JDPS_MAX_SEG];
    complex x[JDPS_MAX_SEG];
    complex x_next[JDPS_MAX_SEG];

#define GROUP_SUM(a, o, k, j) (work->group_sum[(a)][((o) * K + (k)) * K + (j)])

    for (uint8_t k = 0; k < K; k++) {
        seg_phasor[k] = (complex){ 1.0f, 0.0f };
    }
    /* initial intercept from intra-segment pairs only (independent of psi) */
    for (uint8_t a = ant_first; a < ant_first + ant_num; a++) {
        for (uint8_t o = 0; o < cfg->num_orders; o++) {
            complex diag = { 0.0f, 0.0f };
            for (uint8_t k = 0; k < K; k++) {
                diag.r += GROUP_SUM(a, o, k, k).r;
                diag.i += GROUP_SUM(a, o, k, k).i;
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
                complex q = { 0.0f, 0.0f };
                for (uint8_t a = ant_first; a < ant_first + ant_num; a++) {
                    for (uint8_t o = 0; o < cfg->num_orders; o++) {
                        complex t = cplx_mul(GROUP_SUM(a, o, k, j), icpt_conj[a][o]);
                        q.r += t.r;
                        q.i += t.i;
                    }
                }
                hmat[k][j] = q;
            }
        }
        /* H = Q + Q^H, then diagonal := max_k sum_j |H[k][j]| */
        float diag_shift = 0.0f;
        for (uint8_t k = 0; k < K; k++) {
            hmat[k][k] = (complex){ 2.0f * hmat[k][k].r, 0.0f };
            for (uint8_t j = k + 1u; j < K; j++) {
                complex h = { hmat[k][j].r + hmat[j][k].r, hmat[k][j].i - hmat[j][k].i };
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
            hmat[k][k] = (complex){ diag_shift, 0.0f };
        }
        /* power iteration, start from current estimate x = exp(-j*psi) */
        for (uint8_t k = 0; k < K; k++) {
            x[k] = cplx_conj(seg_phasor[k]);
        }
        for (uint8_t it = 0; it < cfg->power_iter_num; it++) {
            for (uint8_t k = 0; k < K; k++) {
                complex s = { 0.0f, 0.0f };
                for (uint8_t j = 0; j < K; j++) {
                    complex t = cplx_mul(hmat[k][j], x[j]);
                    s.r += t.r;
                    s.i += t.i;
                }
                x_next[k] = cplx_unit(s);
            }
            memcpy(x, x_next, K * sizeof(complex));
        }
        /* exp(j*psi_k) = conj(x_k), normalised so that psi_0 = 0 */
        for (uint8_t k = 0; k < K; k++) {
            seg_phasor[k] = cplx_mul(cplx_conj(x[k]), x[0]);
        }
        /* intercept update with all groups */
        for (uint8_t a = ant_first; a < ant_first + ant_num; a++) {
            for (uint8_t o = 0; o < cfg->num_orders; o++) {
                complex s = { 0.0f, 0.0f };
                for (uint8_t k = 0; k < K; k++) {
                    for (uint8_t j = 0; j < K; j++) {
                        /* * exp(-j(psi_j - psi_k)) = * conj(u_j) * u_k */
                        complex rot = cplx_mul_conj(seg_phasor[k], seg_phasor[j]);
                        complex t   = cplx_mul(GROUP_SUM(a, o, k, j), rot);
                        s.r += t.r;
                        s.i += t.i;
                    }
                }
                icpt_conj[a][o] = cplx_conj(cplx_unit(s));
            }
        }
    }
#undef GROUP_SUM
}

/** iq[a][n] *= exp(-j*psi[a][se_n]) */
static void remove_segment_phases(jdps_meas_t *meas,
                                   complex seg_phasor[JDPS_MAX_ANT][JDPS_MAX_SEG])
{
    for (uint8_t a = 0; a < meas->num_ant; a++) {
        for (uint8_t n = 0; n < meas->num_steps; n++) {
            meas->iq[a][n] = cplx_mul_conj(meas->iq[a][n], seg_phasor[a][meas->seg_idx[n]]);
        }
    }
}

/* ------------------------------------------------------------------------- */
/* Public: top level                                                         */
/* ------------------------------------------------------------------------- */
static int check_params(const jdps_cfg_t *cfg, const jdps_meas_t *meas,
                        const jdps_work_t *work)
{
    if (cfg == NULL || meas == NULL || work == NULL || meas->iq == NULL ||
        meas->chan_idx == NULL || meas->seg_idx == NULL || meas->step_time_s == NULL) {
        return 0;
    }
    if (meas->num_ant == 0u || meas->num_ant > JDPS_MAX_ANT ||
        meas->num_steps == 0u || meas->num_steps > JDPS_MAX_STEPS ||
        meas->num_seg == 0u || meas->num_seg > JDPS_MAX_SEG ||
        !(cfg->chan0_freq_hz > 0.0f) || !(cfg->chan_spacing_hz > 0.0f) ||
        cfg->num_orders == 0u || cfg->num_orders > JDPS_MAX_ORDER ||
        !(cfg->v_step_mps > 0.0f) || cfg->v_max_mps < 0.0f ||
        2L * lrintf(cfg->v_max_mps / cfg->v_step_mps) + 1L > (long)JDPS_MAX_V_POINTS) {
        return 0;
    }
    for (uint8_t n = 0; n < meas->num_steps; n++) {
        if (meas->chan_idx[n] >= JDPS_MAX_CHANNELS || meas->seg_idx[n] >= meas->num_seg) {
            return 0;
        }
    }
    return 1;
}

jdps_status_t jdps_process(const jdps_cfg_t *cfg, jdps_meas_t *meas,
                                 jdps_work_t *work, jdps_result_t *res)
{
    complex seg_phasor[JDPS_MAX_ANT][JDPS_MAX_SEG];

    if (!check_params(cfg, meas, work)) {
        return JDPS_ERR_PARAM;
    }
    const uint16_t num_groups = (uint16_t)(cfg->num_orders * meas->num_seg * meas->num_seg);
    const float    t_ref      = calc_ref_time(meas);

    /* 1. pairs + velocity */
    build_pairs(cfg, meas, t_ref, work);
    if (work->num_pairs == 0u) {
        return JDPS_ERR_NO_PAIRS;
    }
    calc_pair_products(meas, work);
    const float v_hat = search_velocity(cfg, meas, work, num_groups);

    /* 2. Doppler + range migration */
    compensate_doppler(cfg, meas, v_hat, t_ref);

    /* 3. segment phases from the compensated data */
    calc_pair_products(meas, work);
    calc_group_sums(meas, work, num_groups);
    if (cfg->per_ant_phase) {
        for (uint8_t a = 0; a < meas->num_ant; a++) {
            sync_segment_phases(cfg, work, meas->num_seg, a, 1u, seg_phasor[a]);
        }
    } else {
        sync_segment_phases(cfg, work, meas->num_seg, 0u, meas->num_ant, seg_phasor[0]);
        for (uint8_t a = 1; a < meas->num_ant; a++) {
            memcpy(seg_phasor[a], seg_phasor[0], meas->num_seg * sizeof(complex));
        }
    }
    remove_segment_phases(meas, seg_phasor);

    if (res != NULL) {
        memset(res, 0, sizeof(*res));
        res->v_mps = v_hat;
        for (uint8_t a = 0; a < meas->num_ant; a++) {
            for (uint8_t k = 0; k < meas->num_seg; k++) {
                res->psi_rad[a][k] = atan2f(seg_phasor[a][k].i, seg_phasor[a][k].r);
            }
        }
    }
    return JDPS_OK;
}
