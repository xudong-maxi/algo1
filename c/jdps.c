/**
 * @file    jdps.c
 * @brief   JDPS: joint Doppler + AGC phase compensation for multi-segment
 *          carrier-phase ranging, streaming over antenna paths.
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

/* call-order state machine */
enum {
    STAGE_IDLE = 0,         /* jdps_begin not called / failed                */
    STAGE_VELOCITY,         /* pass 1: jdps_add_velocity                     */
    STAGE_PHASE,            /* pass 2: jdps_add_phase (v known)              */
    STAGE_DONE              /* jdps_solve_phase done: jdps_apply             */
};

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

/** Doppler + range-migration rotator of step n: exp(j*4*pi/c * f_n * v * (t_n - t_ref)) */
static inline complex doppler_rot(const jdps_ctx_t *ctx, uint8_t n)
{
    const jdps_layout_t *lay = &ctx->layout;
    return cplx_expj(JDPS_K_ROUND_TRIP * chan_freq_hz(&ctx->cfg, lay->chan_idx[n]) * ctx->v_used *
                     (lay->step_time_s[n] - ctx->t_ref));
}

/** Number of grid points on each side of v = 0 (search covers v_max + guard). */
static int16_t velocity_half_points(const jdps_cfg_t *cfg)
{
    return (int16_t)lrintf((cfg->v_max_mps + cfg->v_guard_mps) / cfg->v_step_mps);
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
    cfg.v_guard_mps     = 1.0f;
    cfg.min_v_score     = 5.0f;
    cfg.v_step_mps      = 0.25f;
    cfg.num_orders      = 2u;
    cfg.power_iter_num  = 8u;
    cfg.refine_iter_num = 3u;
    cfg.per_ant_phase   = 0u;
    return cfg;
}

/* ------------------------------------------------------------------------- */
/* Setup: time reference and adjacent-channel pairs (antenna independent)    */
/* ------------------------------------------------------------------------- */
static int check_setup(const jdps_cfg_t *cfg, const jdps_layout_t *lay)
{
    if (cfg == NULL || lay == NULL || lay->chan_idx == NULL || lay->seg_idx == NULL ||
        lay->step_time_s == NULL) {
        return 0;
    }
    if (lay->num_ant == 0u || lay->num_ant > JDPS_MAX_ANT ||
        lay->num_steps == 0u || lay->num_steps > JDPS_MAX_STEPS ||
        lay->num_seg == 0u || lay->num_seg > JDPS_MAX_SEG ||
        !(cfg->chan0_freq_hz > 0.0f) || !(cfg->chan_spacing_hz > 0.0f) ||
        cfg->num_orders == 0u || cfg->num_orders > JDPS_MAX_ORDER ||
        !(cfg->v_step_mps > 0.0f) || cfg->v_max_mps < 0.0f || cfg->v_guard_mps < 0.0f ||
        2L * velocity_half_points(cfg) + 1L > (long)JDPS_MAX_V_POINTS) {
        return 0;
    }
    for (uint8_t n = 0; n < lay->num_steps; n++) {
        if (lay->chan_idx[n] >= JDPS_MAX_CHANNELS || lay->seg_idx[n] >= lay->num_seg) {
            return 0;
        }
    }
    return 1;
}

/** t_ref = mean step time (the output distance refers to it), chan_ref = mean channel index. */
static void calc_references(jdps_ctx_t *ctx)
{
    const jdps_layout_t *lay = &ctx->layout;
    float t_sum = 0.0f;
    float ch_sum = 0.0f;
    for (uint8_t n = 0; n < lay->num_steps; n++) {
        t_sum  += lay->step_time_s[n];
        ch_sum += (float)lay->chan_idx[n];
    }
    ctx->t_ref    = t_sum / (float)lay->num_steps;
    ctx->chan_ref = ch_sum / (float)lay->num_steps;
}

/** Mean of (t - t_ref) per segment; removed from the search phase to keep it small. */
static void calc_seg_mean_tau(const jdps_ctx_t *ctx, float *seg_mean_tau)
{
    const jdps_layout_t *lay = &ctx->layout;
    uint8_t count[JDPS_MAX_SEG] = { 0 };
    for (uint8_t k = 0; k < lay->num_seg; k++) {
        seg_mean_tau[k] = 0.0f;
    }
    for (uint8_t n = 0; n < lay->num_steps; n++) {
        seg_mean_tau[lay->seg_idx[n]] += lay->step_time_s[n] - ctx->t_ref;
        count[lay->seg_idx[n]]++;
    }
    for (uint8_t k = 0; k < lay->num_seg; k++) {
        if (count[k] != 0u) {
            seg_mean_tau[k] /= (float)count[k];
        }
    }
}

/**
 * Pair (lo, hi) with chan[hi] = chan[lo] + order, order = 1..num_orders.
 * Group id = ((order-1)*K + seg[lo])*K + seg[hi].
 * Pair rate alpha = 4*pi/c * (f_hi*tau_hi - f_lo*tau_lo) minus its group-common
 * part 4*pi/c * f_ref * (mean_tau[seg_hi] - mean_tau[seg_lo]) (does not change
 * |group sum|, keeps float32 values small). f_ref = mean measured frequency.
 */
static void build_pairs(jdps_ctx_t *ctx)
{
    const jdps_layout_t *lay = &ctx->layout;
    const jdps_cfg_t    *cfg = &ctx->cfg;
    const uint8_t        num_seg = lay->num_seg;
    const float          freq_ref = cfg->chan0_freq_hz + ctx->chan_ref * cfg->chan_spacing_hz;
    int16_t step_of_chan[JDPS_MAX_CHANNELS];
    float   seg_mean_tau[JDPS_MAX_SEG];

    memset(step_of_chan, 0xFF, sizeof(step_of_chan));           /* -1 = unused */
    for (uint8_t n = 0; n < lay->num_steps; n++) {
        step_of_chan[lay->chan_idx[n]] = (int16_t)n;
    }
    calc_seg_mean_tau(ctx, seg_mean_tau);

    ctx->num_pairs = 0u;
    for (uint8_t order = 1u; order <= cfg->num_orders; order++) {
        for (uint8_t ch = 0u; ch + order < JDPS_MAX_CHANNELS; ch++) {
            int16_t lo = step_of_chan[ch];
            int16_t hi = step_of_chan[ch + order];
            if (lo < 0 || hi < 0) {
                continue;
            }
            uint16_t p      = ctx->num_pairs++;
            uint8_t  seg_lo = lay->seg_idx[lo];
            uint8_t  seg_hi = lay->seg_idx[hi];
            /* frequencies as offsets from f_ref: keeps f*tau products accurate in float32 */
            float df_lo    = ((float)lay->chan_idx[lo] - ctx->chan_ref) * cfg->chan_spacing_hz;
            float df_hi    = ((float)lay->chan_idx[hi] - ctx->chan_ref) * cfg->chan_spacing_hz;
            float tau_lo   = lay->step_time_s[lo] - ctx->t_ref;
            float tau_hi   = lay->step_time_s[hi] - ctx->t_ref;
            float local_dt = (tau_hi - seg_mean_tau[seg_hi]) - (tau_lo - seg_mean_tau[seg_lo]);

            ctx->pair_lo[p]    = (uint8_t)lo;
            ctx->pair_hi[p]    = (uint8_t)hi;
            ctx->pair_group[p] = (uint8_t)(((order - 1u) * num_seg + seg_lo) * num_seg + seg_hi);
            ctx->pair_rate[p]  = JDPS_K_ROUND_TRIP *
                                 (freq_ref * local_dt + df_hi * tau_hi - df_lo * tau_lo);
        }
    }
}

jdps_status_t jdps_begin(jdps_ctx_t *ctx, const jdps_cfg_t *cfg, const jdps_layout_t *layout)
{
    if (ctx == NULL) {
        return JDPS_ERR_PARAM;
    }
    memset(ctx, 0, sizeof(*ctx));                               /* stage = IDLE */
    if (!check_setup(cfg, layout)) {
        return JDPS_ERR_PARAM;
    }
    ctx->cfg          = *cfg;
    ctx->layout       = *layout;
    ctx->num_groups   = (uint16_t)(cfg->num_orders * layout->num_seg * layout->num_seg);
    ctx->num_v_points = (uint16_t)(2 * velocity_half_points(cfg) + 1);
    calc_references(ctx);
    build_pairs(ctx);
    if (ctx->num_pairs == 0u) {
        return JDPS_ERR_NO_PAIRS;
    }
    ctx->stage = STAGE_VELOCITY;
    return JDPS_OK;
}

/* ------------------------------------------------------------------------- */
/* Pass 1: velocity spectrum, accumulated antenna by antenna                 */
/* ------------------------------------------------------------------------- */
/** pair_prod[p] = y[hi] * conj(y[lo]) of the current antenna (amplitude-weighted phase difference). */
static void calc_pair_products(jdps_ctx_t *ctx, const complex *iq)
{
    for (uint16_t p = 0; p < ctx->num_pairs; p++) {
        ctx->pair_prod[p] = cplx_mul_conj(iq[ctx->pair_hi[p]], iq[ctx->pair_lo[p]]);
    }
}

/**
 * Noise statistics of the metric for this antenna. For random pair phases each
 * |group sum| is Rayleigh distributed with mean sqrt(pi/4 * S_g) and variance
 * (1 - pi/4) * S_g, S_g = sum |pair_prod|^2 over the group:
 *   mu += sum_g sqrt(pi/4 * S_g),  var += sum_g (1 - pi/4) * S_g.
 */
static void add_noise_stats(jdps_ctx_t *ctx)
{
    float group_power[JDPS_MAX_GROUPS];
    memset(group_power, 0, ctx->num_groups * sizeof(float));
    for (uint16_t p = 0; p < ctx->num_pairs; p++) {
        const complex z = ctx->pair_prod[p];
        group_power[ctx->pair_group[p]] += z.r * z.r + z.i * z.i;
    }
    for (uint16_t g = 0; g < ctx->num_groups; g++) {
        ctx->noise_mu  += sqrtf(0.25f * JDPS_PI * group_power[g]);
        ctx->noise_var += (1.0f - 0.25f * JDPS_PI) * group_power[g];
    }
}

/**
 * v_metric[iv] += sum_g | sum_{p in g} pair_prod[p] * exp(j*alpha_p*v_iv) |
 * on v = -(v_max + guard) .. +(v_max + guard); the rotator of every pair is
 * advanced by one complex multiply per grid point (no trigonometry in the loop).
 */
static void add_velocity_spectrum(jdps_ctx_t *ctx)
{
    const float v_step  = ctx->cfg.v_step_mps;
    const float v_first = -(float)velocity_half_points(&ctx->cfg) * v_step;

    for (uint16_t p = 0; p < ctx->num_pairs; p++) {
        ctx->pair_rot[p]      = cplx_expj(ctx->pair_rate[p] * v_first);
        ctx->pair_rot_step[p] = cplx_expj(ctx->pair_rate[p] * v_step);
    }
    for (uint16_t iv = 0; iv < ctx->num_v_points; iv++) {
        memset(ctx->group_acc, 0, ctx->num_groups * sizeof(complex));
        for (uint16_t p = 0; p < ctx->num_pairs; p++) {
            complex t = cplx_mul(ctx->pair_prod[p], ctx->pair_rot[p]);
            ctx->group_acc[ctx->pair_group[p]].r += t.r;
            ctx->group_acc[ctx->pair_group[p]].i += t.i;
        }
        float metric = 0.0f;
        for (uint16_t g = 0; g < ctx->num_groups; g++) {
            metric += cplx_abs(ctx->group_acc[g]);
        }
        ctx->v_metric[iv] += metric;
        /* advance every rotator to the next grid point */
        for (uint16_t p = 0; p < ctx->num_pairs; p++) {
            ctx->pair_rot[p] = cplx_mul(ctx->pair_rot[p], ctx->pair_rot_step[p]);
            if ((iv % JDPS_ROT_RENORM_EVERY) == (JDPS_ROT_RENORM_EVERY - 1u)) {
                ctx->pair_rot[p] = cplx_unit(ctx->pair_rot[p]);
            }
        }
    }
}

jdps_status_t jdps_add_velocity(jdps_ctx_t *ctx, const complex *iq)
{
    if (ctx == NULL || iq == NULL || ctx->stage != STAGE_VELOCITY ||
        ctx->ant_count >= ctx->layout.num_ant) {
        return JDPS_ERR_PARAM;
    }
    calc_pair_products(ctx, iq);
    add_noise_stats(ctx);
    add_velocity_spectrum(ctx);
    ctx->ant_count++;
    return JDPS_OK;
}

/**
 * Arg-max of the velocity spectrum with parabolic interpolation, its score
 * (peak - mu) / sigma, and the decision: v is used only if the score is high
 * enough and the peak is not on the first / last grid point (|v| may exceed
 * the range); otherwise v = 0 is applied.
 */
jdps_status_t jdps_solve_velocity(jdps_ctx_t *ctx)
{
    if (ctx == NULL || ctx->stage != STAGE_VELOCITY || ctx->ant_count != ctx->layout.num_ant) {
        return JDPS_ERR_PARAM;
    }
    if (!(ctx->noise_var > 0.0f)) {
        return JDPS_ERR_NO_SIGNAL;                              /* every iq was zero */
    }
    const uint16_t num_points = ctx->num_v_points;
    const float    v_step     = ctx->cfg.v_step_mps;
    uint16_t best = 0u;
    for (uint16_t iv = 1u; iv < num_points; iv++) {
        if (ctx->v_metric[iv] > ctx->v_metric[best]) {
            best = iv;
        }
    }
    float v_est = (float)((int16_t)best - velocity_half_points(&ctx->cfg)) * v_step;
    if (best > 0u && best + 1u < num_points) {
        float m_l = ctx->v_metric[best - 1u];
        float m_c = ctx->v_metric[best];
        float m_r = ctx->v_metric[best + 1u];
        float den = m_l - 2.0f * m_c + m_r;
        if (den < 0.0f) {
            v_est += 0.5f * (m_l - m_r) / den * v_step;
        }
    }
    const uint8_t at_edge = (uint8_t)(best == 0u || best + 1u == num_points);

    ctx->v_est   = v_est;
    ctx->v_score = (ctx->v_metric[best] - ctx->noise_mu) / sqrtf(ctx->noise_var);
    ctx->v_valid = (uint8_t)(ctx->v_score >= ctx->cfg.min_v_score && !at_edge);
    ctx->v_used  = ctx->v_valid ? v_est : 0.0f;
    ctx->stage   = STAGE_PHASE;
    return JDPS_OK;
}

/* ------------------------------------------------------------------------- */
/* Pass 2: segment phase synchronisation                                     */
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
static void sync_segment_phases(const jdps_ctx_t *ctx, uint8_t ant_first, uint8_t ant_num,
                                complex *seg_phasor)
{
    const jdps_cfg_t *cfg = &ctx->cfg;
    const uint8_t     K   = ctx->layout.num_seg;
    complex icpt_conj[JDPS_MAX_ANT][JDPS_MAX_ORDER];
    complex hmat[JDPS_MAX_SEG][JDPS_MAX_SEG];
    complex x[JDPS_MAX_SEG];
    complex x_next[JDPS_MAX_SEG];

#define GROUP_SUM(a, o, k, j) (ctx->group_sum[(a)][((o) * K + (k)) * K + (j)])

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

/**
 * group_sum[ant][g] = sum over pairs of group g of y'_hi * conj(y'_lo), where
 * y' = y * doppler_rot is the Doppler-compensated iq, computed on the fly
 * (no copy of the iq is made).
 */
jdps_status_t jdps_add_phase(jdps_ctx_t *ctx, uint8_t ant, const complex *iq)
{
    if (ctx == NULL || iq == NULL || ctx->stage != STAGE_PHASE || ant >= ctx->layout.num_ant ||
        (ctx->ant_added & (1u << ant)) != 0u) {
        return JDPS_ERR_PARAM;
    }
    complex *acc = ctx->group_sum[ant];
    memset(acc, 0, ctx->num_groups * sizeof(complex));
    for (uint16_t p = 0; p < ctx->num_pairs; p++) {
        const uint8_t lo = ctx->pair_lo[p];
        const uint8_t hi = ctx->pair_hi[p];
        complex y_lo = cplx_mul(iq[lo], doppler_rot(ctx, lo));
        complex y_hi = cplx_mul(iq[hi], doppler_rot(ctx, hi));
        complex t    = cplx_mul_conj(y_hi, y_lo);
        acc[ctx->pair_group[p]].r += t.r;
        acc[ctx->pair_group[p]].i += t.i;
    }
    ctx->ant_added |= (uint8_t)(1u << ant);

    if (ctx->cfg.per_ant_phase) {                   /* phases of this antenna known now */
        sync_segment_phases(ctx, ant, 1u, ctx->seg_phasor[ant]);
        ctx->ant_phase_ready |= (uint8_t)(1u << ant);
    }
    return JDPS_OK;
}

jdps_status_t jdps_solve_phase(jdps_ctx_t *ctx, jdps_result_t *res)
{
    if (ctx == NULL || ctx->stage != STAGE_PHASE) {
        return JDPS_ERR_PARAM;
    }
    const jdps_layout_t *lay      = &ctx->layout;
    const uint8_t        all_ants = (uint8_t)((1u << lay->num_ant) - 1u);
    if (ctx->ant_added != all_ants) {
        return JDPS_ERR_PARAM;
    }
    if (!ctx->cfg.per_ant_phase) {                  /* one phase set shared by all antennas */
        sync_segment_phases(ctx, 0u, lay->num_ant, ctx->seg_phasor[0]);
        for (uint8_t a = 1; a < lay->num_ant; a++) {
            memcpy(ctx->seg_phasor[a], ctx->seg_phasor[0], lay->num_seg * sizeof(complex));
        }
        ctx->ant_phase_ready = all_ants;
    }
    ctx->stage = STAGE_DONE;

    if (res != NULL) {
        memset(res, 0, sizeof(*res));
        res->v_mps     = ctx->v_used;
        res->v_est_mps = ctx->v_est;
        res->v_score   = ctx->v_score;
        res->v_valid   = ctx->v_valid;
        for (uint8_t a = 0; a < lay->num_ant; a++) {
            for (uint8_t k = 0; k < lay->num_seg; k++) {
                res->psi_rad[a][k] = atan2f(ctx->seg_phasor[a][k].i, ctx->seg_phasor[a][k].r);
            }
        }
    }
    return JDPS_OK;
}

/* ------------------------------------------------------------------------- */
/* Pass 3: compensation                                                      */
/* ------------------------------------------------------------------------- */
/** iq[n] *= exp(j*4*pi/c * f_n * v * (t_n - t_ref)) * exp(-j*psi[ant][seg_n]) */
jdps_status_t jdps_apply(const jdps_ctx_t *ctx, uint8_t ant, complex *iq)
{
    if (ctx == NULL || iq == NULL || ant >= ctx->layout.num_ant ||
        (ctx->ant_phase_ready & (1u << ant)) == 0u) {
        return JDPS_ERR_PARAM;
    }
    const jdps_layout_t *lay = &ctx->layout;
    for (uint8_t n = 0; n < lay->num_steps; n++) {
        complex rot = cplx_mul_conj(doppler_rot(ctx, n), ctx->seg_phasor[ant][lay->seg_idx[n]]);
        iq[n] = cplx_mul(iq[n], rot);
    }
    return JDPS_OK;
}
