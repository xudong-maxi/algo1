/**
 * @file    jdps.h
 * @brief   JDPS (Joint Doppler & Phase Stitching): Doppler + AGC phase compensation
 *          for multi-carrier phase-based ranging measured in several segments.
 *
 * Protocol independent: applies to BLE Channel Sounding (segment = subevent),
 * SparkLink (SLE) ranging, or any scheme that measures a round-trip carrier phase
 * on a frequency-hopped channel set. The frequency plan is part of jdps_cfg_t.
 *
 * Terms
 *   step    : one frequency measurement (one channel, one time instant).
 *   segment : a run of consecutive steps measured with one AGC setting; every
 *             segment carries its own unknown constant phase (AGC on both sides).
 *   iq      : round-trip product (initiator tone x reflector tone) of one antenna
 *             path, one value per step in hop (time) order; a missing tone is 0+0j.
 *
 * Streaming interface: only ONE antenna path's iq has to be in memory at a time.
 * Everything shared between antennas is accumulated in jdps_ctx_t, so the iq of
 * each antenna is handed in once per pass:
 *
 *   jdps_begin(&ctx, &cfg, &layout);                          // antenna independent setup
 *   for (a = 0; a < A; a++) jdps_add_velocity(&ctx, iq_a);    // pass 1
 *   if (jdps_solve_velocity(&ctx) != JDPS_OK) -> skip (no signal)
 *   for (a = 0; a < A; a++) jdps_add_phase(&ctx, a, iq_a);    // pass 2
 *   jdps_solve_phase(&ctx, &res);
 *   for (a = 0; a < A; a++) jdps_apply(&ctx, a, iq_a);        // pass 3: in place,
 *                                                             //   then IFFT of path a
 * With cfg.per_ant_phase = 1, jdps_apply(a) may follow jdps_add_phase(a) directly
 * (passes 2 and 3 merge: two passes over the iq of every antenna).
 *
 * Algorithm (see docs/multi_subevent_agc_design.md):
 *   1. Pair every step with the step whose channel index is +1 / +2 higher, group
 *      the pairs by (index spacing, segment of lower tone, segment of higher tone).
 *   2. Velocity search: maximise sum over antennas/groups of
 *      | sum_{pairs in group} y_hi * conj(y_lo) * exp(j*alpha_p*v) |;
 *      reject (v = 0) if the noise-normalised peak score is low or on the edge.
 *   3. Compensate Doppler + range migration: y *= exp(j*4*pi*f*v*tau/c).
 *   4. Estimate one phase per segment from the group sums (power iteration on a
 *      K x K Hermitian matrix) and remove it.
 *
 * Properties: float32 only, no dynamic memory, one jdps_ctx_t per concurrent
 * procedure. One segment degenerates to a plain Doppler compensation.
 */
#ifndef JDPS_H
#define JDPS_H

#include <stdint.h>

/*
 * Complex type: taken from the host project, not defined here.
 * Required: a struct type named `complex` with float members `r` (real) and
 * `i` (imaginary), declared in that order (jdps.c initialises it as { r, i }).
 * Point JDPS_COMPLEX_HEADER at the project header that defines it
 * (default "complex_type.h"), e.g. -DJDPS_COMPLEX_HEADER='"my_types.h"'.
 */
#ifndef JDPS_COMPLEX_HEADER
#define JDPS_COMPLEX_HEADER "complex_type.h"
#endif
#include JDPS_COMPLEX_HEADER

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------------- */
/* Compile-time limits (size jdps_ctx_t)                                     */
/* ------------------------------------------------------------------------- */
/* Override on the compiler command line to save RAM, e.g. -DJDPS_MAX_SEG=4u.     */
#ifndef JDPS_MAX_ANT
#define JDPS_MAX_ANT        4u    /**< antenna paths                         */
#endif
#ifndef JDPS_MAX_STEPS
#define JDPS_MAX_STEPS      80u   /**< ranging steps per procedure           */
#endif
#ifndef JDPS_MAX_SEG
#define JDPS_MAX_SEG        8u    /**< segments per procedure                */
#endif
#ifndef JDPS_MAX_ORDER
#define JDPS_MAX_ORDER      2u    /**< largest channel-index spacing paired  */
#endif
#ifndef JDPS_MAX_V_POINTS
#define JDPS_MAX_V_POINTS   161u  /**< velocity grid points, incl. guard     */
#endif
#define JDPS_MAX_CHANNELS   80u   /**< channel index range 0..MAX_CHANNELS-1 */

#define JDPS_MAX_PAIRS      (JDPS_MAX_STEPS * JDPS_MAX_ORDER)
#define JDPS_MAX_GROUPS     (JDPS_MAX_ORDER * JDPS_MAX_SEG * JDPS_MAX_SEG)

/* ------------------------------------------------------------------------- */
/* Types                                                                     */
/* ------------------------------------------------------------------------- */
typedef enum {
    JDPS_OK = 0,
    JDPS_ERR_PARAM,          /**< NULL pointer, size out of range or wrong call order */
    JDPS_ERR_NO_PAIRS,       /**< no adjacent channel pair found                     */
    JDPS_ERR_NO_SIGNAL       /**< all pair products zero (e.g. remote IQ missing);
                                  skip passes 2/3, iq stays untouched                */
} jdps_status_t;

/** Algorithm configuration (use jdps_default_cfg() for recommended values). */
typedef struct {
    float   chan0_freq_hz;      /**< frequency of channel index 0 [Hz]            */
    float   chan_spacing_hz;    /**< channel spacing [Hz]                         */
    float   v_max_mps;          /**< supported velocity range: [-v_max, +v_max]   */
    float   v_guard_mps;        /**< extra search margin beyond v_max, so a true
                                     |v| close to v_max is not an edge hit       */
    float   v_step_mps;         /**< velocity grid step (0.25..1.0)              */
    float   min_v_score;        /**< velocity accepted if score >= this,
                                     otherwise v = 0 is applied                   */
    uint8_t num_orders;         /**< index spacings paired: 1..num_orders        */
    uint8_t power_iter_num;     /**< power iterations per phase-sync round       */
    uint8_t refine_iter_num;    /**< phase-sync refinement rounds                */
    uint8_t per_ant_phase;      /**< 1: segment phase estimated per antenna path */
} jdps_cfg_t;

/** Antenna-independent description of one ranging procedure.
 *  The arrays are referenced (not copied) and must stay valid until the last call. */
typedef struct {
    uint8_t        num_ant;     /**< antenna paths (<= MAX_ANT)                   */
    uint8_t        num_steps;   /**< ranging steps (<= MAX_STEPS)                 */
    uint8_t        num_seg;     /**< segments (<= MAX_SEG)                        */
    const uint8_t *chan_idx;    /**< [num_steps] channel index,
                                     f = chan0_freq_hz + idx * chan_spacing_hz    */
    const uint8_t *seg_idx;     /**< [num_steps] segment id 0..num_seg-1          */
    const float   *step_time_s; /**< [num_steps] step centre time, any origin [s] */
} jdps_layout_t;

typedef struct {
    float   v_mps;                               /**< velocity applied [m/s]:
                                                      v_est_mps if v_valid, else 0  */
    float   v_est_mps;                           /**< raw search result [m/s]       */
    float   v_score;                             /**< noise-normalised peak score
                                                      (peak - mu) / sigma; <~4.5 for
                                                      noise for any segment count   */
    uint8_t v_valid;                             /**< score >= min_v_score and the
                                                      peak is not on the grid edge  */
    float   psi_rad[JDPS_MAX_ANT][JDPS_MAX_SEG]; /**< removed segment phase [rad]
                                                      (AGC + inter-segment Doppler),
                                                      psi[.][0] = 0                 */
} jdps_result_t;

/**
 * Processing context: 11.1 KB with default limits, 7.1 KB with JDPS_MAX_SEG = 4,
 * 6.5 KB with JDPS_MAX_SEG = 3 (the iq itself is never stored). Place it statically.
 * Fields are internal; read results through jdps_result_t.
 */
typedef struct {
    /* setup */
    jdps_cfg_t    cfg;
    jdps_layout_t layout;
    float         t_ref;                          /**< reference time (mean step time)  */
    float         chan_ref;                       /**< mean channel index (numerics)    */
    uint16_t      num_pairs;
    uint16_t      num_groups;
    uint16_t      num_v_points;
    uint8_t       stage;                          /**< call-order state machine         */
    uint8_t       ant_count;                      /**< antennas added in pass 1         */
    uint8_t       ant_added;                      /**< bit a: antenna a added in pass 2 */
    uint8_t       ant_phase_ready;                /**< bit a: phases of antenna a known */

    /* antenna-independent pair table */
    uint8_t       pair_lo[JDPS_MAX_PAIRS];        /**< lower-frequency step             */
    uint8_t       pair_hi[JDPS_MAX_PAIRS];        /**< higher-frequency step            */
    uint8_t       pair_group[JDPS_MAX_PAIRS];     /**< group id                         */
    float         pair_rate[JDPS_MAX_PAIRS];      /**< alpha_p [rad/(m/s)]              */

    /* single-antenna scratch, reused for every antenna */
    complex       pair_prod[JDPS_MAX_PAIRS];      /**< y_hi * conj(y_lo)                */
    complex       pair_rot[JDPS_MAX_PAIRS];       /**< exp(j*alpha*v), running          */
    complex       pair_rot_step[JDPS_MAX_PAIRS];  /**< exp(j*alpha*dv)                  */
    complex       group_acc[JDPS_MAX_GROUPS];     /**< group sums at one grid point     */

    /* pass 1: accumulated over antennas */
    float         v_metric[JDPS_MAX_V_POINTS];    /**< velocity spectrum                */
    float         noise_mu;                       /**< metric mean for pure noise       */
    float         noise_var;                      /**< metric variance for pure noise   */
    float         v_est;
    float         v_score;
    float         v_used;
    uint8_t       v_valid;

    /* pass 2: per antenna group sums of the Doppler-compensated pairs */
    complex       group_sum[JDPS_MAX_ANT][JDPS_MAX_GROUPS];
    complex       seg_phasor[JDPS_MAX_ANT][JDPS_MAX_SEG]; /**< exp(j*psi)               */
} jdps_ctx_t;

/* ------------------------------------------------------------------------- */
/* API                                                                       */
/* ------------------------------------------------------------------------- */

/**
 * Recommended configuration: BLE CS frequency plan (2402 MHz + idx * 1 MHz),
 * +-10 m/s (+1 m/s guard) search, 0.25 m/s grid, index spacings 1 and 2,
 * velocity accepted for score >= 5.
 * For another protocol overwrite chan0_freq_hz / chan_spacing_hz.
 */
jdps_cfg_t jdps_default_cfg(void);

/** Validate the configuration / layout and build the antenna-independent pair table. */
jdps_status_t jdps_begin(jdps_ctx_t *ctx, const jdps_cfg_t *cfg, const jdps_layout_t *layout);

/** Pass 1: add one antenna path to the velocity spectrum (call once per antenna). */
jdps_status_t jdps_add_velocity(jdps_ctx_t *ctx, const complex *iq);

/** Estimate the velocity from all antennas; JDPS_ERR_NO_SIGNAL if every iq was zero. */
jdps_status_t jdps_solve_velocity(jdps_ctx_t *ctx);

/** Pass 2: add antenna path `ant` to the segment-phase estimation (iq not modified). */
jdps_status_t jdps_add_phase(jdps_ctx_t *ctx, uint8_t ant, const complex *iq);

/** Estimate the segment phases and fill the result (res may be NULL). */
jdps_status_t jdps_solve_phase(jdps_ctx_t *ctx, jdps_result_t *res);

/** Pass 3: compensate the iq of antenna path `ant` in place (Doppler + segment phase). */
jdps_status_t jdps_apply(const jdps_ctx_t *ctx, uint8_t ant, complex *iq);

#ifdef __cplusplus
}
#endif

#endif /* JDPS_H */
