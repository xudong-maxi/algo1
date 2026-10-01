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
 *   iq      : round-trip product (initiator tone x reflector tone) per step and
 *             antenna path, so that its phase is -4*pi*f*d/c + segment phase.
 *
 * Input : iq of one ranging procedure in hop (time) order, all antenna paths.
 * Output: the same iq, compensated in place so that all segments are phase
 *         continuous and referenced to one time instant (mean step time); it can
 *         be placed on the channel grid and fed directly into the IFFT ranging.
 *
 * Algorithm (see docs/multi_subevent_agc_design.md):
 *   1. Pair every step with the step whose channel index is +1 / +2 higher, group
 *      the pairs by (index spacing, segment of lower tone, segment of higher tone).
 *   2. Velocity search: maximise sum over antennas/groups of
 *      | sum_{pairs in group} y_hi * conj(y_lo) * exp(j*alpha_p*v) |.
 *   3. Compensate Doppler + range migration: y *= exp(j*4*pi*f*v*tau/c).
 *   4. Estimate one phase per segment from the group sums (power iteration on a
 *      K x K Hermitian matrix) and remove it.
 *
 * Properties: float32 only, no dynamic memory, not re-entrant per workspace
 * (use one jdps_work_t per concurrent call). One segment degenerates to a plain
 * Doppler compensation.
 */
#ifndef JDPS_H
#define JDPS_H

#include <stdint.h>

/*
 * Complex type: taken from the host project, not defined here.
 * Required: a struct type named `complex` with float members `r` (real) and
 * `i` (imaginary), declared in that order (jdps.c initialises it as { r, i }). Point JDPS_COMPLEX_HEADER at the project header that
 * defines it (default "complex_type.h"), e.g. -DJDPS_COMPLEX_HEADER='"my_types.h"'.
 */
#ifndef JDPS_COMPLEX_HEADER
#define JDPS_COMPLEX_HEADER "complex_type.h"
#endif
#include JDPS_COMPLEX_HEADER

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------------- */
/* Compile-time limits (size the static workspace)                           */
/* ------------------------------------------------------------------------- */
#define JDPS_MAX_ANT        4u    /**< antenna paths                        */
#define JDPS_MAX_STEPS      80u   /**< ranging steps per procedure          */
#define JDPS_MAX_SEG         8u    /**< segments per procedure              */
#define JDPS_MAX_ORDER      2u    /**< largest channel-index spacing paired */
#define JDPS_MAX_V_POINTS   161u  /**< velocity grid points, incl. guard    */
#define JDPS_MAX_CHANNELS   80u   /**< channel index range 0..MAX_CHANNELS-1 */

#define JDPS_MAX_PAIRS      (JDPS_MAX_STEPS * JDPS_MAX_ORDER)
#define JDPS_MAX_GROUPS     (JDPS_MAX_ORDER * JDPS_MAX_SEG * JDPS_MAX_SEG)

/* ------------------------------------------------------------------------- */
/* Types                                                                     */
/* ------------------------------------------------------------------------- */
typedef enum {
    JDPS_OK = 0,
    JDPS_ERR_PARAM,          /**< NULL pointer or size out of range       */
    JDPS_ERR_NO_PAIRS,       /**< no adjacent channel pair found          */
    JDPS_ERR_NO_SIGNAL       /**< all pair products zero (e.g. remote IQ
                                  missing); iq left untouched              */
} jdps_status_t;

/** Algorithm configuration (use jdps_default_cfg() for recommended values). */
typedef struct {
    float   chan0_freq_hz;      /**< frequency of channel index 0 [Hz]         */
    float   chan_spacing_hz;    /**< channel spacing [Hz]                      */
    float   v_max_mps;          /**< supported velocity range: [-v_max, +v_max] */
    float   v_guard_mps;        /**< extra search margin beyond v_max, so a true
                                     |v| close to v_max is not an edge hit      */
    float   v_step_mps;         /**< velocity grid step (0.25..1.0)           */
    float   min_v_score;        /**< velocity accepted if score >= this,
                                     otherwise v = 0 is applied (see result)  */
    uint8_t num_orders;         /**< index spacings paired: 1..num_orders     */
    uint8_t power_iter_num;     /**< power iterations per phase-sync round    */
    uint8_t refine_iter_num;    /**< phase-sync refinement rounds             */
    uint8_t per_ant_phase;      /**< 1: AGC phase estimated per antenna path  */
} jdps_cfg_t;

/** One ranging procedure. iq is compensated in place. */
typedef struct {
    uint8_t        num_ant;     /**< antenna paths used (<= MAX_ANT)          */
    uint8_t        num_steps;   /**< ranging steps (<= MAX_STEPS)             */
    uint8_t        num_seg;     /**< segments (<= MAX_SEG)                    */
    const uint8_t *chan_idx;    /**< [num_steps] channel index,
                                     f = chan0_freq_hz + idx * chan_spacing_hz */
    const uint8_t *seg_idx;     /**< [num_steps] segment id 0..num_seg-1      */
    const float   *step_time_s; /**< [num_steps] step centre time, any origin [s] */
    complex      (*iq)[JDPS_MAX_STEPS]; /**< [num_ant][..] round-trip IQ, hop order;
                                              a missing tone must be 0+0j   */
} jdps_meas_t;

typedef struct {
    float   v_mps;                             /**< velocity applied [m/s]:
                                                    v_est_mps if v_valid, else 0     */
    float   v_est_mps;                         /**< raw search result [m/s]          */
    float   v_score;                           /**< noise-normalised peak score:
                                                    (peak - mu) / sigma, mu/sigma of
                                                    the search metric for pure noise;
                                                    <~4.5 for noise for any segment
                                                    count                            */
    uint8_t v_valid;                           /**< 1: score >= min_v_score and the
                                                    peak is not on the grid edge     */
    float   psi_rad[JDPS_MAX_ANT][JDPS_MAX_SEG]; /**< removed segment phase [rad]
                                                    (AGC + inter-segment Doppler),
                                                    psi[.][0] = 0                  */
} jdps_result_t;

/** Scratch memory (13.9 KB with default limits). Place it statically. */
typedef struct {
    int16_t   pair_lo[JDPS_MAX_PAIRS];            /**< lower-frequency step   */
    int16_t   pair_hi[JDPS_MAX_PAIRS];            /**< higher-frequency step  */
    uint8_t   pair_group[JDPS_MAX_PAIRS];         /**< group id               */
    float     pair_rate[JDPS_MAX_PAIRS];          /**< alpha_i [rad/(m/s)]    */
    complex   pair_prod[JDPS_MAX_ANT][JDPS_MAX_PAIRS]; /**< y_hi*conj(y_lo) */
    complex   pair_rot[JDPS_MAX_PAIRS];           /**< exp(j*alpha*v) running */
    complex   pair_rot_step[JDPS_MAX_PAIRS];      /**< exp(j*alpha*dv)        */
    complex   group_sum[JDPS_MAX_ANT][JDPS_MAX_GROUPS];
    float     v_metric[JDPS_MAX_V_POINTS];
    uint16_t  num_pairs;
} jdps_work_t;

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

/**
 * Estimate velocity and segment phases, then compensate meas->iq in place.
 * @param cfg   configuration
 * @param meas  measurement; iq overwritten with the compensated IQ
 * @param work  scratch memory
 * @param res   estimated velocity and removed phases (may be NULL)
 */
jdps_status_t jdps_process(const jdps_cfg_t *cfg, jdps_meas_t *meas,
                                 jdps_work_t *work, jdps_result_t *res);

#ifdef __cplusplus
}
#endif

#endif /* JDPS_H */
