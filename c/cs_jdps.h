/**
 * @file    cs_jdps.h
 * @brief   BLE Channel Sounding multi-subevent Doppler + AGC phase compensation
 *          (JDPS: Joint Doppler & Phase Stitching).
 *
 * Input : PBR round-trip IQ (initiator tone x reflector tone) of one CS procedure,
 *         in hop (time) order, for every antenna path, whose mode-2 steps are split
 *         over one or more subevents. Each subevent carries an unknown constant
 *         phase (AGC gear on both sides) and the target may move (Doppler).
 * Output: the same IQ, compensated in place so that all subevents are phase
 *         continuous and referenced to one time instant (mean step time); it can
 *         be put on the 1 MHz grid and fed directly into the IFFT ranging.
 *
 * Algorithm (see docs/multi_subevent_agc_design.md):
 *   1. Pair every step with the step whose channel is +1 / +2 MHz higher, group
 *      the pairs by (spacing, subevent of lower tone, subevent of higher tone).
 *   2. Velocity search: maximise sum over antennas/groups of
 *      | sum_{pairs in group} y_m * conj(y_n) * exp(j*alpha_i*v) |.
 *   3. Compensate Doppler + range migration: y *= exp(j*4*pi*f*v*tau/c).
 *   4. Estimate one phase per subevent from the group sums (power iteration on a
 *      K x K Hermitian matrix) and remove it.
 *
 * Properties: float32 only, no dynamic memory, not re-entrant per workspace
 * (use one cs_jdps_work_t per concurrent call). K = 1 degenerates to a
 * single-subevent Doppler compensation.
 */
#ifndef CS_JDPS_H
#define CS_JDPS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------------- */
/* Compile-time limits (size the static workspace)                           */
/* ------------------------------------------------------------------------- */
#define CS_JDPS_MAX_ANT        4u    /**< antenna paths                        */
#define CS_JDPS_MAX_STEPS      80u   /**< mode-2 steps per procedure           */
#define CS_JDPS_MAX_SE         8u    /**< subevents per procedure              */
#define CS_JDPS_MAX_ORDER      2u    /**< largest channel spacing used [MHz]   */
#define CS_JDPS_MAX_V_POINTS   161u  /**< velocity grid points                 */
#define CS_JDPS_NUM_CHANNELS   80u   /**< CS channel index range 0..79         */

#define CS_JDPS_MAX_PAIRS      (CS_JDPS_MAX_STEPS * CS_JDPS_MAX_ORDER)
#define CS_JDPS_MAX_GROUPS     (CS_JDPS_MAX_ORDER * CS_JDPS_MAX_SE * CS_JDPS_MAX_SE)

/* ------------------------------------------------------------------------- */
/* Types                                                                     */
/* ------------------------------------------------------------------------- */
typedef struct {
    float re;
    float im;
} cs_cplx_t;

typedef enum {
    CS_JDPS_OK = 0,
    CS_JDPS_ERR_PARAM,          /**< NULL pointer or size out of range       */
    CS_JDPS_ERR_NO_PAIRS        /**< no adjacent channel pair found          */
} cs_jdps_status_t;

/** Algorithm configuration (use cs_jdps_default_cfg() for recommended values). */
typedef struct {
    float   v_max_mps;          /**< velocity search range: [-v_max, +v_max]  */
    float   v_step_mps;         /**< velocity grid step (0.25..1.0)           */
    uint8_t num_orders;         /**< channel spacings used: 1..num_orders MHz */
    uint8_t power_iter_num;     /**< power iterations per phase-sync round    */
    uint8_t refine_iter_num;    /**< phase-sync refinement rounds             */
    uint8_t per_ant_phase;      /**< 1: AGC phase estimated per antenna path  */
} cs_jdps_cfg_t;

/** One CS procedure. iq is compensated in place. */
typedef struct {
    uint8_t        num_ant;     /**< antenna paths used (<= MAX_ANT)          */
    uint8_t        num_steps;   /**< mode-2 steps (<= MAX_STEPS)              */
    uint8_t        num_se;      /**< subevents (<= MAX_SE)                    */
    const uint8_t *chan_idx;    /**< [num_steps] CS channel index (f = 2402 + idx MHz) */
    const uint8_t *se_idx;      /**< [num_steps] subevent id 0..num_se-1      */
    const float   *step_time_s; /**< [num_steps] step centre time, any origin [s] */
    cs_cplx_t    (*iq)[CS_JDPS_MAX_STEPS]; /**< [num_ant][..] PBR IQ, hop order;
                                                a missing tone must be 0+0j   */
} cs_jdps_meas_t;

typedef struct {
    float v_mps;                                    /**< estimated radial velocity */
    float psi_rad[CS_JDPS_MAX_ANT][CS_JDPS_MAX_SE]; /**< removed subevent phase
                                                         (AGC + inter-subevent Doppler),
                                                         psi[.][0] = 0             */
} cs_jdps_result_t;

/** Scratch memory (13.9 KB with default limits). Place it statically. */
typedef struct {
    int16_t   pair_lo[CS_JDPS_MAX_PAIRS];            /**< lower-frequency step   */
    int16_t   pair_hi[CS_JDPS_MAX_PAIRS];            /**< higher-frequency step  */
    uint8_t   pair_group[CS_JDPS_MAX_PAIRS];         /**< group id               */
    float     pair_rate[CS_JDPS_MAX_PAIRS];          /**< alpha_i [rad/(m/s)]    */
    cs_cplx_t pair_prod[CS_JDPS_MAX_ANT][CS_JDPS_MAX_PAIRS]; /**< y_hi*conj(y_lo) */
    cs_cplx_t pair_rot[CS_JDPS_MAX_PAIRS];           /**< exp(j*alpha*v) running */
    cs_cplx_t pair_rot_step[CS_JDPS_MAX_PAIRS];      /**< exp(j*alpha*dv)        */
    cs_cplx_t group_sum[CS_JDPS_MAX_ANT][CS_JDPS_MAX_GROUPS];
    float     v_metric[CS_JDPS_MAX_V_POINTS];
    uint16_t  num_pairs;
} cs_jdps_work_t;

/* ------------------------------------------------------------------------- */
/* API                                                                       */
/* ------------------------------------------------------------------------- */

/** Recommended configuration: +-10 m/s, 0.25 m/s grid, spacings 1 and 2 MHz. */
cs_jdps_cfg_t cs_jdps_default_cfg(void);

/**
 * Estimate velocity and subevent phases, then compensate meas->iq in place.
 * @param cfg   configuration
 * @param meas  measurement; iq overwritten with the compensated IQ
 * @param work  scratch memory
 * @param res   estimated velocity and removed phases (may be NULL)
 */
cs_jdps_status_t cs_jdps_process(const cs_jdps_cfg_t *cfg, cs_jdps_meas_t *meas,
                                 cs_jdps_work_t *work, cs_jdps_result_t *res);

#ifdef __cplusplus
}
#endif

#endif /* CS_JDPS_H */
