#ifndef SUBEVENT_MOTION_ALG_H
#define SUBEVENT_MOTION_ALG_H

/*
 * 多 subevent 运动补偿 + AGC 相位对齐（JDPS 算法，原理见 docs/multi_subevent_agc_design.md）。
 *
 * 问题：每个 subevent 的 AGC 档位不同，各 subevent 带有一个未知的常数相位；
 *       同时目标运动产生多普勒相位，随机跳频下二者叠加，拼接后的相位不连续。
 * 做法：
 *   1. 信道 (f, f+o) 的 IQ 共轭乘积按 (间隔 o, f 所在 subevent, f+o 所在 subevent) 分组，
 *      o = 1..SUBEVENT_MOTION_PAIR_ORDER；
 *   2. 速度搜索：每组内相干累加、组间和路径间非相干累加，取峰值得到速度；
 *      峰值不可信（分数低或落在搜索边界）时速度按 0 处理；
 *   3. 用速度补偿多普勒和距离迁移；
 *   4. 由补偿后的分组和估计每个 subevent 的相位（K x K 幂迭代），并去除。
 *
 * 使用方式（每次只需要一路 IQ 在内存中，IQ 按信道号排列，可以多次读取）：
 *   subevent_motion_init(&ctx, channel_select_cfg, path_num);
 *   for (p = 0; p < path_num; p++) subevent_motion_add_speed_path(&ctx, iq_p);     // 第 1 遍
 *   subevent_motion_solve_speed(&ctx);
 *   for (p = 0; p < path_num; p++) subevent_motion_add_phase_path(&ctx, p, iq_p);  // 第 2 遍
 *   subevent_motion_solve_phase(&ctx, &res);
 *   for (p = 0; p < path_num; p++) subevent_motion_iq_compensation(&ctx, iq_p);    // 第 3 遍，之后做 IFFT
 *
 * 前提：AGC 按设备设定、每个 subevent 只设一次，所有天线路径共用同一组 subevent 相位；
 *       每个信道在一次测量中最多出现一次。
 *
 * channel_select_t 约定（以 72 个 mode-2 频点、3 个 subevent 为例）：
 *   - 每个 subevent 的第一步是 mode0；ch_num_per_subevent 含 mode0，例如 26, 26, 23；
 *   - time_per_channel 按时间顺序给出每一步的时长 (us)，含每个 subevent 的 mode0，共 75 项；
 *   - ch_hop_orders 只在第 0 项放一个 mode0，其后依次是 72 个 mode-2 信道，ch_num = 73；
 *   - t_mes 为上一个 subevent 结束到下一个 subevent（其 mode0）开始的间隔 (us)。
 *   不满足 ch_num = 1 + sum(ch_num_per_subevent[s] - 1) 时，subevent_motion_init 返回参数错误。
 */

/* 编译期上限，可在编译选项中覆盖 */
#ifndef SUBEVENT_MOTION_MAX_PATH_NUM
#define SUBEVENT_MOTION_MAX_PATH_NUM        4       // 天线路径数上限
#endif
#ifndef SUBEVENT_MOTION_MAX_SUBEVENT_NUM
#define SUBEVENT_MOTION_MAX_SUBEVENT_NUM    4       // subevent 个数上限
#endif
#ifndef SUBEVENT_MOTION_PAIR_ORDER
/* 配对的最大信道间隔：
 *   2（默认）同时使用间隔 1、2 的信道对，估计更稳；
 *   1 只用相邻信道，内存约 -0.5 KB、计算量约减半，低 SNR 时更容易判为速度不可信 */
#define SUBEVENT_MOTION_PAIR_ORDER          2
#endif
#define SUBEVENT_MOTION_SPEED_POINT_NUM     89      // 速度搜索点数：2 * (10 + 1) / 0.25 + 1
#define SUBEVENT_MOTION_GROUP_NUM           (SUBEVENT_MOTION_PAIR_ORDER * SUBEVENT_MOTION_MAX_SUBEVENT_NUM * \
                                             SUBEVENT_MOTION_MAX_SUBEVENT_NUM)

/**
 * @if Eng
 * @brief  Result of the multi-subevent motion and AGC phase correction.
 * @else
 * @brief  多 subevent 运动补偿 + AGC 相位对齐的结果。
 * @endif
 */
typedef struct {
    float speed;                                /* !< @if Eng Speed used for the compensation (m/s), 0 if not valid.
                                                      @else 实际用于补偿的速度 (m/s)，不可信时为 0。 @endif */
    float speed_est;                            /* !< @if Eng Raw speed estimate (m/s).
                                                      @else 速度搜索的原始估计值 (m/s)。 @endif */
    float speed_score;                          /* !< @if Eng Noise-normalised peak score of the speed spectrum.
                                                      @else 速度谱峰值的可信度分数，纯噪声时一般小于 4.5。 @endif */
    uint8_t speed_valid;                        /* !< @if Eng 1: speed estimate is trusted.
                                                      @else 1：速度估计可信。 @endif */
    float subevent_phase[SUBEVENT_MOTION_MAX_SUBEVENT_NUM]; /* !< @if Eng Phase removed from each subevent (rad).
                                                      @else 每个 subevent 被去除的相位 (rad)，第 0 个为 0。 @endif */
} SubeventMotionResult;

/**
 * @if Eng
 * @brief  Measurement time and subevent of every channel, built from channel_select_t.
 * @else
 * @brief  信道表：每个信道的测量时刻和所属 subevent，由 channel_select_t 生成。
 * @endif
 */
typedef struct {
    uint8_t subevent[ALG_CHANNEL_NUM];          /* !< @if Eng Subevent of each channel, 0xFF if not measured.
                                                      @else 每个信道所属的 subevent，未测量为 0xFF。 @endif */
    float tau[ALG_CHANNEL_NUM];                 /* !< @if Eng Measurement time of each channel minus reference time (s).
                                                      @else 每个信道的测量时刻减去参考时刻 (s)。 @endif */
    float subevent_mean_tau[SUBEVENT_MOTION_MAX_SUBEVENT_NUM]; /* !< @if Eng Mean tau of each subevent (s).
                                                      @else 每个 subevent 的平均时刻 (s)。 @endif */
    float mean_channel;                         /* !< @if Eng Mean index of the measured channels.
                                                      @else 已测信道的平均信道号，用于保持 float 精度。 @endif */
} SubeventChannelTable;

/**
 * @if Eng
 * @brief  Working memory of pass 1 (speed search).
 * @else
 * @brief  第 1 遍（速度搜索）的工作内存。
 * @endif
 */
typedef struct {
    float spectrum[SUBEVENT_MOTION_SPEED_POINT_NUM];    /* !< @if Eng Speed spectrum, summed over paths.
                                                      @else 速度谱，所有路径累加。 @endif */
    complex group_acc[SUBEVENT_MOTION_GROUP_NUM];       /* !< @if Eng Group sums at one search point.
                                                      @else 某个速度点上的分组和。 @endif */
    float noise_mean;                           /* !< @if Eng Mean of the speed spectrum for pure noise.
                                                      @else 纯噪声时速度谱的均值。 @endif */
    float noise_var;                            /* !< @if Eng Variance of the speed spectrum for pure noise.
                                                      @else 纯噪声时速度谱的方差。 @endif */
} SpeedSearchWork;

/**
 * @if Eng
 * @brief  Working memory of pass 2 (subevent phase estimation).
 * @else
 * @brief  第 2 遍（subevent 相位估计）的工作内存。
 * @endif
 */
typedef struct {
    complex group_sum[SUBEVENT_MOTION_MAX_PATH_NUM][SUBEVENT_MOTION_GROUP_NUM]; /* !< @if Eng Group sums of the
                                                            Doppler-compensated pairs, per path.
                                                      @else 每条路径多普勒补偿后的分组和。 @endif */
} PhaseEstimateWork;

/**
 * @if Eng
 * @brief  Working context, kept between the calls of one measurement. Allocate it statically.
 * @else
 * @brief  一次测量的处理上下文，在各次调用之间保存中间结果，建议静态分配。
 * @endif
 */
typedef struct {
    uint8_t stage;                              /* !< @if Eng Call-order state.
                                                      @else 调用顺序状态。 @endif */
    uint8_t path_num;                           /* !< @if Eng The number of antenna paths.
                                                      @else 天线路径数。 @endif */
    uint8_t subevent_num;                       /* !< @if Eng The number of subevents.
                                                      @else subevent 个数。 @endif */
    uint8_t path_count;                         /* !< @if Eng Paths added in pass 1.
                                                      @else 第 1 遍已加入的路径数。 @endif */
    uint8_t path_added;                         /* !< @if Eng Bit p: path p added in pass 2.
                                                      @else 第 2 遍已加入路径的位图。 @endif */
    SubeventChannelTable channel;               /* !< @if Eng Channel table.
                                                      @else 信道表。 @endif */
    SubeventMotionResult result;                /* !< @if Eng Speed and subevent phase results.
                                                      @else 速度和 subevent 相位结果。 @endif */
    complex subevent_phasor[SUBEVENT_MOTION_MAX_SUBEVENT_NUM]; /* !< @if Eng exp(j * subevent phase).
                                                      @else 每个 subevent 相位的单位复数，补偿时使用。 @endif */
    union {                                     /* 第 1 遍和第 2 遍不会同时进行，共用内存 */
        SpeedSearchWork speed;
        PhaseEstimateWork phase;
    } work;
} SubeventMotionCtx;

errcode_t subevent_motion_init(SubeventMotionCtx* ctx, channel_select_t* channel_select_cfg, uint8_t path_num);
errcode_t subevent_motion_add_speed_path(SubeventMotionCtx* ctx, complex* iq);
errcode_t subevent_motion_solve_speed(SubeventMotionCtx* ctx);
errcode_t subevent_motion_add_phase_path(SubeventMotionCtx* ctx, uint8_t path, complex* iq);
errcode_t subevent_motion_solve_phase(SubeventMotionCtx* ctx, SubeventMotionResult* res);
errcode_t subevent_motion_iq_compensation(SubeventMotionCtx* ctx, complex* iq);

#endif  /* SUBEVENT_MOTION_ALG_H */
