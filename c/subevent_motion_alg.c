#include <math.h>
#include "common_util.h"
#include "securec.h"
#include "subevent_motion_alg.h"

// 频点规划：信道号 ch 的频率为 CHANNEL0_FREQ + ch * CHANNEL_SPACING
#define CHANNEL0_FREQ       2402.0e6f       // 0 号信道频率 (Hz)
#define CHANNEL_SPACING     1.0e6f          // 信道间隔 (Hz)

// 速度搜索配置
#define MOTION_MAX_SPEED            10.0f   // 支持的最大速度 (m/s)
#define MOTION_SPEED_GUARD          1.0f    // 搜索范围在最大速度外再扩展的保护带 (m/s)
#define MOTION_SPEED_RESOLUTION     0.25f   // 速度搜索步长 (m/s)
#define MOTION_MIN_SPEED_SCORE      5.0f    // 速度谱峰值分数低于该值时认为速度不可信

// subevent 相位估计配置
#define PHASE_POWER_ITER_NUM        8       // 每轮幂迭代次数
#define PHASE_REFINE_NUM            3       // 截距与相位交替估计的轮数

#define SUBEVENT_NOT_MEASURED       0xFF    // 信道未被测量
#define ROUND_TRIP_FACTOR           (4.0f * (float)PI / (float)LIGHT_SPEED)    // 往返相位系数 4 * pi / c
#define SUBEVENT_TWO_PI             (2.0f * (float)PI)

// 调用顺序
#define STAGE_IDLE      0
#define STAGE_SPEED     1   // 第 1 遍：加入路径估计速度
#define STAGE_PHASE     2   // 第 2 遍：加入路径估计 subevent 相位
#define STAGE_DONE      3   // 可以补偿 IQ

// ------------------------------------------------------------------------------------------------
// 复数运算
// ------------------------------------------------------------------------------------------------
static inline complex complex_mul(complex a, complex b)
{
    complex z = { a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r };
    return z;
}

// a * conj(b)
static inline complex complex_mul_conj(complex a, complex b)
{
    complex z = { a.r * b.r + a.i * b.i, a.i * b.r - a.r * b.i };
    return z;
}

static inline complex complex_conj(complex a)
{
    complex z = { a.r, -a.i };
    return z;
}

static inline float complex_abs(complex a)
{
    return sqrtf(a.r * a.r + a.i * a.i);
}

// a / |a|，a 为 0 时返回 1
static inline complex complex_unit(complex a)
{
    float amp = complex_abs(a);
    complex one = { 1.0f, 0.0f };
    if (amp <= 0.0f) {
        return one;
    }
    complex z = { a.r / amp, a.i / amp };
    return z;
}

// exp(j * phase)，先把相位限制在 [-pi, pi] 以保证 sinf/cosf 精度
static inline complex complex_expj(float phase)
{
    phase = remainderf(phase, SUBEVENT_TWO_PI);
    complex z = { cosf(phase), sinf(phase) };
    return z;
}

static inline void complex_acc(complex* acc, complex a)
{
    acc->r += a.r;
    acc->i += a.i;
}

// ------------------------------------------------------------------------------------------------
// 信道表：每个信道的测量时刻和所属 subevent
// ------------------------------------------------------------------------------------------------
static inline float channel_freq(uint8_t ch)
{
    return CHANNEL0_FREQ + (float)ch * CHANNEL_SPACING;
}

static inline bool channel_measured(SubeventMotionCtx* ctx, uint8_t ch)
{
    return ctx->channel_subevent[ch] != SUBEVENT_NOT_MEASURED;
}

// 信道对 (ch, ch + order) 是否都被测量
static inline bool pair_valid(SubeventMotionCtx* ctx, uint8_t ch, uint8_t order)
{
    return (ch + order < ALG_CHANNEL_NUM) && channel_measured(ctx, ch) && channel_measured(ctx, ch + order);
}

// 信道对的分组号：(间隔 order, ch 所在 subevent, ch + order 所在 subevent)
static inline uint8_t pair_group(SubeventMotionCtx* ctx, uint8_t ch, uint8_t order)
{
    uint8_t num = ctx->subevent_num;
    return (uint8_t)(((order - 1) * num + ctx->channel_subevent[ch]) * num + ctx->channel_subevent[ch + order]);
}

static inline uint8_t group_count(SubeventMotionCtx* ctx)
{
    return (uint8_t)(SUBEVENT_MOTION_PAIR_ORDER * ctx->subevent_num * ctx->subevent_num);
}

static errcode_t build_channel_table(SubeventMotionCtx* ctx, channel_select_t* channel_select_cfg)
{
    uint8_t length = channel_select_cfg->ch_num;
    uint32_t current_time_us = 0;
    uint8_t current_step = 0;
    uint16_t channel_count = 0;
    float sum_time = 0.0f;
    float sum_channel = 0.0f;
    uint8_t subevent_count[SUBEVENT_MOTION_MAX_SUBEVENT_NUM] = {0};

    memset_s(ctx->channel_subevent, sizeof(ctx->channel_subevent), SUBEVENT_NOT_MEASURED,
             sizeof(ctx->channel_subevent));

    // 与 motion_effect_analyze 相同的时间约定：subevent 内逐步累加每个信道的测量时长，
    // subevent 之间加上间隔 t_mes；ch_hop_orders[0] 是 mode0，第 step 步的信道是 ch_hop_orders[step + 1]
    for (uint8_t s = 0; s < ctx->subevent_num; ++s) {
        for (uint8_t i = 0; i < channel_select_cfg->ch_num_per_subevent[s]; ++i) {
            if (current_step + 1 >= length) {
                break;
            }
            uint8_t ch = channel_select_cfg->ch_hop_orders[current_step + 1];
            // 只记录第一次被选到的信道
            if ((ch < ALG_CHANNEL_NUM) && !channel_measured(ctx, ch)) {
                float t_sec = (float)current_time_us * 1e-6f;
                ctx->channel_subevent[ch] = s;
                ctx->channel_tau[ch] = t_sec;
                sum_time += t_sec;
                sum_channel += (float)ch;
                subevent_count[s]++;
                channel_count++;
            }
            current_time_us += channel_select_cfg->time_per_channel[current_step];
            current_step++;
        }
        current_time_us += channel_select_cfg->t_mes;
    }
    if (channel_count < TWO) {
        return ERRCODE_RANGING_ALG_NOT_ENOUGH_IQ;
    }

    // 时间以平均测量时刻为参考，输出的距离对应该时刻；信道号以平均值为参考，保持 float 精度
    float t_ref = sum_time / (float)channel_count;
    ctx->channel_ref = sum_channel / (float)channel_count;
    memset_s(ctx->subevent_mean_tau, sizeof(ctx->subevent_mean_tau), 0, sizeof(ctx->subevent_mean_tau));
    for (uint8_t ch = 0; ch < ALG_CHANNEL_NUM; ++ch) {
        if (channel_measured(ctx, ch)) {
            ctx->channel_tau[ch] -= t_ref;
            ctx->subevent_mean_tau[ctx->channel_subevent[ch]] += ctx->channel_tau[ch];
        }
    }
    for (uint8_t s = 0; s < ctx->subevent_num; ++s) {
        if (subevent_count[s] != 0) {
            ctx->subevent_mean_tau[s] /= (float)subevent_count[s];
        }
    }
    return ERRCODE_RANGING_ALG_SUCCESS;
}

// 信道对 (ch, ch + order) 的速度相位率 (rad / (m/s))：
// 4 * pi / c * (f_hi * tau_hi - f_lo * tau_lo)，减去组内公共部分 4 * pi / c * f_ref * (所在 subevent 平均时刻之差)。
// 减去的部分不改变分组和的模，只是让数值变小，保证 float 精度。
static float pair_speed_rate(SubeventMotionCtx* ctx, uint8_t ch, uint8_t order)
{
    uint8_t hi = ch + order;
    float freq_ref = CHANNEL0_FREQ + ctx->channel_ref * CHANNEL_SPACING;
    float df_lo = ((float)ch - ctx->channel_ref) * CHANNEL_SPACING;
    float df_hi = ((float)hi - ctx->channel_ref) * CHANNEL_SPACING;
    float tau_lo = ctx->channel_tau[ch];
    float tau_hi = ctx->channel_tau[hi];
    float local_dt = (tau_hi - ctx->subevent_mean_tau[ctx->channel_subevent[hi]]) -
                     (tau_lo - ctx->subevent_mean_tau[ctx->channel_subevent[ch]]);
    return ROUND_TRIP_FACTOR * (freq_ref * local_dt + df_hi * tau_hi - df_lo * tau_lo);
}

// 多普勒和距离迁移补偿量 exp(j * 4 * pi / c * f * v * tau)
static inline complex doppler_rotation(SubeventMotionCtx* ctx, uint8_t ch)
{
    return complex_expj(ROUND_TRIP_FACTOR * channel_freq(ch) * ctx->speed * ctx->channel_tau[ch]);
}

errcode_t subevent_motion_init(SubeventMotionCtx* ctx, channel_select_t* channel_select_cfg, uint8_t path_num)
{
    if ((ctx == NULL) || (channel_select_cfg == NULL) || (channel_select_cfg->ch_hop_orders == NULL) ||
        (channel_select_cfg->time_per_channel == NULL) || (channel_select_cfg->ch_num_per_subevent == NULL)) {
        return ERRCODE_RANGING_ALG_INVALID_PARAM;
    }
    memset_s(ctx, sizeof(SubeventMotionCtx), 0, sizeof(SubeventMotionCtx));
    if ((path_num == 0) || (path_num > SUBEVENT_MOTION_MAX_PATH_NUM) || (channel_select_cfg->subevent_num == 0) ||
        (channel_select_cfg->subevent_num > SUBEVENT_MOTION_MAX_SUBEVENT_NUM)) {
        return ERRCODE_RANGING_ALG_INVALID_PARAM;
    }
    ctx->path_num = path_num;
    ctx->subevent_num = channel_select_cfg->subevent_num;

    errcode_t err = build_channel_table(ctx, channel_select_cfg);
    if (err != ERRCODE_RANGING_ALG_SUCCESS) {
        return err;
    }
    ctx->stage = STAGE_SPEED;
    return ERRCODE_RANGING_ALG_SUCCESS;
}

// ------------------------------------------------------------------------------------------------
// 第 1 遍：速度搜索
// ------------------------------------------------------------------------------------------------
// 纯噪声时，每个分组和的模服从瑞利分布：均值 sqrt(pi / 4 * S)，方差 (1 - pi / 4) * S，S 为组内 |z|^2 之和。
// 累加得到速度谱在纯噪声下的均值和方差，用于判断速度峰值是否可信（与 subevent 个数无关）。
static void add_noise_statistics(SubeventMotionCtx* ctx, complex* iq)
{
    float group_power[SUBEVENT_MOTION_GROUP_NUM] = {0};

    for (uint8_t order = 1; order <= SUBEVENT_MOTION_PAIR_ORDER; ++order) {
        for (uint8_t ch = 0; ch < ALG_CHANNEL_NUM; ++ch) {
            if (pair_valid(ctx, ch, order)) {
                complex z = complex_mul_conj(iq[ch + order], iq[ch]);
                group_power[pair_group(ctx, ch, order)] += z.r * z.r + z.i * z.i;
            }
        }
    }
    for (uint8_t g = 0; g < group_count(ctx); ++g) {
        ctx->noise_mean += sqrtf(0.25f * (float)PI * group_power[g]);
        ctx->noise_var += (1.0f - 0.25f * (float)PI) * group_power[g];
    }
}

errcode_t subevent_motion_add_speed_path(SubeventMotionCtx* ctx, complex* iq)
{
    if ((ctx == NULL) || (iq == NULL) || (ctx->stage != STAGE_SPEED) || (ctx->path_count >= ctx->path_num)) {
        return ERRCODE_RANGING_ALG_INVALID_PARAM;
    }
    float first_speed = -(MOTION_MAX_SPEED + MOTION_SPEED_GUARD);

    add_noise_statistics(ctx, iq);

    // 速度谱：sum_group | sum_pair iq[ch + order] * conj(iq[ch]) * exp(j * rate * v) |
    // 组内相干累加，组间非相干累加；共轭乘积保留幅度，强信道权重更大
    for (uint8_t step = 0; step < SUBEVENT_MOTION_SPEED_POINT_NUM; ++step) {
        float current_speed = first_speed + (float)step * MOTION_SPEED_RESOLUTION;
        memset_s(ctx->group_acc, sizeof(ctx->group_acc), 0, sizeof(ctx->group_acc));
        for (uint8_t order = 1; order <= SUBEVENT_MOTION_PAIR_ORDER; ++order) {
            for (uint8_t ch = 0; ch < ALG_CHANNEL_NUM; ++ch) {
                if (!pair_valid(ctx, ch, order)) {
                    continue;
                }
                complex z = complex_mul_conj(iq[ch + order], iq[ch]);
                complex rotation = complex_expj(pair_speed_rate(ctx, ch, order) * current_speed);
                complex_acc(&ctx->group_acc[pair_group(ctx, ch, order)], complex_mul(z, rotation));
            }
        }
        float magnitude = 0.0f;
        for (uint8_t g = 0; g < group_count(ctx); ++g) {
            magnitude += complex_abs(ctx->group_acc[g]);
        }
        ctx->work.speed_spectrum[step] += magnitude;
    }
    ctx->path_count++;
    return ERRCODE_RANGING_ALG_SUCCESS;
}

errcode_t subevent_motion_solve_speed(SubeventMotionCtx* ctx)
{
    if ((ctx == NULL) || (ctx->stage != STAGE_SPEED) || (ctx->path_count != ctx->path_num)) {
        return ERRCODE_RANGING_ALG_INVALID_PARAM;
    }
    // 所有路径的 IQ 都为 0（例如远端 IQ 缺失），没有可用信号
    if (ctx->noise_var <= 0.0f) {
        return ERRCODE_RANGING_ALG_NOT_ENOUGH_IQ;
    }

    float* spectrum = ctx->work.speed_spectrum;
    uint8_t best = 0;
    for (uint8_t step = 1; step < SUBEVENT_MOTION_SPEED_POINT_NUM; ++step) {
        if (spectrum[step] > spectrum[best]) {
            best = step;
        }
    }

    // 抛物线插值细化峰值位置
    float speed_est = -(MOTION_MAX_SPEED + MOTION_SPEED_GUARD) + (float)best * MOTION_SPEED_RESOLUTION;
    bool at_edge = (best == 0) || (best == SUBEVENT_MOTION_SPEED_POINT_NUM - 1);
    if (!at_edge) {
        float left = spectrum[best - 1];
        float right = spectrum[best + 1];
        float den = left - 2.0f * spectrum[best] + right;
        if (den < 0.0f) {
            speed_est += 0.5f * (left - right) / den * MOTION_SPEED_RESOLUTION;
        }
    }

    // 分数低或峰值在搜索边界（速度可能超出范围）时，不补偿多普勒
    ctx->speed_est = speed_est;
    ctx->speed_score = (spectrum[best] - ctx->noise_mean) / sqrtf(ctx->noise_var);
    ctx->speed_valid = (ctx->speed_score >= MOTION_MIN_SPEED_SCORE) && !at_edge;
    ctx->speed = ctx->speed_valid ? speed_est : 0.0f;
    ctx->stage = STAGE_PHASE;
    return ERRCODE_RANGING_ALG_SUCCESS;
}

// ------------------------------------------------------------------------------------------------
// 第 2 遍：subevent 相位估计
// ------------------------------------------------------------------------------------------------
errcode_t subevent_motion_add_phase_path(SubeventMotionCtx* ctx, uint8_t path, complex* iq)
{
    if ((ctx == NULL) || (iq == NULL) || (ctx->stage != STAGE_PHASE) || (path >= ctx->path_num) ||
        ((ctx->path_added & (1U << path)) != 0)) {
        return ERRCODE_RANGING_ALG_INVALID_PARAM;
    }
    // 第一条路径加入时，速度谱已经用完，这块内存改存分组和
    if (ctx->path_added == 0) {
        memset_s(&ctx->work, sizeof(ctx->work), 0, sizeof(ctx->work));
    }

    // 多普勒补偿后的信道对共轭乘积，按分组累加（不修改 iq）
    complex* group_sum = ctx->work.group_sum[path];
    for (uint8_t order = 1; order <= SUBEVENT_MOTION_PAIR_ORDER; ++order) {
        for (uint8_t ch = 0; ch < ALG_CHANNEL_NUM; ++ch) {
            if (!pair_valid(ctx, ch, order)) {
                continue;
            }
            complex iq_lo = complex_mul(iq[ch], doppler_rotation(ctx, ch));
            complex iq_hi = complex_mul(iq[ch + order], doppler_rotation(ctx, ch + order));
            complex_acc(&group_sum[pair_group(ctx, ch, order)], complex_mul_conj(iq_hi, iq_lo));
        }
    }
    ctx->path_added |= (uint8_t)(1U << path);
    return ERRCODE_RANGING_ALG_SUCCESS;
}

// 分组和模型：group_sum[p][o][j][k] ~ |.| * exp(j * (c[p][o] + psi_k - psi_j))
//   c[p][o] : 路径 p、信道间隔 o 的截距（局部群时延）
//   psi_k : subevent k 的相位，psi_0 = 0
// 每一轮：
//   Q[j][k] = sum_{p,o} group_sum[p][o][j][k] * exp(-j * c[p][o])  ~ exp(j * (psi_k - psi_j))
//   H = Q + Q^H，对角线取各行 |H| 之和的最大值（K = 2 时不加这个移位，幂迭代会振荡）
//   x = H 的主特征向量（幂迭代）                                 ~ exp(-j * psi)
//   c[p][o] = arg(sum_{j,k} group_sum[p][o][j][k] * exp(-j * (psi_k - psi_j)))
static void estimate_subevent_phase(SubeventMotionCtx* ctx)
{
    uint8_t num = ctx->subevent_num;
    complex intercept_conj[SUBEVENT_MOTION_MAX_PATH_NUM][SUBEVENT_MOTION_PAIR_ORDER];
    complex h[SUBEVENT_MOTION_MAX_SUBEVENT_NUM][SUBEVENT_MOTION_MAX_SUBEVENT_NUM];
    complex x[SUBEVENT_MOTION_MAX_SUBEVENT_NUM];
    complex x_next[SUBEVENT_MOTION_MAX_SUBEVENT_NUM];

#define GROUP_SUM(p, o, j, k) (ctx->work.group_sum[(p)][((o) * num + (j)) * num + (k)])

    for (uint8_t k = 0; k < num; ++k) {
        ctx->subevent_phasor[k].r = 1.0f;
        ctx->subevent_phasor[k].i = 0.0f;
    }
    // 截距初值只用 subevent 内部的信道对，与 subevent 相位无关
    for (uint8_t p = 0; p < ctx->path_num; ++p) {
        for (uint8_t o = 0; o < SUBEVENT_MOTION_PAIR_ORDER; ++o) {
            complex diag = {0.0f, 0.0f};
            for (uint8_t k = 0; k < num; ++k) {
                complex_acc(&diag, GROUP_SUM(p, o, k, k));
            }
            intercept_conj[p][o] = complex_conj(complex_unit(diag));
        }
    }
    if (num == 1) {
        return;
    }

    for (uint8_t round = 0; round < PHASE_REFINE_NUM; ++round) {
        // Q（存放在 h 中）
        for (uint8_t j = 0; j < num; ++j) {
            for (uint8_t k = 0; k < num; ++k) {
                complex q = {0.0f, 0.0f};
                for (uint8_t p = 0; p < ctx->path_num; ++p) {
                    for (uint8_t o = 0; o < SUBEVENT_MOTION_PAIR_ORDER; ++o) {
                        complex_acc(&q, complex_mul(GROUP_SUM(p, o, j, k), intercept_conj[p][o]));
                    }
                }
                h[j][k] = q;
            }
        }
        // H = Q + Q^H，再对角线移位
        for (uint8_t j = 0; j < num; ++j) {
            h[j][j].r = 2.0f * h[j][j].r;
            h[j][j].i = 0.0f;
            for (uint8_t k = j + 1; k < num; ++k) {
                complex sum = { h[j][k].r + h[k][j].r, h[j][k].i - h[k][j].i };
                h[j][k] = sum;
                h[k][j] = complex_conj(sum);
            }
        }
        float diag_shift = 0.0f;
        for (uint8_t j = 0; j < num; ++j) {
            float row_sum = 0.0f;
            for (uint8_t k = 0; k < num; ++k) {
                row_sum += complex_abs(h[j][k]);
            }
            diag_shift = (row_sum > diag_shift) ? row_sum : diag_shift;
        }
        for (uint8_t j = 0; j < num; ++j) {
            h[j][j].r = diag_shift;
            h[j][j].i = 0.0f;
        }
        // 幂迭代，从当前估计 x = exp(-j * psi) 开始
        for (uint8_t k = 0; k < num; ++k) {
            x[k] = complex_conj(ctx->subevent_phasor[k]);
        }
        for (uint8_t iter = 0; iter < PHASE_POWER_ITER_NUM; ++iter) {
            for (uint8_t j = 0; j < num; ++j) {
                complex sum = {0.0f, 0.0f};
                for (uint8_t k = 0; k < num; ++k) {
                    complex_acc(&sum, complex_mul(h[j][k], x[k]));
                }
                x_next[j] = complex_unit(sum);
            }
            memcpy_s(x, sizeof(x), x_next, sizeof(x_next));
        }
        // exp(j * psi_k) = conj(x_k)，并使 psi_0 = 0
        for (uint8_t k = 0; k < num; ++k) {
            ctx->subevent_phasor[k] = complex_mul(complex_conj(x[k]), x[0]);
        }
        // 用所有分组更新截距
        for (uint8_t p = 0; p < ctx->path_num; ++p) {
            for (uint8_t o = 0; o < SUBEVENT_MOTION_PAIR_ORDER; ++o) {
                complex sum = {0.0f, 0.0f};
                for (uint8_t j = 0; j < num; ++j) {
                    for (uint8_t k = 0; k < num; ++k) {
                        // exp(-j * (psi_k - psi_j)) = exp(j * psi_j) * conj(exp(j * psi_k))
                        complex rotation = complex_mul_conj(ctx->subevent_phasor[j], ctx->subevent_phasor[k]);
                        complex_acc(&sum, complex_mul(GROUP_SUM(p, o, j, k), rotation));
                    }
                }
                intercept_conj[p][o] = complex_conj(complex_unit(sum));
            }
        }
    }
#undef GROUP_SUM
}

errcode_t subevent_motion_solve_phase(SubeventMotionCtx* ctx, SubeventMotionResult* res)
{
    if ((ctx == NULL) || (ctx->stage != STAGE_PHASE) || (ctx->path_added != (uint8_t)((1U << ctx->path_num) - 1))) {
        return ERRCODE_RANGING_ALG_INVALID_PARAM;
    }
    estimate_subevent_phase(ctx);
    ctx->stage = STAGE_DONE;

    if (res != NULL) {
        memset_s(res, sizeof(SubeventMotionResult), 0, sizeof(SubeventMotionResult));
        res->speed = ctx->speed;
        res->speed_est = ctx->speed_est;
        res->speed_score = ctx->speed_score;
        res->speed_valid = ctx->speed_valid;
        for (uint8_t k = 0; k < ctx->subevent_num; ++k) {
            res->subevent_phase[k] = atan2f(ctx->subevent_phasor[k].i, ctx->subevent_phasor[k].r);
        }
    }
    return ERRCODE_RANGING_ALG_SUCCESS;
}

// ------------------------------------------------------------------------------------------------
// 第 3 遍：IQ 补偿
// ------------------------------------------------------------------------------------------------
errcode_t subevent_motion_iq_compensation(SubeventMotionCtx* ctx, complex* iq)
{
    if ((ctx == NULL) || (iq == NULL) || (ctx->stage != STAGE_DONE)) {
        return ERRCODE_RANGING_ALG_INVALID_PARAM;
    }
    // iq[ch] *= exp(j * 4 * pi / c * f * v * tau) * exp(-j * psi[subevent])，未测量的信道保持原值
    for (uint8_t ch = 0; ch < ALG_CHANNEL_NUM; ++ch) {
        if (!channel_measured(ctx, ch)) {
            continue;
        }
        complex rotation = complex_mul_conj(doppler_rotation(ctx, ch),
                                            ctx->subevent_phasor[ctx->channel_subevent[ch]]);
        iq[ch] = complex_mul(iq[ch], rotation);
    }
    return ERRCODE_RANGING_ALG_SUCCESS;
}
