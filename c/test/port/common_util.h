/*
 * 主机测试用的工程接口桩：只包含 subevent_motion_alg 用到的工程定义。
 * 集成到工程时不需要这个文件，使用工程自己的 common_util.h。
 */
#ifndef COMMON_UTIL_H
#define COMMON_UTIL_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#define PI              3.14159265358979323846
#define LIGHT_SPEED     299792458.0
#define TWO             2
#define ALG_CHANNEL_NUM 80

typedef uint32_t errcode_t;
#define ERRCODE_RANGING_ALG_SUCCESS         0
#define ERRCODE_RANGING_ALG_INVALID_PARAM   1
#define ERRCODE_RANGING_ALG_NOT_ENOUGH_IQ   2
#define ERRCODE_RANGING_ALG_MALLOC_FAIL     3

typedef struct {
    float r;
    float i;
} complex;

typedef struct {
    uint8_t ch_num;
    uint8_t subevent_num;
    uint16_t t_mes;
    uint8_t* ch_hop_orders;
    uint8_t* antenna_path_permutation_index;
    uint8_t* ch_num_per_subevent;
    uint16_t* time_per_channel;
} channel_select_t;

#endif  /* COMMON_UTIL_H */
