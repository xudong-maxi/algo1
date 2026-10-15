# 固件输入 dump：定位「固件结果与 log_replay 不一致」

思路：把固件**实际传给 `subevent_motion_alg` 的输入**原样打印出来，然后用同一份输入分别跑：
- `python -m cs_agc.replay_dump fw.log`：Python 参考实现；
- `c/replay_dump fw.log`：C 模块。可以把工程里移植后的 `subevent_motion_alg.c` 换进来编译。

判断方法：
| 现象 | 结论 |
|---|---|
| 固件结果 = C 回放 = Python 回放 | 模块和移植都没问题；与 log_replay 的差异来自输入（IQ 合成方式、`channel_select_cfg`、定点缩放） |
| C 回放 = Python 回放 ≠ 固件结果 | 固件里的模块与仓库版本不一致（移植改动、编译宏、`PI` / `LIGHT_SPEED` 等定义） |
| C 回放 ≠ Python 回放 | 把移植后的 `.c` 换进 `c/` 编译时出现差异，说明移植改变了算法 |

## 打印代码（放在 `cs_jdps_preprocess` 中）

```c
#if defined(ALG_TEST) || defined(CI_MODE)
static void jdps_dump_cfg(const channel_select_t* cfg)
{
    uint16_t step_num = 0;
    printf("[jdps_dump] begin\n");
    printf("[jdps_dump] cfg %u %u %u\n", (unsigned)cfg->ch_num, (unsigned)cfg->subevent_num, (unsigned)cfg->t_mes);
    printf("[jdps_dump] hop");
    for (uint16_t i = 0; i < cfg->ch_num; i++) {
        printf(" %u", (unsigned)cfg->ch_hop_orders[i]);
    }
    printf("\n[jdps_dump] per_subevent");
    for (uint8_t s = 0; s < cfg->subevent_num; s++) {
        printf(" %u", (unsigned)cfg->ch_num_per_subevent[s]);
        step_num += cfg->ch_num_per_subevent[s];
    }
    printf("\n[jdps_dump] time");
    for (uint16_t i = 0; i < step_num; i++) {
        printf(" %u", (unsigned)cfg->time_per_channel[i]);
    }
    printf("\n");
}

static void jdps_dump_iq(uint8_t path, const complex* iq)
{
    printf("[jdps_dump] path %u", (unsigned)path);
    for (uint16_t c = 0; c < ALG_CHANNEL_NUM; c++) {
        printf(" %.6e %.6e", (double)iq[c].r, (double)iq[c].i);
    }
    printf("\n");
}
#endif
```

调用位置：
1. `subevent_motion_init` 之前：`jdps_dump_cfg(channel_select_cfg);`
2. pass 1 循环里，`calc_iq` 之后、`subevent_motion_add_speed_path` 之前：`jdps_dump_iq(s, iq_buf);`
3. 现有的 `[jdps] v=...` 打印之后：`printf("[jdps_dump] end\n");`

每次测量约 4 × 80 × 2 个浮点数，log 会比较长，调试时只抓几次测量即可。
