# algo1

BLE Channel Sounding 多 subevent AGC 相位 + 多普勒联合补偿（JDPS）。

- 设计文档：[docs/multi_subevent_agc_design.md](docs/multi_subevent_agc_design.md)
- 仿真代码：`cs_agc/`（`sim_model.py` 测量模型，`algorithms.py` 算法，`ranging.py` IFFT 测距，`run_sim.py` Monte Carlo）
- 仿真结果：`results/`

```bash
pip install -r requirements.txt
python -m cs_agc.run_sim --trials 400 --out results   # 4 核约 3 分钟
```

## 实测 log 回放

```bash
python -m cs_agc.log_replay data/*.txt --out log_plots          # 每个 procedure 一张图 + CSV
# 常用参数：--split 25,25,22  --se-gap-ms 40  --per-ant  --combine mul|conj  --n-mode0 1
```
- 解析 `ch_idx_list`（第 1 个为 mode0，跳过）与 4 路 `[iq_data]`（R=0 本地 / R=1 远端，按信道号索引），PBR 取 local×remote。
- log 打印截断造成的损坏条目会被丢弃（置 0，不参与配对和 IFFT），图标题里的 `missing` 为丢弃数。
- 每张图：每路天线一行，分别为补偿前相位、补偿后相位（去掉最强径斜率，按 subevent 着色）、IFFT 距离谱对比；右上角为速度谱，标题给出 v 是否可信（score、是否落在边界）。
- 输出文件：
  - `summary.csv`：每个 procedure 一行。`d_before`（优化前，不做任何补偿）、`d_after`（JDPS 补偿后）；`v`（JDPS 实际用于补偿的速度）、`v_est / v_valid / v_score`（速度原始估计、是否可信、可信度分数）；`psi`（各 subevent 被去掉的相位）。
  - `per_antenna.csv`：每个 procedure、每路天线一行，保存优化前后的测距值。
  - `trend.png`：多于 1 个 procedure 时生成，画出距离和速度随 procedure 的变化。
- 远端 IQ 按 log 中 `local iq ts / remote iq ts` 的对应关系配对；某路天线缺本地或远端 IQ 时打印 warn，全部缺失（全 0）时跳过该次测量。

## C 模块（MCU，按工程 motion_correct_alg 风格）

- `c/subevent_motion_alg.h` / `c/subevent_motion_alg.c`：多 subevent 运动补偿 + AGC 相位对齐。输入为工程的 `channel_select_t` 和按信道号排列的 IQ（`complex iq[ALG_CHANNEL_NUM]`），返回 `errcode_t`；所有天线路径共用 subevent 相位（AGC 按设备设定）。
- 内存：`SubeventMotionCtx` 1,124 B（默认：最多 4 条路径、4 个 subevent、只用相邻信道）+ 栈 ≤ 368 B；任何时刻只需一路 IQ 在内存中。
- 编译宏：`SUBEVENT_MOTION_MAX_PATH_NUM`、`SUBEVENT_MOTION_MAX_SUBEVENT_NUM`、`SUBEVENT_MOTION_PAIR_ORDER`（1 或 2，2 更稳，内存约 +0.6 KB、耗时约 ×2）。
- 调用顺序（每路 IQ 按顺序提供 3 遍）：
  ```c
  subevent_motion_init(&ctx, channel_select_cfg, path_num);
  for (p = 0; p < path_num; p++) subevent_motion_add_speed_path(&ctx, iq_p);      /* 第 1 遍 */
  if (subevent_motion_solve_speed(&ctx) != ERRCODE_RANGING_ALG_SUCCESS) { /* 无有效信号，跳过 */ }
  for (p = 0; p < path_num; p++) subevent_motion_add_phase_path(&ctx, p, iq_p);   /* 第 2 遍 */
  subevent_motion_solve_phase(&ctx, &res);
  for (p = 0; p < path_num; p++) subevent_motion_iq_compensation(&ctx, iq_p);     /* 第 3 遍，之后做 IFFT */
  ```
- 回归测试（与 Python 参考实现逐点对比，PAIR_ORDER = 1 / 2 两种编译各跑一遍）：`cd c && make test`。`c/test/port/` 下是工程头文件（`common_util.h`、`securec.h`）的主机测试桩，集成时使用工程自己的头文件。
