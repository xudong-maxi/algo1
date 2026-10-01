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

## C 模块（MCU 移植，协议无关）

- `c/jdps.h` / `c/jdps.c`：JDPS 的 C 实现（C99、float32、无动态内存，工作区 13.9 KB 由调用方提供）。
- 复数类型使用工程自带的 `complex`（成员 `r`、`i`），通过 `JDPS_COMPLEX_HEADER` 指定定义它的头文件；`c/test/complex_type.h` 仅供主机测试使用。
- 不绑定具体协议：BLE CS、星闪等都可以用。频点规划由 `cfg.chan0_freq_hz` / `cfg.chan_spacing_hz` 配置（默认值为 BLE CS 的 2402 MHz + idx × 1 MHz）；「segment」指使用同一套 AGC 设置的一段连续 step，对应 BLE CS 的 subevent。
- 入口：`jdps_process(&cfg, &meas, &work, &res)`，就地补偿 `meas.iq`，输出估计速度和各 segment 的相位。
- 回归测试（与 Python 参考实现逐点对比，仿真用例 + 实测 log）：`cd c && make test`
