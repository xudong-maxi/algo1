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
python -m cs_agc.log_replay data/*.txt --out log_plots          # 每个 procedure 一张图 + summary.csv
# 常用参数：--split 25,25,22  --se-gap-ms 40  --per-ant  --combine mul|conj  --n-mode0 1
```
- 解析 `ch_idx_list`（第 1 个为 mode0，跳过）与 4 路 `[iq_data]`（R=0 本地 / R=1 远端，按信道号索引），PBR 取 local×remote。
- log 打印截断造成的损坏条目会被丢弃（置 0，不参与配对和 IFFT），图标题里的 `missing` 为丢弃数。
- 每张图：每路天线一行，分别为补偿前相位、补偿后相位（去掉最强径斜率，按 subevent 着色）、IFFT 距离谱对比。

## C 模块（MCU 移植）

- `c/cs_jdps.h` / `c/cs_jdps.c`：JDPS 的 C 实现（C99、float32、无动态内存，工作区 13.9 KB 由调用方提供）。
- 入口：`cs_jdps_process(&cfg, &meas, &work, &res)`，就地补偿 `meas.iq`，输出估计速度和各 subevent 相位。
- 回归测试（与 Python 参考实现逐点对比，仿真用例 + 实测 log）：`cd c && make test`
