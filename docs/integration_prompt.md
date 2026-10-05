# 集成提示词：把 subevent_motion_alg 接入本地工程

> 用法：把下面「提示词正文」整段发给负责集成的 agent，并让它能访问本地工程和本仓库（`xudong-maxi/algo1`，分支 `claude/bluetooth-handling-d3l8j4`）。

---

## 提示词正文

你要把一个已经验证过的算法模块 **subevent_motion_alg**（多 subevent 的运动补偿 + AGC 相位对齐，算法名 JDPS）集成到我的嵌入式工程里，替换或并行于现有的 `motion_correct_alg`（单 subevent 运动补偿）。这个模块已经按工程的 `motion_correct_alg` 风格编写（`channel_select_t` 输入、按信道号排列的 `complex iq[ALG_CHANNEL_NUM]`、`errcode_t`、`memset_s`、双语 Doxygen 注释），你的工作是：一，接入工程的数据来源、调用流程和构建；二，把剩余与工程规范不一致的地方改过来。**算法的数学过程和数值行为不允许改变**，改完必须通过回归测试证明这一点。

### 1. 参考材料（来自仓库 xudong-maxi/algo1，分支 claude/bluetooth-handling-d3l8j4）
- `c/subevent_motion_alg.h`、`c/subevent_motion_alg.c`：要接入的 C 实现（C99、float32、无 malloc、流式接口）。
- `c/test/subevent_motion_test.c`、`c/test/port/common_util.h`、`c/test/port/securec.h`、`c/Makefile`：主机回归测试。`port/` 下是工程头文件的**测试桩**（只定义了模块用到的 `complex`、`channel_select_t`、`errcode_t`、`ERRCODE_RANGING_ALG_*`、`PI`、`LIGHT_SPEED`、`TWO`、`ALG_CHANNEL_NUM`、`memset_s`、`memcpy_s`）。`make test` 会调用 Python 生成测试向量 `c/test/vectors.txt`，按 `SUBEVENT_MOTION_PAIR_ORDER` = 1 和 2 各编译一次，与 Python 参考实现逐点对比（每种 15 个用例）。
- `cs_agc/export_vectors.py`、`cs_agc/algorithms.py`：Python 参考实现和测试向量生成脚本（需要 numpy）。
- `cs_agc/log_replay.py`：实测 log 的解析方式，用来对照 IQ 和时序的组织方式。
- `docs/multi_subevent_agc_design.md`：算法原理（§2 信号模型，§4 各 Stage，§7 C 实现、内存与复杂度，§9 已确认的硬件条件）。
- 工程中现有的 `motion_correct_alg.h/.c`：风格和数据约定的参照。

开始前先把这些文件读完。

### 2. 已确认的系统条件（按这些条件适配，不要再去猜）
- 协议：BLE Channel Sounding，Mode‑2 PBR，随机跳频，72 个信道（ch 2–22、26–76，f = 2402 + ch MHz；模块里由 `CHANNEL0_FREQ`、`CHANNEL_SPACING` 定义）。
- `channel_select_t`：每个 subevent 的第一步是 mode0，`ch_num_per_subevent` 包含它（例如 26, 26, 23）；`time_per_channel` 按时间顺序给出每一步的时长（us），包含每个 subevent 的 mode0（共 75 项，典型值 mode0 483、mode‑2 715/823）；`ch_hop_orders` 与 `time_per_channel` 逐步对应，每个 subevent 开头为 mode0（`ch_num` = 75；模块也兼容 log 打印格式：只有第 0 项为 mode0，`ch_num` = 73）；`t_mes` 为上一个 subevent 结束到下一个 subevent 开始的间隔（us，典型 40 ms 左右）。其他 `ch_num` 会返回参数错误。注意：原 `motion_correct_alg` 在这个约定下时间戳会错位（见设计文档 §10.5），不要照搬它的时间计算。
- IQ：每条天线路径一组，**本地 IQ × 远端 IQ（直接相乘，不取共轭）**，按信道号排列，长度 `ALG_CHANNEL_NUM`；未测量的信道由跳频表判定，值不参与计算。
- 天线：4 条路径 = 发起端 2 根 × 反射端 2 根（两根天线互相垂直）；path0=(0,0)、path1=(0,1)、path2=(1,0)、path3=(1,1)。
- AGC 按设备设定，每个 subevent 只设一次 → 所有路径共用一组 subevent 相位（模块只支持这种模式）；速度由 4 条路径统一估计。
- 平台：MCU 128 MHz，带单精度 FPU。AGC + 运动补偿处理（不含测距）4 条路径合计必须 < 60 ms；估算默认 `PAIR_ORDER=2` 约 36.4 ms、`PAIR_ORDER=1` 约 18.6 ms（K=3）。
- 内存：**任何时刻只能有一路 IQ 在内存里**，但同一路 IQ 可以多次读取。必须保留流式调用：每路 `add_speed_path` → `solve_speed` → 每路 `add_phase_path` → `solve_phase` → 每路 `iq_compensation`。`SubeventMotionCtx` 1,516 B（默认），栈 ≤ 448 B。
- 复数类型：工程里已有 `complex`，成员**依次**为 `float r`、`float i`；模块用 `{ r, i }` 初始化，依赖这个顺序。

### 3. 必须保持不变的东西（算法契约）
以下任何一项改动都会改变结果，**禁止修改**。如果你认为某项必须改，先停下来问我。
1. 配对与分组：信道 (ch, ch+o)（o = 1..`SUBEVENT_MOTION_PAIR_ORDER`）都被测量才成对；分组号 ((o−1)·K + ch 所在 subevent)·K + (ch+o 所在 subevent)；配对乘积 `iq[ch+o] * conj(iq[ch])`（保留幅度加权，不能只取相位）。
2. 相位率 `pair_speed_rate`：频率减去参考频率（平均信道号），时间减去 t_ref（平均测量时刻）和各 subevent 的平均时刻。这是为了 float32 精度，不能化简回绝对频率乘绝对时间。
3. 速度搜索：范围 ±(10 + 1) m/s、步长 0.25 m/s、89 个点；指标是各组和的模对组、对路径求和；峰值用抛物线插值。
4. 可信度判断：score = (峰值 − μ) / σ，μ = Σ sqrt(π/4·S_g)，σ² = Σ (1 − π/4)·S_g。score < 5 或峰值落在搜索最外侧时速度改为 0；所有 IQ 为 0 时 `subevent_motion_solve_speed` 返回 `ERRCODE_RANGING_ALG_NOT_ENOUGH_IQ`。
5. 多普勒补偿：`exp(+j·4π/c · f · v · τ)`，使用每个信道自己的频率 f（同时补偿距离迁移）；相位先用 `remainderf` 折回 [−π, π] 再算 sin/cos。
6. 相位估计：截距初值只用 subevent 内部的信道对；`H = Q + Q^H`，对角线先取 2·Re(Q_kk)，再整体设为各行 |H| 之和的最大值（**这个对角线移位不能删，否则 K=2 时幂迭代会振荡**）；幂迭代每步逐元素归一化；结果使 ψ_0 = 0；共 3 轮，每轮 8 次幂迭代，每轮更新截距。

### 4. 可以按工程规范调整的东西
- 文件位置、`#include` 列表与顺序、编译宏的默认值放在哪里（模块头文件或工程配置头文件）。
- 工程实际的错误码值、日志宏、断言宏；如果工程有更合适的错误码（例如区分「无有效信号」），可以替换并告诉我。
- `SUBEVENT_MOTION_MAX_PATH_NUM`、`SUBEVENT_MOTION_MAX_SUBEVENT_NUM`（subevent 最多 3 个时设为 3，ctx 1,052 B）；`SUBEVENT_MOTION_PAIR_ORDER` 默认 2，改为 1 前先问我。
- 数学函数如果工程有 CMSIS-DSP 或自研快速实现，可以替换，但替换后必须仍通过回归测试。
- 调用侧胶水代码：从工程的 CS 结果里取出 `channel_select_t` 和每条路径的 IQ（本地 × 远端），按流式顺序调用；补偿后的 IQ 交给工程现有的 IFFT 测距流程；速度、分数、subevent 相位按工程习惯输出到日志或结果结构。
- 只做与集成有关的改动，不要顺手重构工程里其他代码。

### 5. 工作步骤
1. **调研工程，先不要写代码**：找到并列出——`channel_select_t` 的填充位置、每条路径本地/远端 IQ 的存储和读取方式、现有 `motion_effect_analyze` / `motion_iq_compensation` 的调用点、IFFT 测距入口、`common_util.h` 中相关定义的实际值（`ALG_CHANNEL_NUM`、错误码、`PI`、`LIGHT_SPEED`、`TWO`）、构建系统。
2. **给我一份接入方案后等我确认**：模块放在哪、调用侧数据流（每一遍从哪里读 IQ、结果写回哪里）、与现有 `motion_correct_alg` 的关系（替换还是开关切换）、内存预算、需要修改的命名或错误码。**确认后再开始改代码。**
3. 把模块放进工程，按第 4 节调整；算法部分逐行对应原实现。
4. 回归测试：主机上用 `c/test/port/` 的测试桩编译，或改为直接包含工程的真实头文件；`test/vectors.txt` 用 `python -m cs_agc.export_vectors` 从本仓库生成，不要手改。
5. 写调用侧胶水代码，接入测距流程。
6. 在目标板或模拟器上编译；如果可以，用 DWT 周期计数器测量 4 条路径处理的总耗时。

### 6. 验收标准
- 主机回归测试全部通过（`PAIR_ORDER` = 1 和 2 各 15 个用例），容差不放宽：速度差 < 0.01 m/s，subevent 相位差 < 0.01 rad，补偿后 IQ 的相对误差 < 1e-3，返回状态和可信标志完全一致，score 相对误差 < 1e-3。
- 目标平台编译无新增警告；不使用 malloc；IQ 内存峰值为一路。
- 4 条路径处理总耗时 < 60 ms（给出实测值或基于实测的估算）。
- 如果有实测 log，用本仓库的 `python -m cs_agc.log_replay` 和工程实现分别处理同一份数据，速度、subevent 相位和测距结果应一致（log_replay 默认使用间隔 1 和 2 的配对，与默认的 `PAIR_ORDER=2` 一致）。

### 7. 最后交付给我
1. 改动文件清单和简要说明；
2. 回归测试输出（完整表格）；
3. 目标平台的内存占用（ctx、栈、IQ 缓冲）和耗时；
4. 遗留问题和需要我决定的事项。

遇到以下情况先停下来问我，不要自己假设：工程里的数据和第 2 节的条件对不上（例如 IQ 是共轭相乘、时间单位不同、mode0 位置不同）；为满足工程约束必须改动第 3 节的任何一项；回归测试无法在不放宽容差的情况下通过。
