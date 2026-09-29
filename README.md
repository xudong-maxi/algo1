# algo1

BLE Channel Sounding 多 subevent AGC 相位 + 多普勒联合补偿（JDPS）。

- 设计文档：[docs/multi_subevent_agc_design.md](docs/multi_subevent_agc_design.md)
- 仿真代码：`cs_agc/`（`sim_model.py` 测量模型，`algorithms.py` 算法，`ranging.py` IFFT 测距，`run_sim.py` Monte Carlo）
- 仿真结果：`results/`

```bash
pip install -r requirements.txt
python -m cs_agc.run_sim --trials 400 --out results   # 4 核约 3 分钟
```
