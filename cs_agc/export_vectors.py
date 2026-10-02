"""Export test vectors (input + Python JDPS reference output) for c/subevent_motion_alg.c.

Usage: python -m cs_agc.export_vectors c/test/vectors.txt
Cases: random simulated procedures (K = 1..4), robustness cases and the sample log,
for pair order 1 and 2.
"""
import sys

import numpy as np

from .algorithms import JdpsCfg, jdps
from .log_replay import build_measurement, parse_log
from .sim_model import SCENARIOS, SimCfg, generate


def write_case(fp, name, ms, cfg: JdpsCfg):
    y_c, info = jdps(ms, cfg)
    A, N = ms.y.shape
    K = int(ms.se.max()) + 1
    psi = np.broadcast_to(np.atleast_2d(info["psi"]), (A, K))
    fp.write(f"case {name}\n")
    fp.write(f"cfg {cfg.v_max} {cfg.v_guard} {cfg.v_step} {cfg.min_v_score} {len(cfg.orders)} "
             f"{cfg.n_power_iter} {cfg.n_refine_iter} {int(cfg.per_ant_theta)}\n")
    fp.write(f"dims {A} {N} {K}\n")
    fp.write("chan " + " ".join(str(int(c)) for c in ms.ch) + "\n")
    fp.write("se " + " ".join(str(int(s)) for s in ms.se) + "\n")
    fp.write("time " + " ".join(f"{t:.9e}" for t in ms.t) + "\n")
    fp.write("iq " + " ".join(f"{z.real:.9e} {z.imag:.9e}" for z in ms.y.ravel()) + "\n")
    status = 0 if info["status"] == "ok" else 3          # JDPS_OK / JDPS_ERR_NO_SIGNAL
    v_est = info["v_hat"] if np.isfinite(info["v_hat"]) else 0.0
    fp.write(f"exp_status {status} {int(info['v_valid'])} {info['v_score']:.9e} {v_est:.9e}\n")
    fp.write(f"exp_v {info['v_used']:.9e}\n")
    fp.write("exp_psi " + " ".join(f"{p:.9e}" for p in psi.ravel()) + "\n")
    fp.write("exp_iq " + " ".join(f"{z.real:.9e} {z.imag:.9e}" for z in y_c.ravel()) + "\n")


# C 模块 (c/subevent_motion_alg.c) 的配置：所有路径共用 subevent 相位；
# 配对阶数对应编译宏 SUBEVENT_MOTION_PAIR_ORDER（1 或 2），两种都导出，测试程序只跑匹配的用例
C_CFGS = {1: JdpsCfg(orders=(1,), per_ant_theta=False), 2: JdpsCfg(orders=(1, 2), per_ant_theta=False)}


def write_all_cases(fp, order, cfg):
    rng = np.random.default_rng(123)
    for i, K in enumerate([1, 2, 2, 3, 3, 3, 3, 4, 4]):
        sim = SimCfg(n_se=K, snr_db=float(rng.uniform(10, 30)), channel=SCENARIOS["MP_mod"])
        write_case(fp, f"o{order}_sim{i}_K{K}", generate(rng, sim), cfg)
    # robustness cases
    sim = SimCfg(n_se=3)
    ms = generate(rng, sim)
    ms.y = rng.standard_normal(ms.y.shape) + 1j * rng.standard_normal(ms.y.shape)
    write_case(fp, f"o{order}_noise_only_K3", ms, cfg)                 # expect speed not valid
    write_case(fp, f"o{order}_beyond_range_v15", generate(rng, sim, v=15.0), cfg)
    write_case(fp, f"o{order}_low_snr_m3dB", generate(rng, SimCfg(n_se=3, snr_db=-3.0)), cfg)
    write_case(fp, f"o{order}_severe_mp_snr5",
               generate(rng, SimCfg(n_se=3, snr_db=5.0, channel=SCENARIOS["MP_severe"])), cfg)
    tc = dict(t_meas=565e-6, t_gap=150e-6, t_step=715e-6, se_gap=40e-3, mode0=483e-6)
    proc = parse_log("data/sample_log.txt")[0]
    write_case(fp, f"o{order}_log_proc7121", build_measurement(proc, [25, 25, 22], 1, tc), cfg)
    proc["iq"] = {k: v for k, v in proc["iq"].items() if k[1] == 0}   # remote IQ missing
    write_case(fp, f"o{order}_log_no_remote_iq", build_measurement(proc, [25, 25, 22], 1, tc), cfg)


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "c/test/vectors.txt"
    with open(out, "w") as fp:
        for order, cfg in C_CFGS.items():
            write_all_cases(fp, order, cfg)
    print("wrote", out)


if __name__ == "__main__":
    main()
