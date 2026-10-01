"""Export test vectors (input + Python JDPS reference output) for the C module.

Usage: python -m cs_agc.export_vectors c/test/vectors.txt
Cases: random simulated procedures (K = 1..5, both AGC modes) + the sample log.
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


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "c/test/vectors.txt"
    rng = np.random.default_rng(123)
    with open(out, "w") as fp:
        i = 0
        for K in [1, 2, 3, 3, 3, 4, 5]:
            for per_ant in [False, True]:
                cfg = SimCfg(n_se=K, snr_db=float(rng.uniform(10, 30)), channel=SCENARIOS["MP_mod"],
                             agc_shared_across_ant=not per_ant)
                write_case(fp, f"sim{i}_K{K}_perant{int(per_ant)}", generate(rng, cfg),
                           JdpsCfg(per_ant_theta=per_ant))
                i += 1
        # robustness cases
        cfg = SimCfg(n_se=3)
        ms = generate(rng, cfg)
        ms.y = rng.standard_normal(ms.y.shape) + 1j * rng.standard_normal(ms.y.shape)
        write_case(fp, "noise_only_K3", ms, JdpsCfg())                 # expect v_valid = 0
        write_case(fp, "beyond_range_v15", generate(rng, cfg, v=15.0), JdpsCfg())  # expect edge
        write_case(fp, "low_snr_m3dB", generate(rng, SimCfg(n_se=3, snr_db=-3.0)), JdpsCfg())
        tc = dict(t_meas=565e-6, t_gap=150e-6, t_step=715e-6, se_gap=40e-3)
        proc = parse_log("data/sample_log.txt")[0]
        ms = build_measurement(proc, [25, 25, 22], 1, tc)
        write_case(fp, "log_proc7121", ms, JdpsCfg())
        write_case(fp, "log_proc7121_perant", ms, JdpsCfg(per_ant_theta=True))
        proc["iq"] = {k: v for k, v in proc["iq"].items() if k[1] == 0}   # remote IQ missing
        write_case(fp, "log_no_remote_iq", build_measurement(proc, [25, 25, 22], 1, tc), JdpsCfg())
    print("wrote", out)


if __name__ == "__main__":
    main()
