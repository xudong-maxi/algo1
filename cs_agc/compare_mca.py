"""Compare the project's existing motion_correct_alg with JDPS (C module configurations).

Usage: python -m cs_agc.compare_mca [--trials 400] [--out results/compare_mca]
Methods (all ranged with the same IFFT back-end):
    genie      : true speed + true AGC phases (upper bound)
    mca        : existing motion_correct_alg, per antenna path, no AGC handling
    jdps_o1    : JDPS, adjacent channel pairs only   (SUBEVENT_MOTION_PAIR_ORDER = 1)
    jdps_o2    : JDPS, pair spacing 1 and 2          (SUBEVENT_MOTION_PAIR_ORDER = 2, default)
"""
import argparse
import dataclasses
import json
import os
from multiprocessing import Pool

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .algorithms import JdpsCfg, genie, jdps, motion_correct_alg
from .log_replay import build_measurement, parse_log
from .ranging import delay_profile, estimate_distance, to_grid
from .sim_model import SCENARIOS, SimCfg, generate

METHODS = ["genie", "mca", "jdps_o1", "jdps_o2"]
LABELS = {"genie": "Genie (true v, true AGC)", "mca": "Existing motion_correct_alg",
          "jdps_o1": "JDPS, PAIR_ORDER=1", "jdps_o2": "JDPS, PAIR_ORDER=2 (default)"}
COLORS = {"genie": "#222222", "mca": "#d62728", "jdps_o1": "#ff7f0e", "jdps_o2": "#1f77b4"}
CFG_O1, CFG_O2 = JdpsCfg(orders=(1,)), JdpsCfg(orders=(1, 2))


def run_methods(ms):
    out = {"genie": (genie(ms), ms.v)}
    y, sp = motion_correct_alg(ms)
    out["mca"] = (y, sp)                                     # per-path speeds
    for k, cfg in (("jdps_o1", CFG_O1), ("jdps_o2", CFG_O2)):
        y, info = jdps(ms, cfg)
        out[k] = (y, info["v_used"])
    return out


def phase_residual_rms(yc, yg):
    r = yc * np.conj(yg)
    r = r / (np.abs(r) + 1e-30)
    r = r * np.exp(-1j * np.angle(r.sum()))
    return float(np.sqrt(np.mean(np.angle(r) ** 2)))


def one_trial(args):
    seed, cfg = args
    ms = generate(np.random.default_rng(seed), cfg)
    res = run_methods(ms)
    yg = res["genie"][0]
    d_g = estimate_distance(yg, ms.ch)
    pk_g = delay_profile(to_grid(yg, ms.ch))[1].max()
    rows = {}
    for k, (yc, v) in res.items():
        d = estimate_distance(yc, ms.ch)
        v_err = np.abs(np.atleast_1d(v) - ms.v)
        rows[k] = dict(err=abs(d - ms.d0), dev=abs(d - d_g), ph=phase_residual_rms(yc, yg),
                       loss=-10 * np.log10(delay_profile(to_grid(yc, ms.ch))[1].max() / pk_g),
                       verr=float(np.nanmax(v_err)) if np.isfinite(v_err).any() else np.nan)
    return rows


def monte_carlo(cfg, n, pool, seed0):
    rows = pool.map(one_trial, [(seed0 + i, cfg) for i in range(n)], chunksize=8)
    return {k: {m: np.array([r[k][m] for r in rows]) for m in rows[0][k]} for k in METHODS}


def stats(r):
    p = lambda x, q: float(np.nanpercentile(x, q))
    return dict(p50=p(r["err"], 50), p90=p(r["err"], 90), out1m=float(np.mean(r["err"] > 1) * 100),
                dev90=p(r["dev"], 90), dev_out=float(np.mean(r["dev"] > 0.5) * 100), ph90=p(r["ph"], 90),
                loss90=p(r["loss"], 90), v90=p(r["verr"], 90))


def sweep_plot(xs, tabs, xlabel, title, path):
    keys = [("p90", "P90 |range error| [m]"), ("dev90", "P90 |d - d_genie| [m]"),
            ("ph90", "P90 residual phase RMS [rad]"), ("v90", "P90 speed error [m/s] (worst path)")]
    fig, ax = plt.subplots(1, 4, figsize=(18, 4.2))
    for (key, lab), axx in zip(keys, ax):
        for m in METHODS:
            if m == "genie" and key != "p90":
                continue
            axx.plot(xs, [t[m][key] for t in tabs], "o-", color=COLORS[m], label=LABELS[m],
                     lw=2.2 if m == "jdps_o2" else 1.3)
        axx.set_xlabel(xlabel), axx.set_title(lab, fontsize=10), axx.grid(alpha=0.3)
        if key in ("dev90", "ph90", "v90"):
            axx.set_yscale("log")
    ax[0].legend(fontsize=7)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def cdf_plot(r, title, path):
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    for m in METHODS:
        e = np.sort(r[m]["err"])
        ax[0].plot(e, np.arange(1, len(e) + 1) / len(e), color=COLORS[m], label=LABELS[m],
                   lw=2.2 if m == "jdps_o2" else 1.3, ls="--" if m == "genie" else "-")
        if m != "genie":
            e = np.sort(r[m]["dev"])
            ax[1].plot(e, np.arange(1, len(e) + 1) / len(e), color=COLORS[m], label=LABELS[m],
                       lw=2.2 if m == "jdps_o2" else 1.3)
    ax[0].set_xscale("log"), ax[0].set_xlim(0.01, 30), ax[0].set_xlabel("|range error| [m]")
    ax[0].set_ylabel("CDF"), ax[0].set_title(title), ax[0].legend(fontsize=8), ax[0].grid(alpha=0.3)
    ax[1].set_xscale("log"), ax[1].set_xlim(1e-3, 30), ax[1].set_xlabel("|d - d_genie| [m]")
    ax[1].set_title("Deviation from genie (compensation-induced error)"), ax[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def real_log(out):
    tc = dict(t_meas=565e-6, t_gap=150e-6, t_step=715e-6, se_gap=40e-3)
    ms = build_measurement(parse_log("data/sample_log.txt")[0], [25, 25, 22], 1, tc)
    y_mca, sp = motion_correct_alg(ms)
    rows = {"raw": dict(d=estimate_distance(ms.y, ms.ch), v="-"),
            "mca": dict(d=estimate_distance(y_mca, ms.ch), v=" / ".join(f"{v:+.2f}" for v in sp))}
    fig, ax = plt.subplots(figsize=(8, 4))
    for k, y, c in [("raw", ms.y, "#999999"), ("mca", y_mca, COLORS["mca"])]:
        dd, pp = delay_profile(to_grid(y, ms.ch))
        ax.plot(dd, 10 * np.log10(pp / pp.max()), color=c, label=k)
    for k, cfg in (("jdps_o1", CFG_O1), ("jdps_o2", CFG_O2)):
        y, info = jdps(ms, cfg)
        rows[k] = dict(d=estimate_distance(y, ms.ch), v=f"{info['v_used']:+.2f} (score {info['v_score']:.1f})")
        dd, pp = delay_profile(to_grid(y, ms.ch))
        ax.plot(dd, 10 * np.log10(pp / pp.max()), color=COLORS[k], label=LABELS[k])
    ax.set_xlim(0, 30), ax.set_ylim(-30, 1), ax.grid(alpha=0.3), ax.legend(fontsize=8)
    ax.set_xlabel("distance [m] (uncalibrated)"), ax.set_title("Sample log proc 7121: IFFT profile")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "real_log_profile.png"), dpi=120)
    plt.close(fig)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=400)
    ap.add_argument("--out", default="results/compare_mca")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    N, S = a.trials, {}
    base = SimCfg(n_se=3)
    with Pool(os.cpu_count()) as pool:
        r = monte_carlo(base, N, pool, 100)
        S["baseline"] = {m: stats(r[m]) for m in METHODS}
        cdf_plot(r, "K=3, MP_mod, SNR 20 dB, v~U[-10,10]", os.path.join(a.out, "cdf_K3.png"))
        sweeps = {
            "K": ([1, 2, 3, 4], lambda x: dataclasses.replace(base, n_se=x), "number of subevents K",
                  "MP_mod, SNR 20 dB, v~U[-10,10]"),
            "speed": ([0, 2, 5, 8, 10], lambda x: dataclasses.replace(base, v_range=(x, x)), "|v| [m/s]",
                      "K=3, MP_mod, SNR 20 dB"),
            "snr": ([0, 5, 10, 20, 30], lambda x: dataclasses.replace(base, snr_db=x), "SNR per tone [dB]",
                    "K=3, MP_mod, v~U[-10,10]"),
        }
        for name, (xs, mk, xl, title) in sweeps.items():
            tabs = [{m: stats(v) for m, v in monte_carlo(mk(x), N, pool, 1000 + 97 * i).items()}
                    for i, x in enumerate(xs)]
            S[name] = dict(x=xs, table=tabs)
            sweep_plot(xs, tabs, xl, title, os.path.join(a.out, f"sweep_{name}.png"))
        S["scenarios"] = {}
        for i, (name, cfg) in enumerate({
                "LOS": dataclasses.replace(base, channel=SCENARIOS["LOS"]),
                "MP_severe": dataclasses.replace(base, channel=SCENARIOS["MP_severe"]),
                "K=1 (single subevent)": dataclasses.replace(base, n_se=1)}.items()):
            r = monte_carlo(cfg, N, pool, 5000 + 97 * i)
            S["scenarios"][name] = {m: stats(r[m]) for m in METHODS}
    S["real_log"] = real_log(a.out)
    with open(os.path.join(a.out, "summary.json"), "w") as fp:
        json.dump(S, fp, indent=1, default=float)
    write_md(S, os.path.join(a.out, "summary.md"))


def write_md(S, path):
    f = lambda x: f"{x:.2f}"
    L = ["## Baseline: K=3, MP_mod, SNR 20 dB, v~U[-10,10]\n",
         "| method | P50 err [m] | P90 err [m] | >1 m [%] | P90 dev vs genie [m] | dev>0.5 m [%] | "
         "P90 phase RMS [rad] | P90 peak loss [dB] | P90 v err [m/s] |", "|---" * 9 + "|"]
    for m in METHODS:
        t = S["baseline"][m]
        L.append(f"| {LABELS[m]} | " + " | ".join(f(t[k]) for k in
                 ["p50", "p90", "out1m", "dev90", "dev_out", "ph90", "loss90", "v90"]) + " |")
    for name in ["K", "speed", "snr"]:
        L += [f"\n## sweep {name}: P90 range error [m] / P90 dev vs genie [m] / P90 v err [m/s]\n",
              f"| {name} | " + " | ".join(LABELS[m] for m in METHODS) + " |", "|---" * 5 + "|"]
        for x, t in zip(S[name]["x"], S[name]["table"]):
            L.append(f"| {x} | " + " | ".join(f"{f(t[m]['p90'])} / {f(t[m]['dev90'])} / {f(t[m]['v90'])}"
                                              for m in METHODS) + " |")
    L += ["\n## scenarios (K=3 unless noted, SNR 20 dB): P90 range error / P90 dev vs genie / P90 v err\n",
          "| scenario | " + " | ".join(LABELS[m] for m in METHODS) + " |", "|---" * 5 + "|"]
    for name, t in S["scenarios"].items():
        L.append(f"| {name} | " + " | ".join(f"{f(t[m]['p90'])} / {f(t[m]['dev90'])} / {f(t[m]['v90'])}"
                                             for m in METHODS) + " |")
    L += ["\n## Sample log proc 7121 (uncalibrated distance)\n", "| method | distance [m] | speed [m/s] |",
          "|---|---|---|"]
    for k, t in S["real_log"].items():
        L.append(f"| {k} | {t['d']:.3f} | {t['v']} |")
    with open(path, "w") as fp:
        fp.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
