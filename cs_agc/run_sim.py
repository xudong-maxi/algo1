"""Monte-Carlo simulation: multi-subevent AGC phase + Doppler compensation.

Usage:  python -m cs_agc.run_sim [--trials 400] [--out results]
Writes figures (PNG) and summary tables (JSON / Markdown) into --out.
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

from .algorithms import JdpsCfg, doppler_phasor, genie, jdps, legacy, raw, adjacent_pairs
from .ranging import delay_profile, estimate_distance, to_grid
from .sim_model import C, SCENARIOS, SimCfg, TimingCfg, generate

plt.rcParams.update({"font.sans-serif": ["DejaVu Sans"], "axes.grid": True, "grid.alpha": 0.3})

METHODS = ["genie", "raw", "agc_uncomp", "single_se", "legacy_per_se", "legacy_global", "jdps"]
LABELS = {
    "genie": "Genie (true v, true AGC)",
    "raw": "No compensation",
    "agc_uncomp": "True Doppler, AGC not handled",
    "single_se": "Legacy, 1st subevent only",
    "legacy_per_se": "Legacy per subevent + stitch",
    "legacy_global": "Legacy over whole procedure",
    "jdps": "Proposed JDPS",
}
COLORS = {
    "genie": "#222222", "raw": "#bbbbbb", "agc_uncomp": "#e377c2", "single_se": "#8c564b",
    "legacy_per_se": "#d62728", "legacy_global": "#ff7f0e", "jdps": "#1f77b4",
}


# --------------------------------------------------------------------------- #
def run_methods(ms, jcfg: JdpsCfg):
    """Returns {method: (y_comp, channel_idx_used, v_hat or nan)}."""
    vg = jcfg.v_grid
    out = {"genie": (genie(ms), ms.ch, ms.v), "raw": (raw(ms), ms.ch, np.nan)}
    out["agc_uncomp"] = (ms.y * doppler_phasor(ms.f, ms.t - ms.t_ref, ms.v), ms.ch, ms.v)
    sel0 = np.where(ms.se == 0)[0]
    y0, _, v0 = legacy(ms, vg, sel0)
    out["single_se"] = (y0, ms.ch[sel0], v0)
    yp = np.zeros_like(ms.y)
    for k in range(ms.se.max() + 1):
        sel = np.where(ms.se == k)[0]
        yp[:, sel] = legacy(ms, vg, sel)[0]
    out["legacy_per_se"] = (yp, ms.ch, np.nan)
    yl, _, vl = legacy(ms, vg)
    out["legacy_global"] = (yl, ms.ch, vl)
    yj, info = jdps(ms, jcfg)
    out["jdps"] = (yj, ms.ch, info["v_used"])
    return out


def phase_residual_rms(yc, yg):
    """RMS of per-tone phase error vs genie after removing one common phase."""
    r = yc * np.conj(yg)
    r = r / (np.abs(r) + 1e-30)
    r = r * np.exp(-1j * np.angle(r.sum()))
    return float(np.sqrt(np.mean(np.angle(r) ** 2)))


def one_trial(args):
    seed, cfg, jcfg, methods = args
    rng = np.random.default_rng(seed)
    ms = generate(rng, cfg)
    res = run_methods(ms, jcfg)
    yg = res["genie"][0]
    d_g = estimate_distance(yg, ms.ch)
    pk_g = delay_profile(to_grid(yg, ms.ch))[1].max()
    rows = {}
    for k in methods:
        yc, ch, v = res[k]
        d = estimate_distance(yc, ch)
        full = yc.shape[1] == yg.shape[1]
        rows[k] = dict(
            err=d - ms.d0,
            dev=abs(d - d_g),
            ph=phase_residual_rms(yc, yg) if full else np.nan,
            loss=10 * np.log10(delay_profile(to_grid(yc, ch))[1].max() / pk_g) if full else np.nan,
            verr=abs(v - ms.v) if np.isfinite(v) else np.nan,
        )
    return rows


def monte_carlo(cfg, n, jcfg=JdpsCfg(), methods=METHODS, seed0=0, pool=None):
    args = [(seed0 + i, cfg, jcfg, methods) for i in range(n)]
    rows = pool.map(one_trial, args, chunksize=8) if pool else list(map(one_trial, args))
    return {k: {m: np.array([r[k][m] for r in rows]) for m in rows[0][k]} for k in methods}


def stats(r):
    ae = np.abs(r["err"])
    f = lambda x, p: float(np.nanpercentile(x, p)) if np.isfinite(x).any() else float("nan")
    return dict(p50=f(ae, 50), p90=f(ae, 90), p95=f(ae, 95), out1m=float(np.mean(ae > 1.0) * 100),
                dev90=f(r["dev"], 90), dev_out=float(np.mean(r["dev"] > 0.5) * 100),
                ph90=f(r["ph"], 90), loss10=f(-r["loss"], 90), v90=f(r["verr"], 90))


# --------------------------------------------------------------------------- #
# figures
# --------------------------------------------------------------------------- #
def fig_phenomenon(out):
    rng = np.random.default_rng(7)
    cfg = SimCfg(n_se=3, snr_db=25, channel=SCENARIOS["LOS"])
    ms = generate(rng, cfg, d0=6.0, v=5.0)
    res = run_methods(ms, JdpsCfg())
    order = np.argsort(ms.f)
    a = 0
    fig, ax = plt.subplots(2, 2, figsize=(13, 8.5))
    # phase vs frequency, distance slope removed so only impairments remain
    ref = np.exp(1j * 4 * np.pi * ms.f * ms.d0 / C)
    panels = [("raw", "(a) Raw: Doppler + AGC"), ("agc_uncomp", "(b) Doppler removed, AGC not handled"),
              ("jdps", "(c) After proposed JDPS")]
    for (k, title), axx in zip(panels, [ax[0, 0], ax[0, 1], ax[1, 0]]):
        y = res[k][0][a] * ref
        ph = np.angle(y * np.conj(np.exp(1j * np.angle((res["genie"][0][a] * ref).sum()))))
        for s in range(3):
            sel = order[ms.se[order] == s]
            axx.plot(ms.f[sel] / 1e6, ph[sel], "o", ms=4, label=f"subevent {s}")
        axx.set_ylim(-np.pi - 0.2, np.pi + 0.2)
        axx.set_title(title)
        axx.set_xlabel("frequency [MHz]"), axx.set_ylabel("phase (distance removed) [rad]")
        axx.legend(fontsize=8, loc="lower right")
    axx = ax[1, 1]
    for k in ["genie", "agc_uncomp", "legacy_per_se", "legacy_global", "jdps"]:
        dist, p = delay_profile(to_grid(res[k][0], ms.ch))
        axx.plot(dist, 10 * np.log10(p / p.max() + 1e-12), color=COLORS[k], label=LABELS[k],
                 lw=2.2 if k in ("jdps", "genie") else 1.2, ls="--" if k == "genie" else "-")
    axx.axvline(ms.d0, color="k", lw=0.8, ls=":")
    axx.set_xlim(0, 30), axx.set_ylim(-30, 1)
    axx.set_xlabel("distance [m]"), axx.set_ylabel("IFFT power [dB, self-normalised]")
    axx.set_title(f"(d) IFFT profile, d0={ms.d0} m, v={ms.v} m/s, K=3")
    axx.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig1_phenomenon.png"), dpi=130)
    plt.close(fig)


def fig_cdf(res, out, title, fname):
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    for k in METHODS:
        e = np.sort(np.abs(res[k]["err"]))
        ax[0].plot(e, np.arange(1, len(e) + 1) / len(e), color=COLORS[k], label=LABELS[k],
                   lw=2.2 if k == "jdps" else 1.3, ls="--" if k == "genie" else "-")
        if k == "genie" or not np.isfinite(res[k]["ph"]).any():
            continue
        e = np.sort(res[k]["dev"])
        ax[1].plot(e, np.arange(1, len(e) + 1) / len(e), color=COLORS[k], label=LABELS[k],
                   lw=2.2 if k == "jdps" else 1.3)
    ax[0].set_xscale("log"), ax[0].set_xlim(0.01, 30)
    ax[0].set_xlabel("|range error| vs truth [m]"), ax[0].set_ylabel("CDF")
    ax[0].set_title(title)
    ax[0].legend(fontsize=8)
    ax[1].set_xscale("log"), ax[1].set_xlim(1e-3, 30)
    ax[1].set_xlabel("|d - d_genie| [m]  (compensation-induced error only)")
    ax[1].set_title("Deviation from genie")
    fig.tight_layout()
    fig.savefig(os.path.join(out, fname), dpi=130)
    plt.close(fig)


def fig_sweep(xs, table, xlabel, out, fname, methods, title):
    keys = [("p90", "P90 |range error| [m]", True), ("dev90", "P90 |d - d_genie| [m]", True),
            ("ph90", "P90 residual phase RMS [rad]", True), ("loss10", "P90 IFFT peak loss [dB]", False)]
    fig, ax = plt.subplots(1, 4, figsize=(18, 4.2))
    for (key, lab, logy), axx in zip(keys, ax):
        for k in methods:
            ys = [table[i][k][key] for i in range(len(xs))]
            if not np.isfinite(ys).any():
                continue
            if logy:
                ys = np.maximum(ys, 1e-3)          # exact zeros (genie) on a log axis
            axx.plot(xs, ys, "o-", color=COLORS[k], label=LABELS[k], lw=2.2 if k == "jdps" else 1.3)
        axx.set_xlabel(xlabel), axx.set_title(lab, fontsize=10)
        if logy:
            axx.set_yscale("log")
            axx.set_ylim(bottom=5e-4)
    ax[0].legend(fontsize=7)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(os.path.join(out, fname), dpi=130)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# complexity model for the MCU
# --------------------------------------------------------------------------- #
CYC = dict(cmul=8, cmac=8, cabs=20, sincos=60, misc_per_pair=10)   # Cortex-M4F/M33 float32, conservative


def complexity(n_se, n_ant=4, jcfg=JdpsCfg(), n_ch=72, f_cpu=128e6, trials=50):
    """Average op counts / cycle estimate for the proposed algorithm."""
    rng = np.random.default_rng(0)
    P = np.mean([len(adjacent_pairs(rng.permutation(np.r_[2:23, 26:77]), jcfg.orders)[0])
                 for _ in range(trials)])
    G = len(jcfg.v_grid)
    Ng = len(jcfg.orders) * n_se * n_se
    ops = {
        "pair products (y_m*conj(y_n))": P * n_ant * CYC["cmul"],
        "pair phase-rate + start/step phasors": P * (2 * CYC["sincos"] + CYC["misc_per_pair"]),
        "search: phasor recursion": P * G * CYC["cmul"],
        "search: group accumulation": P * G * n_ant * CYC["cmac"],
        "search: |Z| per group": G * n_ant * Ng * CYC["cabs"],
        "Doppler+migration compensation": n_ch * (CYC["sincos"] + n_ant * CYC["cmul"]),
        "re-pair + group sums": P * n_ant * (CYC["cmul"] + CYC["cmac"]),
        "phase sync (K x K power iteration)": jcfg.n_refine_iter * (jcfg.n_power_iter * n_se * n_se
                                                                   + Ng * n_ant) * CYC["cmac"],
        "AGC de-rotation": n_ch * (CYC["cmul"] * n_ant) + n_se * CYC["sincos"],
    }
    total = sum(ops.values())
    return dict(pairs=P, grid=G, groups=Ng, ops=ops, cycles=total, ms=total / f_cpu * 1e3)


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=400)
    ap.add_argument("--out", default="results")
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    N = a.trials
    summary = {}
    fig_phenomenon(a.out)

    with Pool(a.workers) as pool:
        # 1) baseline CDF: K=3, moderate multipath, SNR 20, v ~ U[-10, 10]
        base = SimCfg(n_se=3)
        r = monte_carlo(base, N, pool=pool, seed0=1000)
        summary["baseline_K3_MPmod_snr20"] = {k: stats(v) for k, v in r.items()}
        fig_cdf(r, a.out, "K=3, MP_mod, SNR 20 dB, v~U[-10,10] m/s", "fig2_cdf_K3.png")

        show = ["genie", "agc_uncomp", "legacy_per_se", "legacy_global", "jdps"]
        # 2) speed sweep
        speeds = [0, 2, 5, 8, 10]
        tab = [{k: stats(v) for k, v in monte_carlo(dataclasses.replace(base, v_range=(s, s)), N,
                                                      pool=pool, seed0=2000 + 97 * i).items()}
               for i, s in enumerate(speeds)]
        summary["speed_sweep"] = dict(x=speeds, table=tab)
        fig_sweep(speeds, tab, "|v| [m/s]", a.out, "fig3_speed.png", show, "K=3, MP_mod, SNR 20 dB")

        # 3) SNR sweep
        snrs = [5, 10, 15, 20, 25, 30]
        tab = [{k: stats(v) for k, v in monte_carlo(dataclasses.replace(base, snr_db=s), N,
                                                      pool=pool, seed0=3000 + 97 * i).items()}
               for i, s in enumerate(snrs)]
        summary["snr_sweep"] = dict(x=snrs, table=tab)
        fig_sweep(snrs, tab, "SNR per tone [dB]", a.out, "fig4_snr.png", show, "K=3, MP_mod, v~U[-10,10]")

        # 4) number of subevents
        Ks = [1, 2, 3, 4, 5]
        tab = [{k: stats(v) for k, v in monte_carlo(dataclasses.replace(base, n_se=K), N,
                                                      pool=pool, seed0=4000 + 97 * i).items()}
               for i, K in enumerate(Ks)]
        summary["K_sweep"] = dict(x=Ks, table=tab)
        fig_sweep(Ks, tab, "number of subevents K", a.out, "fig5_K.png", show, "MP_mod, SNR 20 dB, v~U[-10,10]")

        # 5) scenario table
        scen = {
            "LOS": dataclasses.replace(base, channel=SCENARIOS["LOS"]),
            "MP_mod": base,
            "MP_severe": dataclasses.replace(base, channel=SCENARIOS["MP_severe"]),
            "MP_mod + per-path velocity": dataclasses.replace(
                base, channel=dataclasses.replace(SCENARIOS["MP_mod"], per_path_velocity=True)),
            "MP_mod + 100 ms SE gap": dataclasses.replace(base, timing=TimingCfg(se_gap=0.1)),
            "MP_mod + AGC per antenna": dataclasses.replace(base, agc_shared_across_ant=False),
        }
        summary["scenarios"] = {}
        for i, (name, cfg) in enumerate(scen.items()):
            r = monte_carlo(cfg, N, pool=pool, seed0=5000 + 97 * i)
            summary["scenarios"][name] = {k: stats(v) for k, v in r.items()}
        # per-antenna AGC handled by the per-antenna option
        r = monte_carlo(scen["MP_mod + AGC per antenna"], N, jcfg=JdpsCfg(per_ant_theta=True),
                        pool=pool, seed0=5000 + 97 * 5)
        summary["scenarios"]["MP_mod + AGC per antenna"]["jdps_per_ant"] = stats(r["jdps"])

        # 6) complexity knobs: v grid step and pair orders
        summary["knobs"] = {}
        for name, jc in {"step0.25_orders12": JdpsCfg(), "step0.5_orders12": JdpsCfg(v_step=0.5),
                         "step1.0_orders12": JdpsCfg(v_step=1.0), "step0.5_orders1": JdpsCfg(v_step=0.5, orders=(1,))}.items():
            r = monte_carlo(dataclasses.replace(base, snr_db=10), N, jcfg=jc, methods=["genie", "jdps"],
                            pool=pool, seed0=6000)
            summary["knobs"][name] = dict(stats=stats(r["jdps"]), cplx=complexity(3, jcfg=jc)["ms"])

    summary["complexity"] = {f"K={K}": complexity(K) for K in [1, 2, 3, 5]}
    with open(os.path.join(a.out, "summary.json"), "w") as fp:
        json.dump(summary, fp, indent=1, default=float)
    write_markdown(summary, os.path.join(a.out, "summary.md"))


def write_markdown(s, path):
    L = []
    fmt = lambda x: "-" if not np.isfinite(x) else (f"{x:.2f}" if abs(x) < 100 else f"{x:.0f}")
    L.append("## Baseline: K=3, MP_mod, SNR 20 dB, v~U[-10,10] m/s\n")
    L.append("| method | P50 err [m] | P90 err [m] | P95 err [m] | >1 m [%] | P90 dev vs genie [m] | dev>0.5 m [%] | P90 phase RMS [rad] | P90 peak loss [dB] | P90 v err [m/s] |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for k, v in s["baseline_K3_MPmod_snr20"].items():
        L.append(f"| {LABELS[k]} | " + " | ".join(fmt(v[x]) for x in
                 ["p50", "p90", "p95", "out1m", "dev90", "dev_out", "ph90", "loss10", "v90"]) + " |")
    for key, xl in [("speed_sweep", "|v| [m/s]"), ("snr_sweep", "SNR [dB]"), ("K_sweep", "K")]:
        sw = s[key]
        L.append(f"\n## {key}: P90 range error [m] / P90 deviation from genie [m]\n")
        ms = ["genie", "legacy_per_se", "legacy_global", "jdps"]
        L.append(f"| {xl} | " + " | ".join(LABELS[m] for m in ms) + " |")
        L.append("|---" * (len(ms) + 1) + "|")
        for x, t in zip(sw["x"], sw["table"]):
            L.append(f"| {x} | " + " | ".join(f"{fmt(t[m]['p90'])} / {fmt(t[m]['dev90'])}" for m in ms) + " |")
    L.append("\n## Scenarios (K=3, SNR 20 dB): P90 range error [m] / P90 dev vs genie [m] / P90 v err [m/s]\n")
    ms = ["genie", "legacy_per_se", "legacy_global", "jdps"]
    L.append("| scenario | " + " | ".join(LABELS[m] for m in ms) + " |")
    L.append("|---" * (len(ms) + 1) + "|")
    for name, t in s["scenarios"].items():
        row = [f"{fmt(t[m]['p90'])} / {fmt(t[m]['dev90'])} / {fmt(t[m]['v90'])}" for m in ms]
        if "jdps_per_ant" in t:
            j = t["jdps_per_ant"]
            row[-1] += f" (per-ant option: {fmt(j['p90'])} / {fmt(j['dev90'])} / {fmt(j['v90'])})"
        L.append(f"| {name} | " + " | ".join(row) + " |")
    L.append("\n## Complexity knobs (K=3, SNR 10 dB)\n")
    L.append("| config | P90 err [m] | P90 dev vs genie [m] | P90 phase RMS [rad] | P90 v err [m/s] | est. time 4 ant [ms] |")
    L.append("|---|---|---|---|---|---|")
    for name, t in s["knobs"].items():
        st = t["stats"]
        L.append(f"| {name} | {fmt(st['p90'])} | {fmt(st['dev90'])} | {fmt(st['ph90'])} | {fmt(st['v90'])} | {t['cplx']:.2f} |")
    L.append("\n## Estimated MCU load (default config, 4 antennas, 128 MHz)\n")
    L.append("| K | pairs | v grid | groups | cycles | time [ms] |")
    L.append("|---|---|---|---|---|---|")
    for k, c in s["complexity"].items():
        L.append(f"| {k[2:]} | {c['pairs']:.0f} | {c['grid']} | {c['groups']} | {c['cycles']:.0f} | {c['ms']:.2f} |")
    c = s["complexity"]["K=3"]
    L.append("\nBreakdown for K=3:\n")
    L.append("| stage | cycles |")
    L.append("|---|---|")
    for k, v in c["ops"].items():
        L.append(f"| {k} | {v:.0f} |")
    with open(path, "w") as fp:
        fp.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
