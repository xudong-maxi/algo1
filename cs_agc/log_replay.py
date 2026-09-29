"""Replay real CS logs through JDPS and plot phase before / after compensation.

Log format (one CS procedure):
    [cs_local_ch] proc_count = 7121, ch_num = 73
    [cs_local_ch] ch_idx_list: 22,4b,30,...          hop order (hex), 1st = mode-0 step
    key:4,path:P,[iq_data]4,R,7121,[ch]:rssi;i,q;...  P = antenna path, R = 0 local / 1 remote
IQ entries are indexed by channel (0..79), not by hop order.

Usage:
    python -m cs_agc.log_replay LOG [LOG ...] --out log_plots [--split 25,25,22]
"""
import argparse
import csv
import os
import re
from types import SimpleNamespace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .algorithms import JdpsCfg, jdps
from .ranging import delay_profile, estimate_distance, to_grid
from .sim_model import C, CS_CHANNELS, DF, F0

RE_PROC = re.compile(r"proc_count\s*=\s*(\d+),\s*ch_num\s*=\s*(\d+)")
RE_CHLIST = re.compile(r"ch_idx_list:\s*([0-9a-fA-F,\s]+)")
RE_IQ_HDR = re.compile(r"path:(\d+),\[iq_data\](\d+),(\d+),(\d+),")
RE_IQ_ENT = re.compile(r"\[(\d+)\]:(-?\d+);(-?\d+),(-?\d+)")


# --------------------------------------------------------------------------- #
# parsing
# --------------------------------------------------------------------------- #
def parse_log(path):
    """Returns list of procedures: dict(ts, ch_num, hop, iq{(path, role): {ch: complex}})."""
    procs, cur = [], None
    with open(path, errors="replace") as fp:
        for line in fp:
            m = RE_PROC.search(line)
            if m:
                cur = dict(ts=int(m.group(1)), ch_num=int(m.group(2)), hop=None, iq={})
                procs.append(cur)
                continue
            m = RE_CHLIST.search(line)
            if m and cur is not None:
                vals = [int(x, 16) for x in m.group(1).replace(" ", "").strip(",").split(",") if x]
                cur["hop"] = vals[:cur["ch_num"]]
                continue
            m = RE_IQ_HDR.search(line)
            if m:
                ant, role, ts = int(m.group(1)), int(m.group(3)), int(m.group(4))
                tgt = next((p for p in reversed(procs) if p["ts"] == ts), None)
                if tgt is None:
                    continue
                ent = {}
                # corrupted entries (truncated print) simply fail to match and are dropped
                for e in RE_IQ_ENT.finditer(line[m.end() - 1:]):
                    ch, i, q = int(e.group(1)), int(e.group(3)), int(e.group(4))
                    if i or q:
                        ent[ch] = complex(i, q)
                tgt["iq"][(ant, role)] = ent
    return [p for p in procs if p["hop"] and p["iq"]]


def build_measurement(proc, split, n_mode0, tcfg, combine="mul"):
    """Turn one parsed procedure into the (A, N) hop-ordered PBR array JDPS expects."""
    hop = np.array(proc["hop"][n_mode0:])
    if len(hop) != sum(split):
        raise ValueError(f"proc {proc['ts']}: {len(hop)} mode-2 steps but split sums to {sum(split)}")
    ants = sorted({a for a, _ in proc["iq"]})
    y = np.zeros((len(ants), len(hop)), complex)
    missing = []
    for ai, a in enumerate(ants):
        loc, rem = proc["iq"].get((a, 0), {}), proc["iq"].get((a, 1), {})
        for n, ch in enumerate(hop):
            if ch in loc and ch in rem:
                r = rem[ch] if combine == "mul" else np.conj(rem[ch])
                y[ai, n] = loc[ch] * r
            else:
                missing.append((a, int(ch)))       # left as 0 -> ignored by pairs / IFFT
    se = np.repeat(np.arange(len(split)), split)
    idx = np.concatenate([np.arange(c) for c in split])
    se_dur = np.array(split) * tcfg["t_step"] - tcfg["t_gap"]
    se_start = np.concatenate([[0.0], np.cumsum(se_dur + tcfg["se_gap"])[:-1]])
    t = se_start[se] + idx * tcfg["t_step"] + tcfg["t_meas"] / 2
    return SimpleNamespace(y=y, ch=hop, f=F0 + hop * DF, t=t, se=se, t_ref=t.mean(),
                           ants=ants, missing=missing)


# --------------------------------------------------------------------------- #
# plotting
# --------------------------------------------------------------------------- #
SE_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]


def plot_proc(ms, y_c, info, d_before, d_after, title, path):
    A, K = ms.y.shape[0], ms.se.max() + 1
    order = np.argsort(ms.f)
    # remove the strongest-path slope (from the compensated profile) so jumps are visible
    dist, prof = delay_profile(to_grid(y_c, ms.ch))
    d_slope = dist[np.argmax(prof[:len(dist) // 2])]
    ref = np.exp(1j * 4 * np.pi * ms.f * d_slope / C)
    fig, ax = plt.subplots(A, 3, figsize=(17, 3.1 * A + 0.8), squeeze=False)
    for a in range(A):
        for col, (yy, lab) in enumerate([(ms.y[a], "before"), (y_c[a], "after JDPS")]):
            axx = ax[a, col]
            z = yy * ref
            ok = np.abs(z) > 0
            z = z * np.exp(-1j * np.angle(np.sum(z[ok])))  # common rotation for readability
            for k in range(K):
                s = order[(ms.se[order] == k) & ok[order]]
                axx.plot(ms.f[s] / 1e6, np.angle(z[s]), "o-", ms=3.5, lw=0.6, alpha=0.85,
                         color=SE_COLORS[k % len(SE_COLORS)], label=f"SE{k}")
            axx.set_ylim(-np.pi - 0.2, np.pi + 0.2)
            axx.set_title(f"path {ms.ants[a]} – phase {lab} (slope of {d_slope:.2f} m removed)", fontsize=9)
            axx.set_xlabel("freq [MHz]", fontsize=8), axx.set_ylabel("rad", fontsize=8)
            axx.grid(alpha=0.3)
            if a == 0 and col == 0:
                axx.legend(fontsize=7, ncol=K, loc="lower right")
        axx = ax[a, 2]
        for yy, lab, c in [(ms.y[a:a + 1], "before", "#999999"), (y_c[a:a + 1], "after JDPS", "#1f77b4")]:
            dd, pp = delay_profile(to_grid(yy, ms.ch))
            axx.plot(dd, 10 * np.log10(pp / pp.max() + 1e-12), color=c, label=lab)
        axx.set_xlim(0, 40), axx.set_ylim(-30, 1), axx.grid(alpha=0.3)
        axx.set_title(f"path {ms.ants[a]} – IFFT profile", fontsize=9)
        axx.set_xlabel("distance [m] (uncalibrated)", fontsize=8)
        if a == 0:
            axx.legend(fontsize=7)
    psi = np.atleast_1d(info["psi"])
    fig.suptitle(f"{title}   v̂={info['v_hat']:+.2f} m/s   ψ=[{', '.join(f'{p:+.2f}' for p in np.ravel(psi))}] rad   "
                 f"d(before)={d_before:.2f} m  d(after)={d_after:.2f} m   missing={len(ms.missing)}",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=110)
    plt.close(fig)


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("logs", nargs="+")
    ap.add_argument("--out", default="log_plots")
    ap.add_argument("--split", default="25,25,22", help="mode-2 steps per subevent")
    ap.add_argument("--n-mode0", type=int, default=1, help="leading mode-0 entries in ch_idx_list")
    ap.add_argument("--se-gap-ms", type=float, default=40.0)
    ap.add_argument("--t-meas-us", type=float, default=565.0)
    ap.add_argument("--t-gap-us", type=float, default=150.0)
    ap.add_argument("--combine", choices=["mul", "conj"], default="mul",
                    help="PBR combine: local*remote (mul) or local*conj(remote)")
    ap.add_argument("--per-ant", action="store_true", help="estimate AGC phase per antenna path")
    ap.add_argument("--v-max", type=float, default=10.0)
    a = ap.parse_args()

    split = [int(x) for x in a.split.split(",")]
    tcfg = dict(t_meas=a.t_meas_us * 1e-6, t_gap=a.t_gap_us * 1e-6,
                t_step=(a.t_meas_us + a.t_gap_us) * 1e-6, se_gap=a.se_gap_ms * 1e-3)
    jcfg = JdpsCfg(v_max=a.v_max, per_ant_theta=a.per_ant)
    os.makedirs(a.out, exist_ok=True)
    rows = []
    for log in a.logs:
        stem = os.path.splitext(os.path.basename(log))[0]
        for proc in parse_log(log):
            try:
                ms = build_measurement(proc, split, a.n_mode0, tcfg, a.combine)
            except ValueError as e:
                print("skip:", e)
                continue
            y_c, info = jdps(ms, jcfg)
            d_b, d_a = estimate_distance(ms.y, ms.ch), estimate_distance(y_c, ms.ch)
            name = f"{stem}_proc{proc['ts']}"
            plot_proc(ms, y_c, info, d_b, d_a, name, os.path.join(a.out, name + ".png"))
            rows.append(dict(log=stem, proc=proc["ts"], v_hat=round(float(info["v_hat"]), 3),
                             psi=" ".join(f"{p:.3f}" for p in np.ravel(info["psi"])),
                             d_before=round(d_b, 3), d_after=round(d_a, 3), missing=len(ms.missing)))
            print(rows[-1])
    if rows:
        with open(os.path.join(a.out, "summary.csv"), "w", newline="") as fp:
            w = csv.DictWriter(fp, fieldnames=list(rows[0]))
            w.writeheader(), w.writerows(rows)


if __name__ == "__main__":
    main()
