"""Replay the exact inputs that the firmware hands to subevent_motion_alg.

The firmware prints (see docs/firmware_dump.md for the printf snippet), per procedure:
    [jdps_dump] begin
    [jdps_dump] cfg <ch_num> <subevent_num> <t_mes>
    [jdps_dump] hop <ch_hop_orders[0..ch_num-1]>
    [jdps_dump] per_subevent <ch_num_per_subevent[0..subevent_num-1]>
    [jdps_dump] time <time_per_channel[0..sum(per_subevent)-1]>
    [jdps_dump] path <p> <iq[0].r> <iq[0].i> ... <iq[ALG_CHANNEL_NUM-1].r> <iq[ALG_CHANNEL_NUM-1].i>
    [jdps] v=... (est ..., score ..., valid ...), psi[rad]=...      (firmware result, optional)
    [jdps_dump] end
`iq` is the per-path buffer exactly as passed to subevent_motion_add_speed_path (after calc_iq).
The channel table is built with the same rules as build_channel_table() in c/subevent_motion_alg.c,
then the Python reference (cs_agc.algorithms.jdps) runs on identical inputs.

Usage: python -m cs_agc.replay_dump firmware.log [--pair-order 2]
"""
import argparse
import re
from types import SimpleNamespace

import numpy as np

from .algorithms import JdpsCfg, jdps
from .sim_model import DF, F0

RE_FW = re.compile(r"\[jdps\] v=([-+\d.]+) m/s \(est ([-+\d.]+), score ([-+\d.]+), valid (\d)\), psi\[rad\]=(.*)")


def parse_dumps(path):
    procs, cur = [], None
    with open(path, errors="replace") as fp:
        for line in fp:
            if "[jdps_dump]" in line:
                tok = line.split("[jdps_dump]", 1)[1].split()
                if not tok:
                    continue
                key, vals = tok[0], tok[1:]
                if key == "begin":
                    cur = dict(paths={}, fw=None)
                    procs.append(cur)
                elif cur is None:
                    continue
                elif key == "cfg":
                    cur["ch_num"], cur["subevent_num"], cur["t_mes"] = (int(v) for v in vals[:3])
                elif key in ("hop", "per_subevent", "time"):
                    cur[key] = [int(v) for v in vals]
                elif key == "path":
                    x = np.array([float(v) for v in vals[1:]])
                    cur["paths"][int(vals[0])] = x[0::2] + 1j * x[1::2]
                elif key == "end":
                    cur = None
            elif cur is not None:
                m = RE_FW.search(line)
                if m:
                    cur["fw"] = dict(v=float(m.group(1)), est=float(m.group(2)), score=float(m.group(3)),
                                     valid=int(m.group(4)), psi=[float(p) for p in m.group(5).split()])
    return procs


def channel_table(d):
    """Same rules as build_channel_table(): per-step or log hop layout, start time of each mode-2 step."""
    per = d["per_subevent"]
    step_num = sum(per)
    if d["ch_num"] == step_num:
        per_step = True
    elif d["ch_num"] == step_num - len(per) + 1:
        per_step = False
    else:
        raise ValueError(f"ch_num {d['ch_num']} matches neither layout (steps {step_num})")
    ch, se, t = [], [], []
    seen = set()
    t_us, step, mode2 = 0, 0, 0
    for s, n in enumerate(per):
        for i in range(n):
            if i > 0:
                c = d["hop"][step] if per_step else d["hop"][mode2 + 1]
                mode2 += 1
                if c not in seen:
                    seen.add(c)
                    ch.append(c), se.append(s), t.append(t_us * 1e-6)
            t_us += d["time"][step]
            step += 1
        t_us += d["t_mes"]
    return np.array(ch), np.array(se), np.array(t), per_step


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("log")
    ap.add_argument("--pair-order", type=int, default=2, choices=[1, 2])
    a = ap.parse_args()
    cfg = JdpsCfg(orders=(1,) if a.pair_order == 1 else (1, 2))
    for i, d in enumerate(parse_dumps(a.log)):
        ch, se, t, per_step = channel_table(d)
        paths = sorted(d["paths"])
        y = np.stack([d["paths"][p][ch] for p in paths])
        ms = SimpleNamespace(y=y, ch=ch, f=F0 + ch * DF, t=t, se=se, t_ref=t.mean())
        _, info = jdps(ms, cfg)
        print(f"proc #{i}: hop layout {'per-step' if per_step else 'log'}, {len(ch)} channels, "
              f"{len(paths)} paths, subevents {d['per_subevent']}")
        print(f"  python : v={info['v_used']:+.3f} (est {info['v_hat']:+.3f}, score {info['v_score']:.1f}, "
              f"valid {int(info['v_valid'])}), psi={np.round(np.ravel(info['psi']), 3).tolist()}")
        if d["fw"]:
            fw = d["fw"]
            print(f"  firmware: v={fw['v']:+.3f} (est {fw['est']:+.3f}, score {fw['score']:.1f}, "
                  f"valid {fw['valid']}), psi={fw['psi']}")


if __name__ == "__main__":
    main()
