"""Build the interactive algorithm explorer (docs/jdps_explorer.html).

    python -m cs_agc.build_viz                  # writes docs/jdps_explorer.html
    python -m cs_agc.build_viz --bare --out x  # same page without <html>/<head>/<body>
    python -m cs_agc.build_viz --cases out.json # also exports Python reference results for viz/test_core.js

The page runs viz/jdps_core.js (a JS port of algorithms.jdps / sim_model.generate / ranging)
in the browser; the sample log procedure is embedded as the real-data example.
"""
import argparse
import json

import numpy as np

from .algorithms import JdpsCfg, jdps
from .log_replay import build_measurement, parse_log
from .ranging import estimate_distance
from .sim_model import SCENARIOS, SimCfg, generate

TEMPLATE = "viz/explorer.template.html"
CORE = "viz/jdps_core.js"
OUT = "docs/jdps_explorer.html"
LOG_TIMING = dict(t_meas=565e-6, t_gap=150e-6, t_step=715e-6, se_gap=40e-3, mode0=483e-6)


def meas_json(ms, digits=7):
    r = lambda a: [float(f"{x:.{digits}g}") for x in np.ravel(a)]
    return dict(re=[r(row.real) for row in ms.y], im=[r(row.imag) for row in ms.y],
                f=[float(x) for x in ms.f], ch=[int(c) for c in ms.ch], t=[float(x) for x in ms.t],
                se=[int(s) for s in ms.se], tRef=float(ms.t_ref))


def sample_log_dataset():
    proc = parse_log("data/sample_log.txt")[0]
    ms = build_measurement(proc, [25, 25, 22], 1, LOG_TIMING)
    d = meas_json(ms)
    d.update(name=f"sample_log proc {proc['ts']}", note="data/sample_log.txt, mode0 483 us, 715 us steps, 40 ms subevent gap")
    return d


def export_cases(path):
    rng = np.random.default_rng(7)
    cases = []
    sims = [(SimCfg(n_se=K, snr_db=snr, channel=SCENARIOS[sc]), v)
            for K, snr, sc, v in [(1, 20, "MP_mod", 4.0), (2, 15, "LOS", -7.5), (3, 20, "MP_mod", 3.2),
                                  (3, 8, "MP_severe", -9.6), (4, 25, "MP_mod", 1.1)]]
    items = [(f"sim{i}", generate(rng, cfg, v=v)) for i, (cfg, v) in enumerate(sims)]
    proc = parse_log("data/sample_log.txt")[0]
    items.append(("log", build_measurement(proc, [25, 25, 22], 1, LOG_TIMING)))
    for orders in [(1,), (1, 2)]:
        for name, ms in items:
            y_c, info = jdps(ms, JdpsCfg(orders=orders))
            d = meas_json(ms, 17)
            d.update(name=f"{name}_o{len(orders)}", orders=list(orders),
                     exp=dict(vHat=float(info["v_hat"]), score=float(info["v_score"]), vValid=bool(info["v_valid"]),
                              vUsed=float(info["v_used"]), psi=[float(p) for p in np.ravel(info["psi"])],
                              metric=[float(m) for m in info["metric"]],
                              dBefore=float(estimate_distance(ms.y, ms.ch)), dAfter=float(estimate_distance(y_c, ms.ch))))
            cases.append(d)
    with open(path, "w") as fp:
        json.dump(cases, fp)
    print("wrote", path, len(cases), "cases")


def build(out, bare=False):
    page = open(TEMPLATE, encoding="utf-8").read()
    core = open(CORE, encoding="utf-8").read()
    data = json.dumps(sample_log_dataset(), separators=(",", ":"))
    page = page.replace("/*__JDPS_CORE__*/", core).replace("/*__LOG_DATA__*/null", data)
    if not bare:   # standalone file; --bare leaves out the document skeleton (for hosts that add their own)
        head, body = page.split("</style>", 1)
        page = ('<!doctype html>\n<html lang="zh-CN">\n<head>\n<meta charset="utf-8">\n'
                '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
                + head + "</style>\n</head>\n<body>" + body + "</body>\n</html>\n")
    with open(out, "w", encoding="utf-8") as fp:
        fp.write(page)
    print("wrote", out, f"{len(page) / 1024:.0f} KiB")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--bare", action="store_true", help="omit <html>/<head>/<body>")
    ap.add_argument("--cases", help="export Python reference cases for viz/test_core.js")
    a = ap.parse_args()
    if a.cases:
        export_cases(a.cases)
    else:
        build(a.out, a.bare)


if __name__ == "__main__":
    main()
