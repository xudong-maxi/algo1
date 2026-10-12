/*
 * JDPS (Joint Doppler & Phase Stitching) 浏览器端实现 + 仿真模型，逐步保留中间结果用于可视化。
 * 与 cs_agc/algorithms.py::jdps、cs_agc/sim_model.py::generate、cs_agc/ranging.py 一一对应，
 * 由 viz/test_core.js 用 Python 导出的用例逐项比对。
 */
(function (root) {
  "use strict";
  const C = 299792458.0, F0 = 2402e6, DF = 1e6, F_C = 2440e6;
  const K_DOP = 4 * Math.PI / C;
  const N_FFT = 1024, N_BINS = 79, FIRST_PATH_DB = 6.0, D_MAX = 60.0;

  const wrap = (p) => Math.atan2(Math.sin(p), Math.cos(p));

  // ------------------------------------------------------------------ helpers
  function adjacentPairs(ch, orders) {
    const pos = new Map();
    ch.forEach((c, i) => pos.set(c, i));
    const n = [], m = [], o = [];
    for (const od of orders) {
      for (const [c, i] of pos) {
        const j = pos.get(c + od);
        if (j !== undefined) { n.push(i); m.push(j); o.push(od); }
      }
    }
    return { n, m, o };
  }

  function parabolicPeak(grid, metric) {
    let k = 0;
    for (let i = 1; i < metric.length; i++) if (metric[i] > metric[k]) k = i;
    if (k > 0 && k < grid.length - 1) {
      const a = metric[k - 1], b = metric[k], c = metric[k + 1];
      const den = a - 2 * b + c;
      if (den < 0) return grid[k] + 0.5 * (a - c) / den * (grid[1] - grid[0]);
    }
    return grid[k];
  }

  function vGrid(cfg) {
    const n = Math.round((cfg.vMax + cfg.vGuard) / cfg.vStep);
    const g = [];
    for (let i = -n; i <= n; i++) g.push(i * cfg.vStep);
    return g;
  }

  const defaultCfg = () => ({ vMax: 10, vStep: 0.25, vGuard: 1, minScore: 5, orders: [1, 2], nPowerIter: 8, nRefineIter: 3 });

  // y: {re: [A][N], im: [A][N]}; returns a new object multiplied per sample by exp(j*ph[a][n]) or ph[n]
  function rotate(y, phFn) {
    const A = y.re.length, N = y.re[0].length;
    const re = [], im = [];
    for (let a = 0; a < A; a++) {
      const r = new Float64Array(N), i = new Float64Array(N);
      for (let k = 0; k < N; k++) {
        const p = phFn(a, k), c = Math.cos(p), s = Math.sin(p);
        r[k] = y.re[a][k] * c - y.im[a][k] * s;
        i[k] = y.re[a][k] * s + y.im[a][k] * c;
      }
      re.push(r); im.push(i);
    }
    return { re, im };
  }

  function pairProducts(y, pr) {
    const A = y.re.length, P = pr.n.length, zr = [], zi = [];
    for (let a = 0; a < A; a++) {
      const r = new Float64Array(P), i = new Float64Array(P);
      for (let p = 0; p < P; p++) {
        const n = pr.n[p], m = pr.m[p];
        const ar = y.re[a][m], ai = y.im[a][m], br = y.re[a][n], bi = y.im[a][n];
        r[p] = ar * br + ai * bi;          // y[m] * conj(y[n])
        i[p] = ai * br - ar * bi;
      }
      zr.push(r); zi.push(i);
    }
    return { re: zr, im: zi };
  }

  // ------------------------------------------------------------------ JDPS
  /**
   * meas: {y:{re,im}, f[N], ch[N], t[N], se[N], tRef}
   * returns every intermediate result of the pipeline
   */
  function jdps(meas, cfgIn) {
    const cfg = Object.assign(defaultCfg(), cfgIn || {});
    const { y, f, ch, se } = meas;
    const A = y.re.length, N = f.length;
    const tau = meas.t.map((t) => t - meas.tRef);
    let K = 0;
    se.forEach((s) => { if (s + 1 > K) K = s + 1; });
    const orders = cfg.orders, O = orders.length;
    const pr = adjacentPairs(ch, orders);
    const P = pr.n.length;
    const G = O * K * K;
    const g = pr.n.map((n, p) => (orders.indexOf(pr.o[p]) * K + se[n]) * K + se[pr.m[p]]);

    // ---- stage 1: speed search
    const z = pairProducts(y, pr);
    const seMid = [];
    for (let k = 0; k < K; k++) {
      let s = 0, c = 0;
      for (let i = 0; i < N; i++) if (se[i] === k) { s += tau[i]; c++; }
      seMid.push(c ? s / c : 0);
    }
    const alphaRed = pr.n.map((n, p) => {
      const m = pr.m[p];
      const alpha = K_DOP * (f[m] * tau[m] - f[n] * tau[n]);
      return alpha - K_DOP * F_C * (seMid[se[m]] - seMid[se[n]]);
    });
    const grid = vGrid(cfg);
    const metric = new Float64Array(grid.length);
    const groupMetric = Array.from({ length: G }, () => new Float64Array(grid.length));
    const accR = new Float64Array(G), accI = new Float64Array(G);
    for (let vi = 0; vi < grid.length; vi++) {
      const v = grid[vi];
      for (let a = 0; a < A; a++) {
        accR.fill(0); accI.fill(0);
        for (let p = 0; p < P; p++) {
          const ph = v * alphaRed[p], c = Math.cos(ph), s = Math.sin(ph);
          accR[g[p]] += c * z.re[a][p] - s * z.im[a][p];
          accI[g[p]] += c * z.im[a][p] + s * z.re[a][p];
        }
        for (let q = 0; q < G; q++) {
          const mag = Math.hypot(accR[q], accI[q]);
          metric[vi] += mag;
          groupMetric[q][vi] += mag;
        }
      }
    }
    const vHat = parabolicPeak(grid, metric);
    let kPeak = 0;
    for (let i = 1; i < metric.length; i++) if (metric[i] > metric[kPeak]) kPeak = i;

    const groupPairs = new Array(G).fill(0);
    g.forEach((q) => groupPairs[q]++);
    let mu = 0, var_ = 0;
    for (let a = 0; a < A; a++) {
      const S = new Float64Array(G);
      for (let p = 0; p < P; p++) S[g[p]] += z.re[a][p] ** 2 + z.im[a][p] ** 2;
      for (let q = 0; q < G; q++) { mu += Math.sqrt(Math.PI / 4 * S[q]); var_ += (1 - Math.PI / 4) * S[q]; }
    }
    const sigma = Math.sqrt(var_);
    const base = { cfg, K, A, N, O, G, tau, pairs: pr, g, groupPairs, z, alphaRed, seMid, grid, metric, groupMetric, vHat, kPeak, mu, sigma };
    if (P === 0 || sigma === 0) {
      return Object.assign(base, { status: "no_signal", score: 0, vValid: false, atEdge: false, vUsed: 0, y1: y, y2: y, psi: new Array(K).fill(0), trace: [] });
    }
    const score = (metric[kPeak] - mu) / sigma;
    const atEdge = kPeak === 0 || kPeak === grid.length - 1;
    const vValid = score >= cfg.minScore && !atEdge;
    const vUsed = vValid ? vHat : 0;

    // ---- stage 2: Doppler + range migration compensation
    const y1 = rotate(y, (a, n) => K_DOP * f[n] * vUsed * tau[n]);

    // ---- stage 3: subevent phase
    const z1 = pairProducts(y1, pr);
    // Z[a][oi][k][j] as flat {re,im} arrays of length A*G
    const Zr = new Float64Array(A * G), Zi = new Float64Array(A * G);
    for (let a = 0; a < A; a++) for (let p = 0; p < P; p++) { Zr[a * G + g[p]] += z1.re[a][p]; Zi[a * G + g[p]] += z1.im[a][p]; }
    const sync = phaseSync(Zr, Zi, A, O, K, cfg.nPowerIter, cfg.nRefineIter);
    const psi = sync.psi;
    const y2 = rotate(y1, (a, n) => -psi[se[n]]);
    return Object.assign(base, { status: "ok", score, vValid, atEdge, vUsed, y1, z1, Zr, Zi, y2, psi, c: sync.c, cInit: sync.cInit, trace: sync.trace });
  }

  function phaseSync(Zr, Zi, A, O, K, nPower, nRefine) {
    const G = O * K * K, idx = (a, o, k, j) => a * G + (o * K + k) * K + j;
    const c = [];
    for (let a = 0; a < A; a++) {
      const row = [];
      for (let o = 0; o < O; o++) {
        let r = 0, i = 0;
        for (let k = 0; k < K; k++) { r += Zr[idx(a, o, k, k)]; i += Zi[idx(a, o, k, k)]; }
        row.push(Math.atan2(i, r));
      }
      c.push(row);
    }
    const cInit = c.map((r) => r.slice());
    let psi = new Array(K).fill(0);
    const trace = [];
    if (K === 1) return { psi, c, cInit, trace };
    for (let it = 0; it < nRefine; it++) {
      const Qr = Array.from({ length: K }, () => new Float64Array(K)), Qi = Array.from({ length: K }, () => new Float64Array(K));
      for (let a = 0; a < A; a++) for (let o = 0; o < O; o++) {
        const cc = Math.cos(-c[a][o]), ss = Math.sin(-c[a][o]);
        for (let k = 0; k < K; k++) for (let j = 0; j < K; j++) {
          const zr = Zr[idx(a, o, k, j)], zi = Zi[idx(a, o, k, j)];
          Qr[k][j] += zr * cc - zi * ss; Qi[k][j] += zr * ss + zi * cc;
        }
      }
      const Hr = [], Hi = [];
      for (let k = 0; k < K; k++) {
        Hr.push(new Float64Array(K)); Hi.push(new Float64Array(K));
        for (let j = 0; j < K; j++) { Hr[k][j] = Qr[k][j] + Qr[j][k]; Hi[k][j] = Qi[k][j] - Qi[j][k]; }
      }
      let shift = 0;
      for (let k = 0; k < K; k++) { let s = 0; for (let j = 0; j < K; j++) s += Math.hypot(Hr[k][j], Hi[k][j]); shift = Math.max(shift, s); }
      const H = { re: Hr.map((r) => Array.from(r)), im: Hi.map((r) => Array.from(r)) };
      for (let k = 0; k < K; k++) { Hr[k][k] = shift; Hi[k][k] = 0; }
      let xr = psi.map((p) => Math.cos(-p)), xi = psi.map((p) => Math.sin(-p));
      const iters = [];
      for (let t = 0; t < nPower; t++) {
        const nr = new Array(K).fill(0), ni = new Array(K).fill(0);
        for (let k = 0; k < K; k++) for (let j = 0; j < K; j++) {
          nr[k] += Hr[k][j] * xr[j] - Hi[k][j] * xi[j];
          ni[k] += Hr[k][j] * xi[j] + Hi[k][j] * xr[j];
        }
        for (let k = 0; k < K; k++) { const mag = Math.hypot(nr[k], ni[k]) + 1e-30; nr[k] /= mag; ni[k] /= mag; }
        xr = nr; xi = ni;
        const p = xr.map((r, k) => -Math.atan2(xi[k], r));
        iters.push(p.map((v) => wrap(v - p[0])));
      }
      psi = xr.map((r, k) => -Math.atan2(xi[k], r));
      psi = psi.map((v) => wrap(v - psi[0]));
      for (let a = 0; a < A; a++) for (let o = 0; o < O; o++) {
        let r = 0, i = 0;
        for (let k = 0; k < K; k++) for (let j = 0; j < K; j++) {
          const ph = -(psi[j] - psi[k]), cc = Math.cos(ph), ss = Math.sin(ph);
          const zr = Zr[idx(a, o, k, j)], zi = Zi[idx(a, o, k, j)];
          r += zr * cc - zi * ss; i += zr * ss + zi * cc;
        }
        c[a][o] = Math.atan2(i, r);
      }
      trace.push({ H, shift, iters, psi: psi.slice(), c: c.map((r) => r.slice()) });
    }
    return { psi, c, cInit, trace };
  }

  // ------------------------------------------------------------------ IFFT ranging
  function fftInPlace(re, im, inverse) {
    const n = re.length;
    for (let i = 1, j = 0; i < n; i++) {
      let bit = n >> 1;
      for (; j & bit; bit >>= 1) j ^= bit;
      j ^= bit;
      if (i < j) { [re[i], re[j]] = [re[j], re[i]]; [im[i], im[j]] = [im[j], im[i]]; }
    }
    for (let len = 2; len <= n; len <<= 1) {
      const ang = (inverse ? 2 : -2) * Math.PI / len, wr = Math.cos(ang), wi = Math.sin(ang);
      for (let i = 0; i < n; i += len) {
        let cr = 1, ci = 0;
        for (let j = 0; j < len / 2; j++) {
          const ur = re[i + j], ui = im[i + j];
          const vr = re[i + j + len / 2] * cr - im[i + j + len / 2] * ci;
          const vi = re[i + j + len / 2] * ci + im[i + j + len / 2] * cr;
          re[i + j] = ur + vr; im[i + j] = ui + vi;
          re[i + j + len / 2] = ur - vr; im[i + j + len / 2] = ui - vi;
          const t = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = t;
        }
      }
    }
  }

  function delayProfile(y, ch, paths) {
    const prof = new Float64Array(N_FFT);
    const use = paths || y.re.map((_, a) => a);
    for (const a of use) {
      const re = new Float64Array(N_FFT), im = new Float64Array(N_FFT);
      ch.forEach((c, n) => {
        const w = 0.5 - 0.5 * Math.cos(2 * Math.PI * (c + 1) / (N_BINS + 1));   // np.hanning(81)[1:-1]
        re[c] = y.re[a][n] * w; im[c] = y.im[a][n] * w;
      });
      fftInPlace(re, im, true);
      for (let k = 0; k < N_FFT; k++) prof[k] += (re[k] ** 2 + im[k] ** 2) / (N_FFT * N_FFT);
    }
    const dist = Array.from({ length: N_FFT }, (_, k) => k * C / (2 * N_FFT * DF));
    return { dist, prof };
  }

  function estimateDistance(y, ch, paths) {
    const { prof } = delayProfile(y, ch, paths);
    const binM = C / (2 * N_FFT * DF);
    const nMax = Math.ceil(D_MAX / binM);
    const p = prof.slice(0, nMax);
    let pk = 0;
    for (let i = 1; i < nMax; i++) if (p[i] > p[pk]) pk = i;
    const thr = p[pk] * Math.pow(10, -FIRST_PATH_DB / 10);
    let k = pk;
    for (let i = 1; i < nMax - 1; i++) if (p[i] >= thr && p[i] >= p[i - 1] && p[i] >= p[i + 1]) { k = i; break; }
    if (k > 0 && k < nMax - 1) {
      const a = Math.log(p[k - 1] + 1e-30), b = Math.log(p[k] + 1e-30), c = Math.log(p[k + 1] + 1e-30);
      const den = a - 2 * b + c;
      k = k + (den !== 0 ? 0.5 * (a - c) / den : 0);
    }
    return k * binM;
  }

  // ------------------------------------------------------------------ simulation
  function rng(seed) {
    let s = seed >>> 0;
    const u = () => { s = (s + 0x6D2B79F5) >>> 0; let t = s; t = Math.imul(t ^ (t >>> 15), t | 1); t ^= t + Math.imul(t ^ (t >>> 7), t | 61); return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
    const normal = () => { let a = u(); while (a === 0) a = u(); return Math.sqrt(-2 * Math.log(a)) * Math.cos(2 * Math.PI * u()); };
    return { u, normal, uniform: (lo, hi) => lo + (hi - lo) * u(), exponential: (mean) => -mean * Math.log(1 - u()) };
  }

  const SCENARIOS = {
    LOS: { nNlos: 2, kDb: 10, meanExcess: 3, decay: 6 },
    MP_mod: { nNlos: 4, kDb: 3, meanExcess: 4, decay: 6 },
    MP_severe: { nNlos: 6, kDb: 0, meanExcess: 3, decay: 6 },
  };
  const CS_CHANNELS = [];
  for (let c = 2; c < 77; c++) if (c < 23 || c > 25) CS_CHANNELS.push(c);

  function splitSteps(n, K) {
    const out = [];
    for (let k = 0; k < K; k++) out.push(Math.floor(n / K) + (k < n % K ? 1 : 0));
    return out;
  }

  /** p: {seed, K, A, snrDb, d0, v, scenario, tMeas, tGap, seGap} */
  function simulate(p) {
    const R = rng(p.seed);
    const tMeas = p.tMeas ?? 565e-6, tGap = p.tGap ?? 150e-6, seGap = p.seGap ?? 40e-3, tStep = tMeas + tGap;
    const ch = CS_CHANNELS.slice();
    for (let i = ch.length - 1; i > 0; i--) { const j = Math.floor(R.u() * (i + 1)); [ch[i], ch[j]] = [ch[j], ch[i]]; }
    const counts = splitSteps(ch.length, p.K);
    const se = [], idx = [];
    counts.forEach((c, k) => { for (let i = 0; i < c; i++) { se.push(k); idx.push(i); } });
    const seStart = [0];
    for (let k = 1; k < p.K; k++) seStart.push(seStart[k - 1] + counts[k - 1] * tStep - tGap + seGap);
    const t = se.map((k, n) => seStart[k] + idx[n] * tStep + tMeas / 2);
    const tRef = t.reduce((s, v) => s + v, 0) / t.length;
    const f = ch.map((c) => F0 + c * DF);

    const sc = SCENARIOS[p.scenario] || SCENARIOS.MP_mod;
    const excess = [0];
    for (let i = 0; i < sc.nNlos; i++) excess.push(R.exponential(sc.meanExcess));
    const nl = excess.slice(1).sort((a, b) => a - b);
    const ex = [0].concat(nl);
    let pn = nl.map((e) => Math.exp(-e / sc.decay));
    const tot = pn.reduce((s, v) => s + v, 0) || 1;
    pn = pn.map((v) => v / tot * Math.pow(10, -sc.kDb / 10));
    const amp = [1].concat(pn).map(Math.sqrt);
    const dP = ex.map((e) => p.d0 + e);

    const A = p.A, N = ch.length;
    const thI = [], thR = [];
    for (let k = 0; k < p.K; k++) thI.push(R.uniform(0, 2 * Math.PI));
    for (let k = 0; k < p.K; k++) thR.push(R.uniform(0, 2 * Math.PI));
    const theta = thI.map((v, k) => wrap(v + thR[k]));
    const h = { re: [], im: [] }, h0 = { re: [], im: [] };
    for (let a = 0; a < A; a++) {
      const ph = dP.map(() => R.uniform(0, 2 * Math.PI));
      const hr = new Float64Array(N), hi = new Float64Array(N), h0r = new Float64Array(N), h0i = new Float64Array(N);
      for (let n = 0; n < N; n++) {
        for (let q = 0; q < dP.length; q++) {
          const d = dP[q] + p.v * (t[n] - tRef);
          const pp = ph[q] - 2 * Math.PI * f[n] * d / C, p0 = ph[q] - 2 * Math.PI * f[n] * dP[q] / C;
          hr[n] += amp[q] * Math.cos(pp); hi[n] += amp[q] * Math.sin(pp);
          h0r[n] += amp[q] * Math.cos(p0); h0i[n] += amp[q] * Math.sin(p0);
        }
      }
      h.re.push(hr); h.im.push(hi);
      h0.re.push(h0r.map((r, n) => r * r - h0i[n] * h0i[n])); h0.im.push(h0r.map((r, n) => 2 * r * h0i[n]));
    }
    let pSig = 0;
    for (let a = 0; a < A; a++) for (let n = 0; n < N; n++) pSig += h.re[a][n] ** 2 + h.im[a][n] ** 2;
    pSig /= A * N;
    const sn = Math.sqrt(pSig / Math.pow(10, p.snrDb / 10) / 2);
    const y = { re: [], im: [] };
    for (let a = 0; a < A; a++) {
      const r = new Float64Array(N), i = new Float64Array(N);
      for (let n = 0; n < N; n++) {
        const k = se[n];
        const ci = Math.cos(thI[k]), si = Math.sin(thI[k]), cr = Math.cos(thR[k]), sr = Math.sin(thR[k]);
        const ar = h.re[a][n] * ci - h.im[a][n] * si + sn * R.normal(), ai = h.re[a][n] * si + h.im[a][n] * ci + sn * R.normal();
        const br = h.re[a][n] * cr - h.im[a][n] * sr + sn * R.normal(), bi = h.re[a][n] * sr + h.im[a][n] * cr + sn * R.normal();
        r[n] = ar * br - ai * bi; i[n] = ar * bi + ai * br;
      }
      y.re.push(r); y.im.push(i);
    }
    return { y, f, ch, t, se, tRef, counts, truth: { d0: p.d0, v: p.v, theta, psi: theta.map((v) => wrap(v - theta[0])), hClean: h0, paths: { d: dP, amp } } };
  }

  const api = { C, F0, DF, F_C, K_DOP, N_FFT, wrap, adjacentPairs, parabolicPeak, vGrid, defaultCfg, rotate, jdps, phaseSync, delayProfile, estimateDistance, simulate, SCENARIOS, splitSteps };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.JDPS = api;
})(typeof window !== "undefined" ? window : globalThis);
