// 用 Python 参考实现导出的用例检查 jdps_core.js：node viz/test_core.js cases.json
const J = require("./jdps_core.js");
const cases = require(require("path").resolve(process.argv[2]));
let fail = 0;
for (const c of cases) {
  const meas = { y: { re: c.re, im: c.im }, f: c.f, ch: c.ch, t: c.t, se: c.se, tRef: c.tRef };
  const r = J.jdps(meas, { orders: c.orders });
  const e = c.exp;
  const dpsi = Math.max(...r.psi.map((p, k) => Math.abs(J.wrap(p - e.psi[k]))));
  const dmet = Math.max(...e.metric.map((m, i) => Math.abs(m - r.metric[i]) / Math.max(...e.metric)));
  const dB = J.estimateDistance(meas.y, c.ch), dA = J.estimateDistance(r.y2, c.ch);
  const ok = Math.abs(r.vHat - e.vHat) < 1e-6 && Math.abs(r.score - e.score) < 1e-6 && r.vValid === e.vValid &&
    dpsi < 1e-6 && dmet < 1e-9 && Math.abs(dB - e.dBefore) < 1e-6 && Math.abs(dA - e.dAfter) < 1e-6;
  if (!ok) fail++;
  console.log(`${ok ? "PASS" : "FAIL"} ${c.name.padEnd(8)} v ${r.vHat.toFixed(4)}/${e.vHat.toFixed(4)} score ${r.score.toFixed(3)}/${e.score.toFixed(3)} ` +
    `valid ${r.vValid}/${e.vValid} dpsi ${dpsi.toExponential(1)} dmetric ${dmet.toExponential(1)} d ${dA.toFixed(3)}/${e.dAfter.toFixed(3)} (before ${dB.toFixed(3)}/${e.dBefore.toFixed(3)})`);
}
console.log(fail ? `${fail} case(s) failed` : `all ${cases.length} cases match the Python reference`);
process.exit(fail ? 1 : 0);
