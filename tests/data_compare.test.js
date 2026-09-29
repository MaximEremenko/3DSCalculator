"use strict";

// Calculation-versus-data comparison (js/data_compare.js) and the sigma
// column of intensity lists (js/intensity_list_reader.js).
const test = require("node:test");
const assert = require("node:assert/strict");
const { loadScripts } = require("./helpers/browser_env");

const ctx = loadScripts(["js/intensity_list_reader.js", "js/data_compare.js"]);
const C = ctx.DataCompare, R = ctx.IntensityListReader;

const A = 4.4344;
const Bq = [[2 * Math.PI / A, 0, 0], [0, 2 * Math.PI / A, 0], [0, 0, 2 * Math.PI / A]];
const qOf = (h, k, l) => (2 * Math.PI / A) * Math.hypot(h, k, l);

function axis(min, step, n) {
  return Array.from({ length: n }, (_, i) => Math.round((min + i * step) * 1e9) / 1e9);
}

// A grid { h, k, l, shape, I } with I = f(h, k, l).
function grid(h, k, l, f) {
  const I = new Float64Array(h.length * k.length * l.length);
  let p = 0;
  for (const a of h) for (const b of k) for (const c of l) I[p++] = f(a, b, c);
  return { h, k, l, shape: [h.length, k.length, l.length], I };
}

const calcValue = (h, k, l) => 1 + Math.cos(Math.PI * h) * Math.cos(Math.PI * k) + 0.5 * Math.sin(2 * h + k - l) ** 2;
const H = axis(-1, 0.25, 9), K = axis(-1, 0.25, 9), L = axis(-0.5, 0.25, 5);

test("recovers scale, flat and linear background exactly; Q = 0 is left out", () => {
  const calc = grid(H, K, L, calcValue);
  const data = grid(H, K, L, (h, k, l) => 2.5 * calcValue(h, k, l) + 0.7 + 0.3 * qOf(h, k, l));
  data.I[((4 * 9 + 4) * 5) + 2] = 1e6; // (0 0 0): would ruin the fit if used
  const match = C.matchGrids(data, calc);
  assert.equal(match.same, true);
  const f = C.fit(data, calc, match, { scale: true, background: "linear", Bq });
  assert.ok(Math.abs(f.scale - 2.5) < 1e-10 && Math.abs(f.background - 0.7) < 1e-9 && Math.abs(f.linear - 0.3) < 1e-10, JSON.stringify(f));
  assert.equal(f.n, data.I.length - 1);
  assert.equal(f.skipped.origin, 1);
  assert.equal(f.parameters, 3);
  assert.ok(f.chi2 < 1e-15 && f.rwp < 1e-6 && f.r < 1e-6);
  assert.equal(f.weights, "equal");
});

test("flat background and scale agree with Spinteract's closed form (weights 1/sigma^2)", () => {
  const calc = grid(H, K, L, calcValue);
  const data = grid(H, K, L, (h, k, l) => 3 * calcValue(h, k, l) + 1.1 + 0.4 * Math.sin(7 * h + 3 * k + 5 * l));
  data.sigma = data.I.map((v, i) => 0.05 + 0.02 * (i % 7));
  const match = C.matchGrids(data, calc);
  const f = C.fit(data, calc, match, { background: "flat", Bq });
  // spinteract.f90 calc_chi_sq with refine_scale and refine_flat_bgr
  let sum_ss = 0, sum_fs = 0, sum_s = 0, sum_e = 0, sum_f = 0, wyy = 0, n = 0;
  const pts = [];
  for (let i = 0; i < 9; i++) for (let j = 0; j < 9; j++) for (let m = 0; m < 5; m++) {
    if (H[i] === 0 && K[j] === 0 && L[m] === 0) continue;
    const at = (i * 9 + j) * 5 + m, e = 1 / data.sigma[at] ** 2, x = calc.I[at], y = data.I[at];
    sum_ss += e * x * x; sum_fs += e * x * y; sum_s += e * x; sum_e += e; sum_f += e * y; wyy += e * y * y; n++;
    pts.push([e, x, y]);
  }
  const scale = (sum_e * sum_fs - sum_s * sum_f) / (sum_e * sum_ss - sum_s ** 2);
  const flat = (sum_f - sum_s * scale) / sum_e;
  const chi2 = pts.reduce((a, [e, x, y]) => a + e * (scale * x + flat - y) ** 2, 0);
  assert.equal(f.n, n);
  assert.ok(Math.abs(f.scale / scale - 1) < 1e-12 && Math.abs(f.background / flat - 1) < 1e-12, `${f.scale} ${scale} ${f.background} ${flat}`);
  assert.ok(Math.abs(f.chi2 / chi2 - 1) < 1e-9);
  assert.ok(Math.abs(f.rwp - 100 * Math.sqrt(chi2 / wyy)) < 1e-9);
  assert.ok(Math.abs(f.chi2nu - chi2 / (n - 2)) < 1e-9);
  assert.equal(f.weights, "sigma");
});

test("held scale and no background; masked, sigma, calculation and outside points are left out", () => {
  const calc = grid(H, K, L.slice(0, 4), calcValue);
  const data = grid(H, K, L, (h, k, l) => calcValue(h, k, l) + 0.25);
  data.sigma = new Float64Array(data.I.length).fill(0.1);
  data.I[0] = NaN;                  // masked
  data.sigma[1] = 0;                // no uncertainty: left out, as Spinteract does
  data.sigma[2] = NaN;
  calc.I[(1 * 9 + 1) * 4 + 1] = NaN; // calculation not finite at data point (1,1,1)
  const match = C.matchGrids(data, calc);
  assert.equal(match.same, false);
  assert.equal(match.shared, 9 * 9 * 4);
  const held = C.fit(data, calc, match, { scale: false, background: "flat", Bq });
  assert.equal(held.scale, 1);
  assert.ok(Math.abs(held.background - 0.25) < 1e-12);
  assert.deepEqual({ ...held.skipped }, { masked: 1, outside: 81, sigma: 2, calc: 1, origin: 1 });
  assert.equal(held.n, 9 * 9 * 5 - 1 - 81 - 2 - 1 - 1);
  const none = C.fit(data, calc, match, { scale: false, background: "none", Bq });
  assert.equal(none.parameters, 0);
  assert.ok(Math.abs(none.chi2 - none.n * (0.25 / 0.1) ** 2) < 1e-6);
});

test("a data grid that is a sub-grid of the calculation, model and difference grids", () => {
  const calc = grid(axis(-1, 0.125, 17), axis(-1, 0.125, 17), axis(-0.5, 0.125, 9), calcValue);
  const data = grid(H, K, L, (h, k, l) => 2 * calcValue(h, k, l) + 0.5 + (h === 0.5 && k === 0.25 && l === 0 ? 0.3 : 0));
  const match = C.matchGrids(data, calc);
  assert.deepEqual(Array.from(match.ih), H.map((_, i) => 2 * i));
  assert.equal(match.shared, data.I.length);
  const f = C.fit(data, calc, match, { background: "flat", Bq });
  assert.ok(Math.abs(f.scale - 2) < 1e-3 && Math.abs(f.background - 0.5) < 1e-2);
  const model = C.modelGrid(calc, f, Bq), diff = C.differenceGrid(data, calc, match, f, Bq);
  assert.equal(model.length, calc.I.length);
  const at = (i, j, m) => (i * 9 + j) * 5 + m;
  assert.ok(Number.isNaN(diff[at(4, 4, 2)]), "Q = 0 is blank");
  const bump = diff[at(6, 5, 2)];
  for (let i = 0; i < 9; i++) for (let j = 0; j < 9; j++) for (let m = 0; m < 5; m++) {
    if (i === 4 && j === 4 && m === 2) continue;
    const ci = 2 * i, cj = 2 * j, cm = 2 * m, cAt = (ci * 17 + cj) * 9 + cm;
    assert.ok(Math.abs(diff[at(i, j, m)] - (data.I[at(i, j, m)] - model[cAt])) < 1e-12);
    if (i !== 6 || j !== 5 || m !== 2) assert.ok(Math.abs(diff[at(i, j, m)]) < Math.abs(bump));
  }
  const lv = C.symmetricLevels(diff, 1);
  assert.ok(Math.abs(lv.max - Math.abs(bump)) < 1e-12);
});

test("errors: nothing shared, and a calculation that does not vary", () => {
  const calc = grid(axis(5, 1, 3), K, L, calcValue);
  const data = grid(H, K, L, calcValue);
  assert.match(C.fit(data, calc, C.matchGrids(data, calc), { Bq }).error, /no data point/);
  const flat = grid(H, K, L, () => 2);
  const f = C.fit(data, flat, C.matchGrids(data, flat), { background: "flat", Bq });
  assert.match(f.error, /does not vary/);
});

test("intensity lists: sigma from a fifth column", () => {
  const rows = (fn) => {
    const out = [];
    for (const h of [0, 0.5]) for (const k of [0, 0.5]) for (const l of [0, 1]) out.push(fn(h, k, l));
    return out.join("\n") + "\n";
  };
  const varying = R.parseText(rows((h, k, l) => `${h} ${k} ${l} ${10 + h + k + l} ${0.1 + h}`));
  assert.ok(varying.sigma && Math.abs(varying.sigma[0] - 0.1) < 1e-7 && Math.abs(varying.sigma[7] - 0.6) < 1e-7);
  assert.equal(varying.sigmaConstant, null);
  const constant = R.parseText(rows((h, k, l) => `${h} ${k} ${l} ${10 + h} 1.00000000`));
  assert.equal(constant.sigma, null);
  assert.equal(constant.sigmaConstant, 1);
  const none = R.parseText(rows((h, k, l) => `${h} ${k} ${l} ${10 + h}`));
  assert.equal(none.sigma, null);
  assert.equal(none.sigmaConstant, null);
  const six = R.parseText(rows((h, k, l) => `${h} ${k} ${l} ${10 + h} 0.5 7`));
  assert.equal(six.sigma, null);
  assert.equal(six.I[0], 10);
  // duplicates average to sigma sqrt(sum sigma^2) / n; sigma <= 0 marks the point
  const dup = R.parseText(rows((h, k, l) => `${h} ${k} ${l} 5 0.3`) + "0 0 0 7 0.4\n0.5 0.5 1 5 0\n");
  assert.ok(Math.abs(dup.I[0] - 6) < 1e-12);
  assert.ok(Math.abs(dup.sigma[0] - 0.25) < 1e-7);
  assert.ok(Number.isNaN(dup.sigma[7]));
  assert.equal(dup.sigmaInvalid, 1);
});
