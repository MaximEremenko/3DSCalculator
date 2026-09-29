"use strict";

// Modulated structures against the exact satellite theory of a monatomic
// lattice (commensurate q, so the grid hits the satellites exactly):
//   displacive u = A e sin(2 pi q.R):  |A(G + m q)|^2 = N^2 b^2 J_m(Q.e A)^2
//   proper helix (axis || q):          I_M(G +- q) = N^2 p^2 f^2 A^2 / 2
const test = require("node:test");
const assert = require("node:assert/strict");
const { loadScripts } = require("./helpers/browser_env");

const ctx = loadScripts(["js/rmc6f_reader.js", "js/magnetic_form_factors.js", "js/diffuse_core.js", "js/diffuse_amplitude.js", "js/modulation.js"]);
const core = ctx.DiffuseCore, M = ctx.Modulation;

async function exactSum(spec) {
  const s = spec.sourcesPacked, t = spec.targetsPacked, c = spec.strengths, m = s.length / 3, n = t.length / 3, out = new Float64Array(n * 2);
  for (let k = 0; k < n; k++) {
    let re = 0, im = 0;
    for (let j = 0; j < m; j++) {
      const p = t[k * 3] * s[j * 3] + t[k * 3 + 1] * s[j * 3 + 1] + t[k * 3 + 2] * s[j * 3 + 2];
      re += c[j * 2] * Math.cos(p) - c[j * 2 + 1] * Math.sin(p);
      im += c[j * 2] * Math.sin(p) + c[j * 2 + 1] * Math.cos(p);
    }
    out[k * 2] = re;
    out[k * 2 + 1] = im;
  }
  return { out };
}

// Bessel J_m by its series (small arguments)
function besselJ(m, x) {
  let sum = 0, term = Math.pow(x / 2, m);
  for (let k = 1; k <= m; k++) term /= k;
  for (let k = 0; k < 40; k++) {
    sum += term;
    term *= -(x * x / 4) / ((k + 1) * (k + 1 + m));
  }
  return sum;
}

// simple cubic Fe, a = 3 Angstrom, N x N x N cells
function cubic(N, element = "Fe") {
  const frac = [], cellIndex = [], elements = [];
  for (let a = 0; a < N; a++) for (let b = 0; b < N; b++) for (let c = 0; c < N; c++) {
    frac.push(a / N, b / N, c / N);
    cellIndex.push(a, b, c);
    elements.push(element);
  }
  const s = ctx.RMC6fReader.buildSupercellStructure({ parentCellDeg: [3, 3, 3, 90, 90, 90], supercell: [N, N, N], fractional: Float64Array.from(frac), cellIndex: Int32Array.from(cellIndex), elements });
  return core.attachNeutronCoefficients({ ...s, cellDeg: [3 * N, 3 * N, 3 * N, 90, 90, 90], atoms: N * N * N }, 10);
}

const intensity = (parsed, h, k, l, mode = "nuclear", sub = false) => ctx.DiffuseAmplitude.computeIntensity({
  parsed, h, k, l, Bq: parsed.Bq, backend: "test", sub, deltaOnLattice: true,
  scattering: { type: "neutron", model: "fast" }, runType3: exactSum, magnetic: { mode, qZero: "average" },
});

test("displacive wave: satellites follow N^2 b^2 J_m(Q.e A)^2", async () => {
  const N = 8, A = 0.12;
  const base = cubic(N);
  const mod = M.applyModulation(base, { kind: "displacive", q: [0.25, 0, 0], amplitude: A, direction: [1, 0, 0] });
  assert.deepEqual(Array.from(mod.modulation.commensurate), [true, true, true]);
  const b = base.fca[0], n = base.atoms;
  // G = (1 0 0) with satellites m = -2..2 at h = 1 + m/4
  const hs = [0.5, 0.75, 1, 1.25, 1.5];
  const out = await intensity(mod, hs, [0], [0]);
  // every order m with h + m q an integer contributes J_m (J_-m = (-1)^m J_m),
  // e.g. h = 1/2 is both m = -2 of (1 0 0) and m = +2 of (0 0 0)
  const J = (m, x) => (m < 0 && (-m) % 2 ? -1 : 1) * besselJ(Math.abs(m), x);
  hs.forEach((h, i) => {
    const Q = 2 * Math.PI * h / 3;
    let amp = 0;
    for (let m = -40; m <= 40; m++) if (Math.abs(h + m * 0.25 - Math.round(h + m * 0.25)) < 1e-9) amp += J(m, Q * A);
    const expected = n * n * b * b * amp * amp;
    assert.ok(Math.abs(out.I[i] - expected) / expected < 1e-9, `h=${h}: ${out.I[i]} vs ${expected}`);
  });
  // average subtraction keeps the satellites and removes the main peak
  const sub = await intensity(mod, [1, 1.25], [0], [0], "nuclear", true);
  assert.ok(sub.I[0] < 1e-9 * out.I[2]);
  assert.ok(Math.abs(sub.I[1] - out.I[3]) / out.I[3] < 1e-9);
});

test("proper helix: magnetic satellites N^2 p^2 f^2 A^2 / 2 at G +- q, none at G", async () => {
  const N = 8, A = 2.2;
  const mod = M.applyModulation(cubic(N), { kind: "helix", q: [0, 0, 0.125], amplitude: A, direction: [0, 0, 1] });
  const sp = mod.magnetic.species[0], n = mod.atoms, p = core.MAGNETIC_LENGTH;
  const moments = Array.from({ length: n }, (_, j) => Math.hypot(...mod.magnetic.vectors.subarray(j * 3, j * 3 + 3)));
  assert.ok(moments.every((m) => Math.abs(m - A) < 1e-12));
  const ls = [0.875, 1, 1.125];
  const out = await intensity(mod, [0], [0], ls, "magnetic");
  const f = (l) => ctx.MagneticFormFactors.evaluate(sp.formFactor, 2 * Math.PI * l / 3);
  assert.ok(Math.abs(out.I[0] - n * n * p * p * f(0.875) ** 2 * A * A / 2) / out.I[0] < 1e-9);
  assert.ok(out.I[1] < 1e-18 * out.I[0]);
  assert.ok(Math.abs(out.I[2] - n * n * p * p * f(1.125) ** 2 * A * A / 2) / out.I[2] < 1e-9);
});

test("incommensurate q: satellites at G + q with the box's finite-size width", async () => {
  const N = 8, q = 0.2873;
  const mod = M.applyModulation(cubic(N), { kind: "sdw", q: [q, 0, 0], amplitude: 1, direction: [0, 1, 0] });
  assert.deepEqual(Array.from(mod.modulation.commensurate), [false, true, true]);
  // exact finite-box result along (h 0 0), moments along y _|_ Q:
  //   M = p f A N^2 (D(h + q) + D(h - q)) / 2,  D(x) = sum_n exp(2 pi i x n), n < N
  const hs = Array.from({ length: 81 }, (_, i) => 0.6 + i * 0.01);
  const out = await intensity(mod, hs, [0], [0], "magnetic");
  const sp = mod.magnetic.species[0], p = core.MAGNETIC_LENGTH;
  const D = (x) => { let re = 0, im = 0; for (let n = 0; n < N; n++) { re += Math.cos(2 * Math.PI * x * n); im += Math.sin(2 * Math.PI * x * n); } return [re, im]; };
  let worst = 0, peak = 0;
  hs.forEach((h, i) => {
    const a = D(h + q), b = D(h - q), f = ctx.MagneticFormFactors.evaluate(sp.formFactor, 2 * Math.PI * h / 3);
    const re = (a[0] + b[0]) / 2, im = (a[1] + b[1]) / 2;
    const expected = (p * f * N * N) ** 2 * (re * re + im * im);
    peak = Math.max(peak, expected);
    worst = Math.max(worst, Math.abs(out.I[i] - expected));
  });
  assert.ok(worst / peak < 1e-9, `finite-box formula: ${worst / peak}`);
  // local maxima within one scan step of the satellites at 1 - q and 1 + q
  for (const target of [1 - q, 1 + q]) {
    let best = -1;
    hs.forEach((h, i) => { if (Math.abs(h - target) < 0.1 && (best < 0 || out.I[i] > out.I[best])) best = i; });
    assert.ok(Math.abs(hs[best] - target) < 0.011, `peak near ${target} at ${hs[best]}`);
  }
});
