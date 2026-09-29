"use strict";

// The 3D smoothing kernels of index.html: one bin per box reciprocal-lattice
// spacing (1/supercell r.l.u.), resampled to the grid step.
const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");

const NAMES = ["acosh1", "dftComplex", "fftMagShiftNorm", "lanczosWindow", "chebwinWindow", "resampleKernel", "smoothKernels",
  "kernelEdgeWeights", "normalizeSmoothEdges", "convAxis0Same", "convAxis1Same", "convAxis2Same", "minMaxArray", "smoothIntensity3d"];
const html = fs.readFileSync(path.join(__dirname, "..", "index.html"), "utf8").split("\n");
const ctx = vm.createContext({ Math, Float64Array, Number, Array });
for (const name of NAMES) {
  const line = html.find((l) => l.startsWith(`function ${name}(`));
  assert.ok(line, `index.html defines ${name} on one line`);
  vm.runInContext(line, ctx);
}
const { smoothKernels, fftMagShiftNorm, chebwinWindow, lanczosWindow, smoothIntensity3d } = ctx;

// width of a kernel in r.l.u.: sqrt of the second moment
function width(k, step) {
  const c = Math.floor(k.length / 2);
  let m2 = 0;
  for (let i = 0; i < k.length; i++) m2 += k[i] * ((i - c) * step) ** 2;
  return Math.sqrt(m2);
}
const sum = (k) => k.reduce((a, b) => a + b, 0);

test("on the default grid (step 1/supercell) the kernels are the window transforms", () => {
  for (const type of ["chebyshev", "lanczos"]) {
    const ks = smoothKernels([16, 12, 6], [1 / 16, 1 / 12, 1 / 6], { type, chebDb: 100, scale: 1 });
    [16, 12, 6].forEach((n, d) => {
      const w = type === "chebyshev" ? chebwinWindow(n, 100) : lanczosWindow(n);
      assert.deepEqual(Array.from(ks[d]), Array.from(fftMagShiftNorm(w, n)));
    });
  }
});

test("on other steps the kernel keeps its width in r.l.u. and unit sum", () => {
  const ref = smoothKernels([16, 16, 16], [1 / 16, 1 / 16, 1 / 16], { type: "chebyshev", chebDb: 100, scale: 1 })[0];
  const w0 = width(ref, 1 / 16);
  for (const step of [0.04, 0.02, 0.1]) {
    const k = smoothKernels([16, 16, 16], [step, step, step], { type: "chebyshev", chebDb: 100, scale: 1 })[0];
    assert.ok(Math.abs(sum(k) - 1) < 1e-12);
    assert.equal(k.length % 2, 1);
    assert.ok(Math.abs(width(k, step) / w0 - 1) < 0.08, `step ${step}: width ${width(k, step)} vs ${w0}`);
  }
  // the old behaviour, one bin per grid point, was 36% narrower on a 0.04 grid
  assert.ok(width(ref, 0.04) / w0 < 0.65);
});

test("LAMMPS boxes (supercell 1) at Scale 1: a one-point kernel", () => {
  const ks = smoothKernels([1, 1, 1], [1, 1, 1], { type: "chebyshev", chebDb: 100, scale: 1 });
  for (const k of ks) assert.deepEqual(Array.from(k), [1]);
  const ks2 = smoothKernels([1, 1, 1], [0.5, 0.5, 0.5], { type: "lanczos", scale: 1 });
  for (const k of ks2) assert.deepEqual(Array.from(k), [1]);
});

test("smoothing keeps a flat field flat, edges included", () => {
  const shape = [9, 7, 5], n = 9 * 7 * 5;
  const flat = new Float64Array(n).fill(3.5);
  const out = smoothIntensity3d(flat, shape, smoothKernels([6, 6, 6], [0.1, 0.1, 0.1], { type: "chebyshev", chebDb: 100, scale: 1 }));
  for (const v of out.data) assert.ok(Math.abs(v - 3.5) < 1e-12);
});

test("kernels peak at their centre for odd and even lengths, and smoothing does not move a peak", () => {
  for (const n of [5, 6, 20, 21]) {
    for (const w of [chebwinWindow(n, 100), lanczosWindow(n)]) {
      const k = fftMagShiftNorm(w, n), c = Math.floor(n / 2);
      assert.equal(k.indexOf(Math.max(...k)), c, `n = ${n}: peak at ${k.indexOf(Math.max(...k))}, centre ${c}`);
      for (let d = 1; c - d >= (n % 2 ? 0 : 1); d++) assert.ok(Math.abs(k[c - d] - k[c + d]) < 1e-12, `n = ${n}: asymmetric at +-${d}`);
    }
  }
  // a single peak keeps its position (odd supercell, as the 21-cell Cr example)
  const shape = [41, 1, 1], I = new Float64Array(41);
  I[20] = 1;
  const out = smoothIntensity3d(I, shape, smoothKernels([21, 1, 1], [1 / 21, 1, 1], { type: "chebyshev", chebDb: 100, scale: 1 })).data;
  let m0 = 0, m1 = 0;
  out.forEach((v, i) => { m0 += v; m1 += v * i; });
  assert.ok(Math.abs(m1 / m0 - 20) < 1e-9, `centroid ${m1 / m0}`);
});
