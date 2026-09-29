"use strict";

// Laue symmetrization (js/laue_symmetry.js): Scatty's 12 classes and their
// operations on h, k, l, the mean over the equivalents on the grid.
const test = require("node:test");
const assert = require("node:assert/strict");
const { loadScripts } = require("./helpers/browser_env");

const ctx = loadScripts(["js/laue_symmetry.js"]);
const L = ctx.LaueSymmetry;

const apply = (M, v) => M.map((row) => row[0] * v[0] + row[1] * v[1] + row[2] * v[2]);
const key = (v) => v.join(",");
const cubic = [[1, 0, 0], [0, 1, 0], [0, 0, 1]].map((r) => r.map((x) => x * 1.4));
// hexagonal: a = b, gamma = 120 deg; rows are the reciprocal vectors (2 pi / a units)
const hexagonal = [[1, 1 / Math.sqrt(3), 0], [0, 2 / Math.sqrt(3), 0], [0, 0, 0.7]];
const tetragonal = [[1, 0, 0], [0, 1, 0], [0, 0, 0.6]];

test("group orders, closure, identity and inversion", () => {
  const ORDER = { "-1": 2, "2/m": 4, "mmm": 8, "4/m": 8, "4/mmm": 16, "-3": 6, "-31m": 12, "-3m1": 12, "6/m": 12, "6/mmm": 24, "m-3": 24, "m-3m": 48 };
  assert.deepEqual([...L.CLASSES].sort(), Object.keys(ORDER).sort());
  for (const laue of L.CLASSES) {
    const ops = L.operations(laue), keys = new Set(ops.map((M) => key(M.flat())));
    assert.equal(ops.length, ORDER[laue], laue);
    assert.equal(keys.size, ops.length, `${laue}: repeated operations`);
    assert.ok(keys.has("1,0,0,0,1,0,0,0,1") && keys.has("-1,0,0,0,-1,0,0,0,-1"), `${laue}: identity and inversion`);
    for (const A of ops) for (const B of ops) {
      const C = A.map((row) => [0, 1, 2].map((c) => row[0] * B[0][c] + row[1] * B[1][c] + row[2] * B[2][c]));
      assert.ok(keys.has(key(C.flat())), `${laue}: not closed`);
    }
  }
});

test("the operations keep |Q| for cells of their system", () => {
  for (const laue of ["-1", "2/m", "mmm", "4/m", "4/mmm", "m-3", "m-3m"]) assert.ok(L.metricDeviation(L.operations(laue), cubic) < 1e-12, laue);
  for (const laue of ["4/m", "4/mmm"]) assert.ok(L.metricDeviation(L.operations(laue), tetragonal) < 1e-12, laue);
  for (const laue of ["-3", "-31m", "-3m1", "6/m", "6/mmm"]) assert.ok(L.metricDeviation(L.operations(laue), hexagonal) < 1e-12, laue);
  assert.ok(L.metricDeviation(L.operations("m-3m"), tetragonal) > 0.1);
  assert.ok(L.metricDeviation(L.operations("6/mmm"), cubic) > 0.1);
});

test("the trigonal classes are Scatty's: -3m1 has (-k,-h,l), -31m has (-k,-h,-l)", () => {
  const eq = (laue) => new Set(L.operations(laue).map((M) => key(apply(M, [1, 2, 3]))));
  const i = -3;
  assert.ok(eq("-3m1").has(key([-2, -1, 3])) && !eq("-3m1").has(key([-2, -1, -3])));
  assert.ok(eq("-31m").has(key([-2, -1, -3])) && !eq("-31m").has(key([-2, -1, 3])));
  assert.ok(eq("-3").has(key([2, i, 3])) && eq("-3").has(key([i, 1, 3])));
});

function grid(axis) {
  return { h: axis, k: axis, l: axis, shape: [axis.length, axis.length, axis.length] };
}

test("symmetrizing: invariant result, unchanged when already symmetric, partial equivalents off the grid", () => {
  const ax = Array.from({ length: 9 }, (_, i) => (i - 4) * 0.25), g = grid(ax), n = ax.length;
  const at = (i, j, m) => (i * n + j) * n + m;
  const noisy = new Float64Array(n ** 3).map((_, p) => Math.sin(p * 12.9898) * 43758.5453 % 1 + 2);
  const sym = L.symmetrize(noisy, g, "m-3m", { Bq: cubic });
  assert.equal(sym.meta.operations, 48);
  assert.ok(Math.abs(sym.meta.coverage - 1) < 1e-12);
  const idx = (v) => Math.round(v / 0.25) + 4;
  for (const M of L.operations("m-3m")) {
    for (const p of [[1, 2, 3], [0, 4, 1], [-3, 2, 2]]) {
      const hkl = p.map((v) => v * 0.25), q = apply(M, hkl);
      assert.ok(Math.abs(sym.I[at(...p.map((v) => v + 4))] - sym.I[at(idx(q[0]), idx(q[1]), idx(q[2]))]) < 1e-12);
    }
  }
  // a function of |Q| is already cubic-symmetric
  const radial = new Float64Array(n ** 3);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) for (let m = 0; m < n; m++) radial[at(i, j, m)] = Math.exp(-(ax[i] ** 2 + ax[j] ** 2 + ax[m] ** 2));
  const same = L.symmetrize(radial, g, "m-3m", {});
  for (let p = 0; p < radial.length; p++) assert.ok(Math.abs(same.I[p] - radial[p]) < 1e-12);
  // a one-sided grid: -1 sends (h,k,l) off it, so only the point itself counts
  const half = grid([0, 0.5, 1]);
  const own = L.symmetrize(new Float64Array(27).map((_, p) => p), half, "-1", {});
  assert.ok(Math.abs(own.I[26] - 26) < 1e-12);
  assert.ok(Math.abs(own.I[0] - 0) < 1e-12);
  assert.ok(own.meta.coverage < 0.6);
  // masked points are left out of the mean
  const masked = new Float64Array(n ** 3).fill(1);
  masked[at(4 + 1, 4 + 2, 4 + 3)] = NaN;
  const filled = L.symmetrize(masked, g, "m-3m", {});
  assert.equal(filled.I[at(4 + 1, 4 + 2, 4 + 3)], 1);
});
