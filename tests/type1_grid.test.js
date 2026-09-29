"use strict";

const test = require("node:test");
const assert = require("node:assert/strict");
const { loadScripts } = require("./helpers/browser_env");

const ctx = loadScripts(["js/rmc6f_reader.js", "js/type1_grid.js"]);
const T1 = ctx.DiffuseType1;

// A small triclinic supercell with displaced atoms, parsed by the real reader.
function structure() {
  const supercell = [3, 2, 2];
  const rows = [];
  let seed = 3;
  const rnd = () => ((seed = (seed * 16807) % 2147483647) / 2147483647) - 0.5;
  for (let cx = 0; cx < 3; cx++)
    for (let cy = 0; cy < 2; cy++)
      for (let cz = 0; cz < 2; cz++)
        for (const [el, s] of [["Fe", [0, 0, 0]], ["O", [0.5, 0.4, 0.6]]]) {
          const f = [cx, cy, cz].map((c, d) => (c + s[d] + 0.1 * rnd()) / supercell[d]);
          rows.push(`${rows.length + 1} ${el} [1] ${f.join(" ")} 1 ${cx} ${cy} ${cz}`);
        }
  return ctx.RMC6fReader.parse("t.rmc6f", [
    "(Version 6f format configuration file)", `Number of atoms: ${rows.length}`,
    "Supercell dimensions: 3 2 2", "Cell (Ang/deg): 13.2 7.8 6.2 84 95 103", "Atoms:", ...rows,
  ].join("\n"));
}

// Type-1 sum by definition over centred modes, dimension 0 (l) fastest.
function type1ByDefinition(points, strengths, modes) {
  const [n0, n1, n2] = modes;
  const out = new Float64Array(n0 * n1 * n2 * 2);
  for (let i2 = 0; i2 < n2; i2++)
    for (let i1 = 0; i1 < n1; i1++)
      for (let i0 = 0; i0 < n0; i0++) {
        const m = [i0 - Math.floor(n0 / 2), i1 - Math.floor(n1 / 2), i2 - Math.floor(n2 / 2)];
        let re = 0, im = 0;
        for (let j = 0; j < points.length / 3; j++) {
          const p = m[0] * points[j * 3] + m[1] * points[j * 3 + 1] + m[2] * points[j * 3 + 2];
          re += strengths[j * 2] * Math.cos(p) - strengths[j * 2 + 1] * Math.sin(p);
          im += strengths[j * 2] * Math.sin(p) + strengths[j * 2 + 1] * Math.cos(p);
        }
        const idx = i0 + n0 * (i1 + n1 * i2);
        out[idx * 2] = re;
        out[idx * 2 + 1] = im;
      }
  return out;
}

const axis = (start, step, n) => Array.from({ length: n }, (_, i) => Number((start + i * step).toFixed(12)));

test("uniformAxis recognises uniform axes and rejects others", () => {
  assert.deepEqual(JSON.parse(JSON.stringify(T1.uniformAxis(axis(-1.3, 0.1, 5)))), { start: -1.3, step: 0.1, count: 5 });
  assert.equal(T1.uniformAxis([0, 0.1, 0.25]), null);
  assert.equal(T1.uniformAxis([2]).count, 1);
});

test("type-1 slabs reproduce the direct sum on offset, odd and even grids", () => {
  const parsed = structure();
  const n = parsed.atoms;
  const cart = new Float64Array(n * 3);
  for (let j = 0; j < n; j++) { cart[j * 3] = parsed.x[j]; cart[j * 3 + 1] = parsed.y[j]; cart[j * 3 + 2] = parsed.z[j]; }
  const strengths = new Float64Array(n * 2);
  for (let j = 0; j < n; j++) { strengths[j * 2] = parsed.elements[j] === "Fe" ? 0.945 : 0.5803; strengths[j * 2 + 1] = 0.1 * (j % 3); }
  const u = T1.fractionalCoordinates(cart, parsed.Bq);
  for (const [h, k, l] of [
    [axis(-1.3, 0.23, 7), axis(0.4, 0.5, 4), axis(-0.25, 1 / 2, 5)],
    [axis(2.1, 0.07, 6), axis(-1, 1 / 2, 5), axis(0.3, 0.9, 3)],
  ]) {
    const axes = T1.gridFromAxes(h, k, l);
    const points = T1.type1Points(u, axes);
    const nk = k.length, nl = l.length;
    // plain slabs, then padded mode counts (extra modes extend each axis)
    for (const [planes, padK, padL] of [[h.length, 0, 0], [3, 0, 0], [2, 0, 0], [4, 3, 2]]) {
      const mk = nk + padK, ml = nl + padL;
      for (let ih0 = 0; ih0 < h.length; ih0 += planes) {
        const modes = [ml, mk, planes];
        const hc = T1.slabCentre(axes, ih0, modes);
        const phased = T1.phasedStrengths(u, strengths, hc);
        const out = type1ByDefinition(points, phased, modes);
        for (let p = 0; p < planes && ih0 + p < h.length; p++)
          for (let ik = 0; ik < nk; ik++)
            for (let il = 0; il < nl; il++) {
              const idx = il + ml * (ik + mk * p);
              const want = T1.directAmplitude(u, strengths, [h[ih0 + p], k[ik], l[il]]);
              const err = Math.hypot(out[idx * 2] - want[0], out[idx * 2 + 1] - want[1]);
              assert.ok(err < 2e-5, `slab ${ih0}+${p} (${ik},${il}) error ${err}`);
            }
      }
    }
  }
});

test("fineGridLength follows the 2,3,5,7,11,13-smooth even rule", () => {
  assert.equal(T1.fineGridLength(241, 2, 1e-6), 484);
  assert.equal(T1.fineGridLength(10, 2, 1e-6), 20);
  const len = T1.fineGridLength(201, 1.25, 1e-6);
  assert.ok(len >= 252 && len % 2 === 0);
});
