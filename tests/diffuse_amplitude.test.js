"use strict";

const test = require("node:test");
const assert = require("node:assert/strict");
const { loadScripts } = require("./helpers/browser_env");

const ctx = loadScripts(["js/diffuse_core.js", "js/diffuse_amplitude.js", "js/rmc6f_reader.js"]);
const core = ctx.DiffuseCore;

function random(seed) {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// RMC6f text for a supercell of a two-site cell (Fe at 0,0,0 and O at the
// body centre) with random displacements; `skip(i, cell)` removes atoms.
function rmc6f(supercell, skip = () => false, seed = 7) {
  const rnd = random(seed);
  const sites = [["Fe", [0, 0, 0]], ["O", [0.5, 0.5, 0.5]]];
  const rows = [];
  for (let cx = 0; cx < supercell[0]; cx++)
    for (let cy = 0; cy < supercell[1]; cy++)
      for (let cz = 0; cz < supercell[2]; cz++)
        sites.forEach(([element, s], site) => {
          if (skip(rows.length, [cx, cy, cz], site)) return;
          const cell = [cx, cy, cz];
          const f = s.map((v, d) => (cell[d] + v + 0.08 * (rnd() - 0.5)) / supercell[d]);
          rows.push(`${rows.length + 1} ${element} [${site + 1}] ${f.join(" ")} ${site + 1} ${cell.join(" ")}`);
        });
  return [
    "(Version 6f format configuration file)",
    `Number of atoms: ${rows.length}`,
    `Supercell dimensions: ${supercell.join(" ")}`,
    `Cell (Ang/deg): ${4 * supercell[0]} ${3.5 * supercell[1]} ${3 * supercell[2]} 90 95 100`,
    "Atoms:",
    ...rows,
  ].join("\n");
}

// Exact sum used in place of a NUFFT: f_k = sum_j c_j exp(i*isign*(t_k . s_j)).
async function directSum(spec) {
  const s = spec.sourcesPacked, t = spec.targetsPacked, c = spec.strengths;
  const m = s.length / 3, n = t.length / 3, out = new Float64Array(n * 2);
  for (let k = 0; k < n; k++) {
    let re = 0, im = 0;
    for (let j = 0; j < m; j++) {
      const p = spec.isign * (t[k * 3] * s[j * 3] + t[k * 3 + 1] * s[j * 3 + 1] + t[k * 3 + 2] * s[j * 3 + 2]);
      const cr = c[j * 2], ci = c[j * 2 + 1], cos = Math.cos(p), sin = Math.sin(p);
      re += cr * cos - ci * sin;
      im += cr * sin + ci * cos;
    }
    out[k * 2] = re;
    out[k * 2 + 1] = im;
  }
  return { out };
}

function grid(lo, hi, step) {
  const out = [];
  for (let i = 0; lo + i * step <= hi + 1e-9; i++) out.push(Number((lo + i * step).toFixed(12)));
  return out;
}

function parse(text) {
  return core.attachNeutronCoefficients(ctx.RMC6fReader.parse("test.rmc6f", text), 10);
}

async function compute(parsed, axes, extra = {}) {
  return ctx.DiffuseAmplitude.computeIntensity({
    parsed, h: axes[0], k: axes[1], l: axes[2], Bq: parsed.Bq, backend: "test",
    sub: true, chunkSize: 37, runType3: directSum, deltaOnLattice: true, ...extra,
  });
}

// |sum b e^{iq.x}|^2 minus nothing, and the previous formula
// |A - Aavg * Adelta / N| with Aavg transformed from the cell origins.
async function reference(parsed, axes) {
  const n = axes[0].length * axes[1].length * axes[2].length;
  const trg = core.targetsChunk(axes[0], axes[1], axes[2], parsed.Bq, 0, n);
  const b = core.complexReal(parsed.fca);
  const A = (await directSum({ isign: 1, sourcesPacked: core.interleave3(parsed.x, parsed.y, parsed.z), targetsPacked: trg, strengths: b })).out;
  const Aa = (await directSum({ isign: 1, sourcesPacked: core.interleave3(parsed.xa, parsed.ya, parsed.za), targetsPacked: trg, strengths: core.complexOnes(parsed.atoms) })).out;
  const Ad = (await directSum({ isign: 1, sourcesPacked: core.interleave3(parsed.dx, parsed.dy, parsed.dz), targetsPacked: trg, strengths: b })).out;
  const total = new Float64Array(n), previous = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    total[i] = A[i * 2] ** 2 + A[i * 2 + 1] ** 2;
    const br = (Aa[i * 2] * Ad[i * 2] - Aa[i * 2 + 1] * Ad[i * 2 + 1]) / parsed.atoms;
    const bi = (Aa[i * 2] * Ad[i * 2 + 1] + Aa[i * 2 + 1] * Ad[i * 2]) / parsed.atoms;
    previous[i] = (A[i * 2] - br) ** 2 + (A[i * 2 + 1] - bi) ** 2;
  }
  return { total, previous };
}

function maxRelDiff(a, b) {
  let peak = 0, diff = 0;
  for (let i = 0; i < a.length; i++) {
    peak = Math.max(peak, Math.abs(b[i]));
    diff = Math.max(diff, Math.abs(a[i] - b[i]));
  }
  return diff / peak;
}

test("lattice-sum average matches the cell-origin transform for complete cells", async () => {
  const parsed = parse(rmc6f([3, 2, 2]));
  const axes = [grid(-1.3, 1.1, 0.23), grid(-0.5, 1.5, 0.37), grid(0.2, 1.4, 0.3)];
  const res = await compute(parsed, axes);
  const ref = await reference(parsed, axes);
  assert.ok(maxRelDiff(res.I, ref.previous) < 1e-10);
  assert.equal(res.profile.averageCellCount, 12);
});

test("vacancies: Bragg points vanish and commensurate points keep |A|^2", async () => {
  const supercell = [4, 3, 2];
  const parsed = parse(rmc6f(supercell, (i, cell, site) => site === 1 && (cell[0] + 2 * cell[1] + cell[2]) % 3 === 0));
  const axes = supercell.map((n) => grid(-1, 1, 1 / n));
  const res = await compute(parsed, axes);
  const ref = await reference(parsed, axes);
  const [nh, nk, nl] = axes.map((a) => a.length);
  let bragg = 0, peak = 0, between = 0, previousError = 0;
  for (let ih = 0; ih < nh; ih++)
    for (let ik = 0; ik < nk; ik++)
      for (let il = 0; il < nl; il++) {
        const i = (ih * nk + ik) * nl + il;
        const integer = [axes[0][ih], axes[1][ik], axes[2][il]].every((v) => Math.abs(v - Math.round(v)) < 1e-9);
        peak = Math.max(peak, ref.total[i]);
        if (integer) bragg = Math.max(bragg, res.I[i]);
        else {
          between = Math.max(between, Math.abs(res.I[i] - ref.total[i]));
          previousError = Math.max(previousError, Math.abs(ref.previous[i] - ref.total[i]));
        }
      }
  assert.ok(bragg < 1e-9 * peak, `Bragg residual ${bragg}`);
  assert.ok(between < 1e-9 * peak, `commensurate mismatch ${between}`);
  assert.ok(previousError > 1e-3 * peak, "the old cell-origin formula should differ with vacancies");
});

test("A_delta is evaluated only on integer hkl when the steps are commensurate", async () => {
  const supercell = [4, 3, 2];
  const parsed = parse(rmc6f(supercell, () => false, 5));
  const cases = [
    { axes: supercell.map((n) => grid(-1, 1, 1 / n)), points: 3 * 3 * 3 },
    // h commensurate, k and l not: only integer h rows are needed
    { axes: [grid(-1, 1, 1 / 4), grid(-0.5, 0.7, 0.3), grid(0.1, 0.8, 0.35)], points: 3 * 5 * 3 },
  ];
  for (const { axes, points } of cases) {
    const res = await compute(parsed, axes);
    const ref = await reference(parsed, axes);
    assert.equal(res.profile.deltaPoints, points);
    assert.ok(maxRelDiff(res.I, ref.previous) < 1e-10);
    const full = await compute(parsed, axes, { deltaOnLattice: false });
    assert.equal(full.profile.deltaPoints, null);
    assert.ok(maxRelDiff(full.I, res.I) < 1e-10);
  }
});

test("average subtraction requires per-atom cell indices", async () => {
  const parsed = parse(rmc6f([2, 2, 2]));
  const noCells = { ...parsed, cellIndex: undefined };
  await assert.rejects(compute(noCells, [[0], [0], [0.5]]), /unit-cell indices/);
  const lammpsLike = { ...parsed, hasSupercell: false };
  await assert.rejects(compute(lammpsLike, [[0], [0], [0.5]]), /unit-cell indices/);
});

test("fast and grouped neutron engines agree", async () => {
  const parsed = parse(rmc6f([2, 3, 2], () => false, 11));
  const axes = [grid(-0.7, 0.9, 0.4), grid(0, 1, 0.25), grid(-0.5, 0.5, 0.5)];
  const fast = await compute(parsed, axes, { scattering: { type: "neutron", model: "fast" } });
  const grouped = await compute(parsed, axes, { scattering: { type: "neutron", model: "grouped_exact" } });
  assert.ok(maxRelDiff(grouped.I, fast.I) < 1e-10);
});

test("laueAxis is exact at integers and zero at other supercell points", () => {
  const L = core.laueAxis([0, 1, -2, 0.25, 0.5, 0.75, 1.25], 4, 0);
  for (const i of [0, 1, 2]) {
    assert.ok(Math.abs(L[i * 2] - 4) < 1e-12 && Math.abs(L[i * 2 + 1]) < 1e-12);
  }
  for (const i of [3, 4, 5, 6]) {
    assert.ok(Math.hypot(L[i * 2], L[i * 2 + 1]) < 1e-12);
  }
});
