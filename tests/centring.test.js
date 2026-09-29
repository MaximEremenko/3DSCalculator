"use strict";

const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const { loadScripts } = require("./helpers/browser_env");

const ctx = loadScripts(["js/rmc6f_reader.js", "js/magnetic_form_factors.js", "js/diffuse_core.js", "js/diffuse_amplitude.js", "js/scatty_reader.js"]);
const core = ctx.DiffuseCore;
const FIX = path.join(__dirname, "fixtures", "scatty_nuclear");
const plain = (v) => JSON.parse(JSON.stringify(v));

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

test("centring factor: 1 at allowed reflections, 0 at forbidden ones", () => {
  const at = (letter, h, k, l) => Array.from(core.centringChunk([h], [k], [l], 0, 1, letter)).map((v) => Math.round(v * 1e12) / 1e12 + 0);
  assert.deepEqual(at("F", 1, 1, 1), [1, 0]);
  assert.deepEqual(at("F", 2, 0, 0), [1, 0]);
  assert.deepEqual(at("F", 1, 0, 0), [0, 0]);
  assert.deepEqual(at("F", 2, 1, 0), [0, 0]);
  assert.deepEqual(at("I", 1, 1, 0), [1, 0]);
  assert.deepEqual(at("I", 1, 0, 0), [0, 0]);
  assert.deepEqual(at("C", 1, 1, 7), [1, 0]);
  assert.deepEqual(at("C", 1, 0, 7), [0, 0]);
  assert.deepEqual(at("R", 1, 1, 0), [1, 0]); // obverse: -h + k + l = 3n
  assert.deepEqual(at("R", -1, 1, 0), [0, 0]);
  assert.deepEqual(at("R", 1, 0, 0), [0, 0]);
  assert.deepEqual(at("P", 0.3, 0.1, 5), [1, 0]);
});

test("nuclear intensities match Scatty, with F-centred Bragg removal", async () => {
  const parsed = core.attachNeutronCoefficients(ctx.ScattyReader.parse("mnfeo_atoms_01.txt", fs.readFileSync(path.join(FIX, "mnfeo_atoms_01.txt"), "utf8")), 10);
  const ref = JSON.parse(fs.readFileSync(path.join(FIX, "scatty_mnfeo_reference.json"), "utf8"));
  const b = ref["scatty_b_1e-12cm"];
  parsed.elements.forEach((e, i) => { parsed.fca[i] = b[e]; });
  assert.deepEqual(plain(core.structureCentring(parsed)), { nuclear: "F", magnetic: null });
  const axis = Array.from({ length: 21 }, (_, i) => Number((-2 + i * 0.2).toFixed(12)));
  const run = (sub, centring) => ctx.DiffuseAmplitude.computeIntensity({
    parsed, h: axis, k: axis, l: axis, Bq: parsed.Bq, backend: "test", sub, deltaOnLattice: true,
    scattering: { type: "neutron", model: "fast" }, runType3: exactSum, centring,
  });
  const [total, diffuse, primitive] = await Promise.all([run(false, "F"), run(true, "F"), run(true, "P")]);
  const index = (h, k, l) => (Math.round((h + 2) / 0.2) * 21 + Math.round((k + 2) / 0.2)) * 21 + Math.round((l + 2) / 0.2);
  let worstTotal = 0, worstDiffuse = 0, peak = 0, forbiddenGap = 0;
  for (const [h, k, l, iTotal, iDiffuse] of ref.points) {
    const i = index(h, k, l);
    peak = Math.max(peak, iDiffuse);
    worstTotal = Math.max(worstTotal, Math.abs(total.I[i] / parsed.atoms - iTotal) / iTotal);
    worstDiffuse = Math.max(worstDiffuse, Math.abs(diffuse.I[i] / parsed.atoms - iDiffuse));
    const lattice = [h, k, l].every((v) => Math.abs(v - Math.round(v)) < 1e-9);
    const parity = [h, k, l].map((v) => ((Math.round(v) % 2) + 2) % 2);
    if (lattice && !parity.every((v) => v === parity[0])) forbiddenGap = Math.max(forbiddenGap, Math.abs(primitive.I[i] / parsed.atoms - iDiffuse));
  }
  assert.ok(worstTotal < 1e-10, `total ${worstTotal}`);
  assert.ok(worstDiffuse / peak < 1e-10, `diffuse ${worstDiffuse / peak}`);
  // subtracting at every lattice point of the conventional cell would remove
  // real diffuse intensity at the F-forbidden points
  assert.ok(forbiddenGap / peak > 0.1, `P subtraction gap ${forbiddenGap / peak}`);
});

test("detects the centring of simple average structures", () => {
  const build = (sites, box, elementOf, jitter = 0) => {
    const frac = [], cellIndex = [], elements = [];
    let seed = 7;
    const rnd = () => ((seed = (seed * 16807) % 2147483647) / 2147483647 - 0.5) * jitter;
    for (let a = 0; a < box; a++) for (let b = 0; b < box; b++) for (let c = 0; c < box; c++) {
      sites.forEach((s, i) => {
        frac.push((a + s[0] + rnd()) / box, (b + s[1] + rnd()) / box, (c + s[2] + rnd()) / box);
        cellIndex.push(a, b, c);
        elements.push(elementOf(i));
      });
    }
    const parsed = ctx.RMC6fReader.buildSupercellStructure({ parentCellDeg: [4, 4, 4, 90, 90, 90], supercell: [box, box, box], fractional: Float64Array.from(frac), cellIndex: Int32Array.from(cellIndex), elements });
    return { ...parsed, cellDeg: [4 * box, 4 * box, 4 * box, 90, 90, 90] };
  };
  const fcc = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]];
  const nuclear = (p) => core.structureCentring(p).nuclear;
  assert.equal(nuclear(build(fcc, 4, () => "Cu", 0.02)), "F");
  assert.equal(nuclear(build([[0, 0, 0], [0.5, 0.5, 0.5]], 4, () => "Fe", 0.02)), "I");
  assert.equal(nuclear(build([[0, 0, 0], [0.5, 0.5, 0.5]], 4, (i) => (i ? "Cl" : "Cs"))), "P"); // CsCl
  assert.equal(nuclear(build([[0, 0, 0], [0.5, 0.5, 0]], 4, () => "Cu")), "C");
  assert.equal(nuclear(build(fcc, 4, (i) => (i ? "Cu" : "Au"))), "P"); // ordered Cu3Au (L1_2)
  // magnetic: moments along z, + on the (0 0 0) and (1/2 1/2 0) sites, - on
  // the others (type-I order): it keeps C centring but breaks the nuclear F
  const af = build(fcc, 4, () => "Mn");
  af.magnetic = {
    species: [{ scale: 5 }],
    speciesOfAtom: new Int32Array(af.atoms),
    vectors: Float64Array.from({ length: af.atoms * 3 }, (_, i) => (i % 3 === 2 ? (Math.floor(i / 3) % 4 < 2 ? 1 : -1) : 0)),
  };
  assert.deepEqual(plain(core.structureCentring(af)), { nuclear: "F", magnetic: "C" });
});
