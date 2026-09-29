"use strict";

const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const { loadScripts } = require("./helpers/browser_env");

const ctx = loadScripts([
  "js/rmc6f_reader.js", "js/magnetic_form_factors.js", "js/diffuse_core.js",
  "js/diffuse_amplitude.js", "js/scatty_reader.js",
]);
const core = ctx.DiffuseCore;
const FF = ctx.MagneticFormFactors;
const FIXTURES = path.join(__dirname, "fixtures", "magnetic");

// Exact sum in place of a NUFFT: f_k = sum_j c_j exp(i t_k . s_j).
async function directSum(spec) {
  const s = spec.sourcesPacked, t = spec.targetsPacked, c = spec.strengths;
  const m = s.length / 3, n = t.length / 3, out = new Float64Array(n * 2);
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

function compute(parsed, axes, magnetic, extra = {}) {
  return ctx.DiffuseAmplitude.computeIntensity({
    parsed, h: axes[0], k: axes[1], l: axes[2], Bq: parsed.Bq, backend: "test", sub: false,
    scattering: { type: "neutron", model: "fast" }, runType3: directSum, magnetic, ...extra,
  });
}

function scatty(name) {
  const text = fs.readFileSync(path.join(FIXTURES, name), "utf8");
  return core.attachNeutronCoefficients(ctx.ScattyReader.parse(name, text), 10);
}

const P2 = core.MAGNETIC_LENGTH ** 2;
const AXIS = Array.from({ length: 13 }, (_, i) => Number((-1 + i / 6).toFixed(12)));
const REFERENCE = JSON.parse(fs.readFileSync(path.join(FIXTURES, "scatty_mno_reference.json"), "utf8"));

test("magnetic length squared equals Scatty's (gamma r0 / 2)^2", () => {
  assert.ok(Math.abs(P2 - 0.07265289) < 1e-8);
});

test("form factors reproduce reference values", () => {
  const mn = { j0: FF.find("Mn2+").j0, j2: null, c2: 0 };
  // Mn2+ <j0> (Brown), as used for MnO
  [[0.5, 0.98193], [1, 0.93231], [2, 0.76297], [3, 0.55899], [4, 0.3767], [5, 0.23838]].forEach(([q, f]) => {
    assert.ok(Math.abs(FF.evaluate(mn, q) - f) < 2e-5, `Mn2+ at ${q}`);
  });
  // Gd3+ <j0> values from RMCProfile's magnetic_form_factors_test.cpp
  const gd = { j0: FF.coefficients([0.0186, 25.3867, 0.2895, 11.1421, 0.7135, 3.752, -0.0217]), j2: null, c2: 0 };
  assert.ok(Math.abs(FF.evaluate(gd, 12) - 0.001618544492) < 1e-11);
  assert.ok(Math.abs(FF.evaluate(gd, 14) + 0.01492458694) < 1e-11);
  assert.deepEqual(Array.from(FF.find("Gd3+").j0.slice(0, 6)), [0.0186, 25.3867, 0.2895, 11.1421, 0.7135, 3.752]);
});

test("ion labels, Lande g and C2 suggestions", () => {
  assert.equal(FF.ionKey("Mn2+"), "Mn2");
  assert.equal(FF.ionKey("Mn+2"), "Mn2");
  assert.equal(FF.ionKey("ho3"), "Ho3");
  assert.equal(FF.ionKey("Fe"), "Fe0");
  assert.ok(Math.abs(FF.landeG("5I8") - 1.25) < 1e-12);
  assert.ok(Math.abs(FF.landeG("6H15/2") - 4 / 3) < 1e-12);
  assert.ok(Math.abs(FF.suggestedC2(FF.find("Ho3+")) - 0.6) < 1e-12);
  assert.equal(FF.suggestedC2(FF.find("Mn2+")), 0);
  // Scatty's spin-ice example gives Ho3+ <j0> as 7 coefficients
  assert.equal(FF.identify([0.0566, 18.3176, 0.3365, 7.688, 0.6317, 2.9427, -0.0248]).ion, "Ho3");
});

test("Scatty spin file: geometry, species and moments", () => {
  const parsed = scatty("mno_order_spins_01.txt");
  assert.equal(parsed.atoms, 864);
  assert.deepEqual(Array.from(parsed.super), [6, 6, 6]);
  assert.equal(parsed.magnetic.species.length, 1);
  assert.equal(parsed.magnetic.species[0].label, "Mn2+");
  assert.equal(parsed.elements[0], "Mn");
  assert.equal(parsed.magneticOnly, true);
  for (let i = 0; i < parsed.atoms; i++) {
    const m = parsed.magnetic.moments;
    assert.ok(Math.abs(Math.hypot(m[i * 3], m[i * 3 + 1], m[i * 3 + 2]) - 5) < 1e-9);
  }
});

for (const [stem, key] of [["mno_random", "mno_random"], ["mno_order", "mno_order"]]) {
  for (const paramagnet of [false, true]) {
    test(`${stem}${paramagnet ? " minus ideal paramagnet" : ""} matches scatty.exe`, async () => {
      const parsed = scatty(`${stem}_spins_01.txt`);
      const res = await compute(parsed, [AXIS, AXIS, AXIS], { mode: "magnetic", qZero: "average", subtractParamagnet: paramagnet });
      const ref = REFERENCE[paramagnet ? `${key}_temp_subtract` : key];
      let peak = 0, diff = 0;
      ref.forEach((v, i) => {
        peak = Math.max(peak, Math.abs(v));
        diff = Math.max(diff, Math.abs(res.I[i] / parsed.atoms - v));
      });
      assert.ok(diff < 1e-6 * peak, `max deviation ${diff} of peak ${peak}`);
    });
  }
}

test("random spins follow the 2/3 law and cancel it with the paramagnet term", async () => {
  const parsed = scatty("mno_random_spins_01.txt");
  const ff = parsed.magnetic.species[0].formFactor;
  const [total, net] = await Promise.all([false, true].map((p) =>
    compute(parsed, [AXIS, AXIS, AXIS], { mode: "magnetic", subtractParamagnet: p })));
  let ratio = 0, residual = 0, n = 0;
  AXIS.forEach((h, ih) => AXIS.forEach((k, ik) => AXIS.forEach((l, il) => {
    const q = (2 * Math.PI * Math.hypot(h, k, l)) / 4.4344;
    if (q === 0) return;
    const ideal = (2 / 3) * P2 * 25 * FF.evaluate(ff, q) ** 2;
    const i = (ih * 13 + ik) * 13 + il;
    ratio += total.I[i] / parsed.atoms / ideal;
    residual += net.I[i] / parsed.atoms / ideal;
    n++;
  })));
  assert.ok(Math.abs(ratio / n - 1) < 0.05, `mean ratio ${ratio / n}`);
  assert.ok(Math.abs(residual / n) < 0.05, `mean residual ${residual / n}`);
});

test("type-II order: (1/2 1/2 1/2) peak and nuclear-cell averaging keeps it", async () => {
  const parsed = scatty("mno_order_spins_01.txt");
  const half = [[0.5], [0.5], [0.5]];
  const q = (2 * Math.PI * Math.sqrt(0.75)) / 4.4344;
  const f = FF.evaluate(parsed.magnetic.species[0].formFactor, q);
  const expect = parsed.atoms * P2 * 25 * f * f; // S is perpendicular to [111]
  const plain = await compute(parsed, half, { mode: "magnetic" });
  const averaged = await compute(parsed, half, { mode: "magnetic" }, { sub: true });
  assert.ok(Math.abs(plain.I[0] / parsed.atoms - expect) < 1e-9 * expect);
  assert.ok(Math.abs(averaged.I[0] - plain.I[0]) < 1e-9 * plain.I[0]);
});

test("RMCProfile's analytic magnetic case: one spin in a 2x1x1 cell", async () => {
  // Cubic 1 A cell, supercell 2x1x1 with a single atom (b = 2, mu = 1, f = 1)
  // at the origin of cell 0; Q = (pi/2, 0, 0) is h = 1/4. RMCProfile's
  // diffuse_scattering_magnetic_test.cpp expects I_N = 2 and I_M = p^2 / 2.
  const s = ctx.RMC6fReader.buildSupercellStructure({
    parentCellDeg: [1, 1, 1, 90, 90, 90], supercell: [2, 1, 1],
    fractional: Float64Array.from([0, 0, 0]), cellIndex: Int32Array.from([0, 0, 0]), elements: ["X"],
  });
  const one = { j0: [0, 0, 0, 0, 0, 0, 0, 0, 1], j2: null, c2: 0 };
  const parsed = (spin) => ({
    file: "analytic", ...s, fca: Float64Array.from([2]),
    magnetic: { species: [{ label: "m", formFactor: one }], speciesOfAtom: Int32Array.from([0]), moments: Float64Array.from(spin) },
  });
  const axes = [[0.25], [0], [0]];
  const both = await compute(parsed([0, 1, 0]), axes, { mode: "both" }, { sub: true });
  const nuclear = await compute(parsed([0, 1, 0]), axes, { mode: "nuclear" }, { sub: true });
  const magnetic = await compute(parsed([0, 1, 0]), axes, { mode: "magnetic" }, { sub: true });
  assert.ok(Math.abs(nuclear.I[0] - 2) < 1e-12);
  assert.ok(Math.abs(magnetic.I[0] - P2 / 2) < 1e-12);
  assert.ok(Math.abs(both.I[0] - (2 + P2 / 2)) < 1e-12);
  // a spin along Q does not scatter
  const along = await compute(parsed([1, 0, 0]), axes, { mode: "magnetic" }, { sub: true });
  assert.ok(Math.abs(along.I[0]) < 1e-15);
  // Q = 0: 0 by default, the orientational average (2/3)|M|^2 on request
  const origin = [[0], [0], [0]];
  const zero = await compute(parsed([0, 1, 0]), origin, { mode: "magnetic" });
  const avg = await compute(parsed([0, 1, 0]), origin, { mode: "magnetic", qZero: "average" });
  assert.equal(zero.I[0], 0);
  assert.ok(Math.abs(avg.I[0] - (2 / 3) * P2) < 1e-15);
});

test("magnetic modes are refused without moments or for X-rays", async () => {
  const parsed = scatty("mno_order_spins_01.txt");
  await assert.rejects(compute({ ...parsed, magnetic: undefined }, [[0.5], [0.5], [0.5]], { mode: "magnetic" }), /magnetic moments/);
  await assert.rejects(compute(parsed, [[0.5], [0.5], [0.5]], { mode: "magnetic" }, { scattering: { type: "xray" } }), /neutrons/);
});
