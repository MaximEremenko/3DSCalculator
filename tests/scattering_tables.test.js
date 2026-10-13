"use strict";

// X-ray and electron scattering factors: every table against independent
// references. X-ray: Waasmaier-Kirfel, with f(0) = Z. Electron: the
// Mott-Bethe transform of Waasmaier-Kirfel, f_e = (Z - f_x) / (8 pi^2 a0 s^2).
const test = require("node:test");
const assert = require("node:assert/strict");
const { loadScripts } = require("./helpers/browser_env");

const ctx = loadScripts(["js/scattering_coeff_db.js", "js/scattering_coeff_sources.js"]);
const S = ctx.ScatteringCoeffSource;
const A0 = 0.529177210903; // Bohr radius in Angstrom
const ELEMENTS = [2, 6, 8, 14, 26, 42, 56, 79];
const QS = [1, 2, 4, 6, 8, 12]; // |Q| in 1/Angstrom

function mottBethe(z, q) {
  const s = q / (4 * Math.PI);
  return (z - S.evalXray(z, q, "waasmaier")) / (8 * Math.PI * Math.PI * A0 * s * s);
}

// Elements a table does not cover (Doyle-Turner lacks 44) must raise, not
// scatter as a point charge f_x = Z.
function covered(model, z) {
  return model === "waasmaier" || ctx.ScatteringCoeffDB.scTables[model][z].a.some((x) => x !== 0);
}

test("X-ray scattering factors: f(0) = Z and every table follows Waasmaier-Kirfel", () => {
  for (const model of ["waasmaier", "lobato", "peng_0_4", "doyle", "weickenmeier", "kirkland"]) {
    for (let z = 2; z <= 92; z++) {
      if (!covered(model, z)) {
        assert.throws(() => S.evalXray(z, 0, model), new RegExp(`no parameters for Z=${z};`), `${model} Z=${z}`);
        continue;
      }
      assert.ok(Math.abs(S.evalXray(z, 0, model) - z) / z < 2e-3, `${model} f(0) for Z=${z}: ${S.evalXray(z, 0, model)}`);
    }
    if (model === "waasmaier") continue;
    for (const z of ELEMENTS) for (const q of QS) {
      // The fits differ by a few percent. X-ray values derived from electron
      // fits (Z - 8 pi^2 a0 s^2 f_e) magnify the fit error at large Q, so where
      // f is small the tolerance is 5% of Z. The old s/q and prefactor errors
      // were 13-300% off.
      const f = S.evalXray(z, q, model), ref = S.evalXray(z, q, "waasmaier");
      assert.ok(Math.abs(f - ref) <= Math.max(0.03 * ref, 0.05 * z), `${model} Z=${z} Q=${q}: ${f} vs ${ref}`);
    }
  }
});

test("electron scattering factors follow the Mott-Bethe transform", () => {
  const tables = { lobato: 1, peng_0_4: 2, doyle: 3, weickenmeier: 4, kirkland: 5 };
  for (const [name, num] of Object.entries(tables)) {
    for (const z of ELEMENTS) for (const q of QS) {
      // Weickenmeier-Kohl's He fit is the least accurate here (4.5% at Q = 12),
      // as in Weickenmeier's own code (V = 0.5, B as tabulated)
      const ratio = S.evalElectronNeutral(z, q, num) / mottBethe(z, q);
      assert.ok(Math.abs(ratio - 1) < (name === "weickenmeier" ? 0.06 : 0.03), `${name} Z=${z} Q=${q}: ratio ${ratio}`);
    }
  }
});

test("Weickenmeier-Kohl and Doyle-Turner have no hydrogen; Doyle-Turner lacks Y-Ru and La", () => {
  assert.throws(() => S.evalElectronNeutral(1, 2, 4), /Weickenmeier-Kohl has no parameters for Z=1/);
  assert.throws(() => S.evalXray(1, 2, "weickenmeier"), /Weickenmeier-Kohl has no parameters for Z=1/);
  for (const z of [1, 39, 40, 44, 57]) {
    assert.throws(() => S.evalElectronNeutral(z, 2, 3), new RegExp(`Doyle-Turner has no parameters for Z=${z};`));
    assert.throws(() => S.evalXray(z, 2, "doyle"), new RegExp(`Doyle-Turner has no parameters for Z=${z};`));
  }
  // a vacancy (Z = 0) still scatters nothing
  for (const model of ["lobato", "peng_0_4", "doyle", "kirkland"]) assert.equal(S.evalXray("Va", 2, model), 0);
});

// Every Doyle-Turner row against Waasmaier-Kirfel and its Mott-Bethe transform.
// The table came from MULTEM with twelve transcription errors (Mg a1, Al a4,
// Si b3, K a4, Sc b1, Fe b2, Ga a2, As b2, Ag b2, Cd a2, Sb a3, Rn b1) and Eu's
// row filed under La. The values used now are those two independent copies
// (Felix; pymatgen and EMsoft) agree on, or, where the copies differ, the one
// that follows WK95 and Peng (Sc b1 88.96, Rn b1 38.41). As, Sc and K were
// 8-13% of Z off in X-ray and 16-27% in electron; the Eu row is 6% of Z and
// 12% off as La, 1.4% and 2.3% as Eu. The four-Gaussian fits reach 3% of Z
// near Q = 11 (Li, Na, K, Ca), in every copy.
test("Doyle-Turner rows follow Waasmaier-Kirfel and Mott-Bethe", () => {
  let checked = 0;
  for (let z = 1; z <= 98; z++) {
    if (!covered("doyle", z)) continue;
    checked++;
    for (let q = 0.5; q <= 12; q += 0.5) {
      const ref = S.evalXray(z, q, "waasmaier"), fx = S.evalXray(z, q, "doyle");
      assert.ok(Math.abs(fx - ref) <= 0.035 * z, `Doyle Z=${z} Q=${q}: f_x ${fx} vs ${ref}`);
      const ratio = S.evalElectronNeutral(z, q, 3) / mottBethe(z, q);
      assert.ok(Math.abs(ratio - 1) < 0.04, `Doyle Z=${z} Q=${q}: f_e ratio ${ratio}`);
    }
  }
  assert.equal(checked, 54);
  assert.ok(covered("doyle", 63) && !covered("doyle", 57));
});

// Hydrogen 1s, exact: f_x(Q) = 1 / (1 + (Q a0 / 2)^2)^2, 0.2225 at Q = 4. A table
// evaluated in the wrong variable (Q/4pi for 1/d, or the reverse) gives 0.610
// (Lobato, Kirkland) or 0.03 (Waasmaier, Peng).
function hydrogen(q) {
  return 1 / (1 + (q * A0 / 2) ** 2) ** 2;
}
const H_TABLES = { waasmaier: 0, lobato: 1, peng_0_4: 2, kirkland: 5 };

test("hydrogen follows the exact 1s form factor in every table that has it", () => {
  assert.ok(Math.abs(hydrogen(4) - 0.2225) < 1e-4);
  for (const [name, num] of Object.entries(H_TABLES)) {
    assert.ok(Math.abs(S.evalXray(1, 4, name) - 0.2225) < 1e-3, `${name} f(4) = ${S.evalXray(1, 4, name)}`);
    for (let q = 0.5; q <= 12; q += 0.5) {
      const s = q / (4 * Math.PI), exact = hydrogen(q);
      assert.ok(Math.abs(S.evalXray(1, q, name) - exact) < 0.01, `${name} Q=${q}: ${S.evalXray(1, q, name)} vs ${exact}`);
      if (num === 0) continue;
      const fe = S.evalElectronNeutral(1, q, num), exactFe = (1 - exact) / (8 * Math.PI * Math.PI * A0 * s * s);
      assert.ok(Math.abs(fe / exactFe - 1) < 0.012, `${name} electron Q=${q}: ${fe} vs ${exactFe}`);
    }
  }
});

// Peng (1998) ions: f_e = sum a exp(-b s^2) + 0.023934 dZ / s^2. The inverse
// Mott-Bethe transform gives f_x(0) = Z - dZ, and at large Q the ion and the
// neutral atom scatter alike. H- is left out: its fit is poor beyond Q = 10.
test("Peng ions: f_x(0) = Z - charge and the large-Q limit of the neutral atom", () => {
  const keys = Object.keys(ctx.ScatteringCoeffDB.pengIonByZValence);
  assert.ok(keys.length > 100);
  for (const key of keys) {
    const [z, dz] = key.split(":").map(Number), q0 = 1e-3, s0 = q0 / (4 * Math.PI);
    const fx0 = z - 8 * Math.PI * Math.PI * A0 * s0 * s0 * S.evalElectronIonPeng(z, dz, q0);
    assert.ok(Math.abs(fx0 - (z - dz)) < 1e-3 * z, `${key}: f_x(0) = ${fx0}`);
    if (z === 1) continue;
    for (const q of [6, 8, 10, 12]) {
      const ratio = S.evalElectronIonPeng(z, dz, q) / mottBethe(z, q);
      assert.ok(Math.abs(ratio - 1) < 0.04, `${key} Q=${q}: ratio ${ratio}`);
    }
  }
});

const FULL_SCRIPTS = [
  "js/scattering_coeff_db.js", "js/scattering_coeff_sources.js", "js/diffuse_core.js",
  "js/diffuse_amplitude.js", "js/rmc6f_reader.js",
];

// Exact direct sum in place of a NUFFT: f_k = sum_j c_j exp(i t_k . s_j).
async function directSum(spec) {
  const src = spec.sourcesPacked, t = spec.targetsPacked, c = spec.strengths, out = new Float64Array(t.length / 3 * 2);
  for (let k = 0; k < t.length / 3; k++) for (let j = 0; j < src.length / 3; j++) {
    const p = spec.isign * (t[k * 3] * src[j * 3] + t[k * 3 + 1] * src[j * 3 + 1] + t[k * 3 + 2] * src[j * 3 + 2]);
    out[k * 2] += c[j * 2] * Math.cos(p) - c[j * 2 + 1] * Math.sin(p);
    out[k * 2 + 1] += c[j * 2] * Math.sin(p) + c[j * 2 + 1] * Math.cos(p);
  }
  return { out };
}

// Hydrogen on a cubic cell (a = 3 A), 2 x 2 x 2 supercell, slightly displaced.
function hydrogenCell(full) {
  const rows = [];
  for (let x = 0; x < 2; x++) for (let y = 0; y < 2; y++) for (let z = 0; z < 2; z++) {
    const f = [x + 0.03 * (y - z), y + 0.02 * (x - 0.5), z - 0.04 * x].map((v) => v / 2);
    rows.push(`${rows.length + 1} H [1] ${f.join(" ")} 1 ${x} ${y} ${z}`);
  }
  const text = ["(Version 6f format configuration file)", "Number of atoms: 8", "Supercell dimensions: 2 2 2",
    "Cell (Ang/deg): 6 6 6 90 90 90", "Atoms:", ...rows].join("\n");
  return full.DiffuseCore.attachNeutronCoefficients(full.RMC6fReader.parse("h.rmc6f", text), 10);
}

// End to end: for one species, sqrt(I / I_neutron) * |b| = |f(Q)|. Hydrogen on
// a cubic cell (a = 3 A), Q = 4 1/A on the h axis.
test("one-species intensity over the neutron intensity gives f(Q)/b", async () => {
  const full = loadScripts(FULL_SCRIPTS);
  const parsed = hydrogenCell(full);
  const run = async (scattering) => (await full.DiffuseAmplitude.computeIntensity({
    parsed, h: [4 * 3 / (2 * Math.PI)], k: [0], l: [0], Bq: parsed.Bq, backend: "test", sub: false, runType3: directSum, scattering,
  })).I[0];
  const b = Math.abs(parsed.fca[0]), neutron = await run({ type: "neutron", model: "fast" });
  const exactFe = (1 - 0.22248) / (8 * Math.PI * Math.PI * A0 / (Math.PI * Math.PI)); // s = 1/pi at Q = 4
  for (const [name, num] of Object.entries(H_TABLES)) {
    const fx = b * Math.sqrt((await run({ type: "xray", model: num ? "table" : "waasmaier", tableNum: num })) / neutron);
    assert.ok(Math.abs(fx - 0.2225) < 1e-3, `${name} X-ray f(4) = ${fx}`);
    if (num === 0) continue;
    const fe = b * Math.sqrt((await run({ type: "electron", model: "neutral_table", tableNum: num })) / neutron);
    assert.ok(Math.abs(fe / exactFe - 1) < 0.012, `${name} electron f(4) = ${fe} vs ${exactFe}`);
  }
});

// Published fit ranges: Waasmaier & Kirfel (1995) 0 <= s <= 6 1/A; Peng et al. (1996)
// Table 1 (this table, identical to dynasor's "kmax2" copy; GSAS-II and dynasor's
// "kmax6" hold the s <= 6 fit of Table 3) s <= 2; Doyle & Turner (1968) s <= 2;
// Kirkland (2010) and Lobato & Van Dyck (2014) q = 2s <= 12. Beyond its range a
// Gaussian electron fit falls below the Rutherford tail and its Mott-Bethe X-ray
// value climbs back towards Z: Fe at Q = 30 gave 6.0 (Peng) and 9.6 (Doyle-Turner)
// against 2.37 (Dirac-Fock, arXiv:2608.14934), 2.39 (WK95), 2.44 (Lobato) and 2.38
// (Kirkland); at Q = 40, 15.4 and 19.8 against 1.53. The Peng (1998) ion fits lose
// the tail beyond s = 2 as well (84 of 104 ions by 10% before Q = 40). The evaluators
// keep returning the fits; the intensity calculation warns (RMCProfile's choice too).
const S_MAX = { waasmaier: 6, lobato: 6, peng_0_4: 2, doyle: 2, weickenmeier: 6, kirkland: 6 };
// Dirac-Fock f_x at Q = 30 and 40 (arXiv:2608.14934) for C, Fe and Au.
const DF_FX = { 6: [0.2434, 0.1085], 26: [2.3744, 1.5340], 79: [12.7221, 9.0497] };

test("tables publish their fit range; beyond it the Gaussian fits climb back towards Z", () => {
  const nums = { lobato: 1, peng_0_4: 2, doyle: 3, weickenmeier: 4, kirkland: 5 };
  for (const [name, sMax] of Object.entries(S_MAX)) {
    const r = S.rangeFor("xray", name);
    assert.equal(r.sMax, sMax, name);
    assert.equal(r.qMax, 4 * Math.PI * sMax, name);
    assert.ok(Number.isFinite(S.evalXray(26, r.qMax + 1, name)), name);
    if (name === "waasmaier") continue;
    assert.equal(S.rangeFor("electron", nums[name]).qMax, r.qMax, name);
  }
  assert.equal(S.rangeFor("ion").sMax, 2);
  assert.ok(S.evalXray(26, 30, "peng_0_4") > 5.9 && S.evalXray(26, 30, "doyle") > 9.5);
  assert.ok(S.evalXray(26, 40, "peng_0_4") > 15 && S.evalXray(26, 40, "doyle") > 19);
  // The tables that reach s = 6 follow Dirac-Fock there to 1% of Z.
  for (const [z, ref] of Object.entries(DF_FX)) for (const name of ["waasmaier", "lobato", "weickenmeier", "kirkland"]) {
    [30, 40].forEach((q, i) => {
      const f = S.evalXray(Number(z), q, name);
      assert.ok(Math.abs(f - ref[i]) < 0.01 * z, `${name} Z=${z} Q=${q}: ${f} vs ${ref[i]}`);
    });
  }
});

test("a grid reaching beyond the table's range computes with one warning", async () => {
  const full = loadScripts(FULL_SCRIPTS);
  const parsed = hydrogenCell(full);
  // a = 3 A: the corner (-13, 1, 0) reaches |Q| = 2 pi sqrt(170) / 3 = 27.31 1/A, past
  // s = 2 (25.1) but within s = 6.
  const compute = (scattering, h, k = [0, 1]) => full.DiffuseAmplitude.computeIntensity({
    parsed, h, k, l: [0], Bq: parsed.Bq, backend: "test", sub: false, runType3: directSum, scattering,
  });
  const xray = await compute({ type: "xray", model: "table", tableNum: 2 }, [-13, 0, 5]);
  assert.ok(xray.I.every(Number.isFinite));
  assert.equal(xray.profile.warnings.length, 1);
  assert.match(xray.profile.warnings[0], /\|Q\| = 27\.31.*Peng table \(s <= 2 1\/A, \|Q\| <= 25\.1 1\/A\).*Lobato, Kirkland or Waasmaier-Kirfel/);
  assert.match((await compute({ type: "electron", model: "neutral_table", tableNum: 2 }, [-13, 0, 5])).profile.warnings[0], /Peng table/);
  assert.match((await compute({ type: "electron", model: "ion_peng", valenceMapText: "H:-1" }, [-13, 0, 5])).profile.warnings[0], /Peng ion table/);
  for (const scattering of [{ type: "xray", model: "waasmaier" }, { type: "xray", model: "table", tableNum: 1 }, { type: "electron", model: "neutral_table", tableNum: 5 }]) {
    const res = await compute(scattering, [-13, 0, 5]);
    assert.ok(res.I.every(Number.isFinite));
    assert.deepEqual(Array.from(res.profile.warnings), []);
  }
  // (12, 0, 0) is |Q| = 8 pi, the edge of s = 2: no warning
  assert.deepEqual(Array.from((await compute({ type: "xray", model: "table", tableNum: 2 }, [12], [0])).profile.warnings), []);
});

// f_e(0) = <r^2>/(3 a0) against the Dirac-Fock reference of arXiv:2608.14934 (atomx
// B-spline Dirac-Fock densities; its parameterisation fixes f_e(0) exactly, and H is
// the exact Dirac 1s). The tables rest on two older data sets: Peng, Doyle-Turner and
// Weickenmeier-Kohl on the International Tables (Doyle-Turner, Cromer-Mann), Lobato
// and Kirkland on Rez, Rez & Grant (1994). Up to Hg they agree with the reference to
// 2% (Weickenmeier-Kohl 4.5%); the data sets part for Tl, Pb and Bi (Peng, Doyle-Turner,
// Weickenmeier-Kohl 4-6% low) and for Po, Fr and Ra (Lobato, Kirkland 3-9% low), which
// only matters below |Q| ~ 2.
const DF_FE0 = {
  1: 0.529161, 2: 0.417951, 6: 2.4790, 8: 1.9875, 14: 5.7793, 26: 7.1671, 29: 5.6017, 42: 10.2989, 47: 8.6727,
  56: 18.2749, 79: 10.5786, 81: 12.6120, 82: 13.2651, 84: 13.5260, 87: 18.7262, 88: 20.5725, 92: 19.3034,
};

test("f_e(0) = <r^2>/(3 a0) follows the Dirac-Fock reference", () => {
  const nums = { lobato: 1, peng_0_4: 2, doyle: 3, weickenmeier: 4, kirkland: 5 };
  for (const [zs, ref] of Object.entries(DF_FE0)) {
    const z = Number(zs);
    for (const [name, num] of Object.entries(nums)) {
      if (!covered(name, z)) continue;
      const f = S.evalElectronNeutral(z, 0, num);
      const tol = z > 80 ? 0.1 : name === "weickenmeier" ? 0.045 : 0.02;
      assert.ok(Math.abs(f / ref - 1) < tol, `${name} Z=${z}: f_e(0) ${f} vs ${ref}`);
    }
  }
  const dev = (z, num) => S.evalElectronNeutral(z, 0, num) / DF_FE0[z] - 1;
  for (const z of [81, 82]) for (const num of [1, 5]) assert.ok(Math.abs(dev(z, num)) < 0.02, `Z=${z} table ${num}`);
  for (const z of [84, 87, 88]) assert.ok(Math.abs(dev(z, 2)) < 0.015, `Peng Z=${z}`);
  for (const num of [1, 5]) assert.ok(dev(84, num) < -0.09 && dev(84, num) > -0.1, `Po table ${num}: ${dev(84, num)}`);
});

// Neutron scattering lengths (Sears 1992 / NIST, 1e-12 cm): signs and isotopes, and
// elements without a tabulated natural value (Po, At, Rn, Fr, Ac, Pu, Cm) are reported
// as unknown instead of scattering as b = 0.
test("neutron scattering lengths follow Sears 1992 and elements without one are unknown", () => {
  const core = loadScripts(["js/diffuse_core.js"]).DiffuseCore;
  const sears = { h: -0.3739, d: 0.6671, li: -0.19, "7l": -0.222, ti: -0.3438, v: -0.03824, mn: -0.373, fe: 0.945, o: 0.5803, cl: 0.9577, ra: 1, u: 0.8417 };
  for (const [k, b] of Object.entries(sears)) assert.equal(core.NEUTRON[k], b, k);
  for (const k of ["po", "at", "rn", "fr", "ac", "pu", "cm"]) assert.ok(!(k in core.NEUTRON), k);
  const parsed = core.attachNeutronCoefficients({ elements: ["Po", "Fe", "Xx", "Va"] }, 10);
  assert.deepEqual(Array.from(parsed.unknown), ["Po", "Xx"]);
  assert.deepEqual(Array.from(parsed.fca), [10, 0.945, 10, 0]);
});
