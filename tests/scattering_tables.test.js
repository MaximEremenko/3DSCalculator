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

test("X-ray scattering factors: f(0) = Z and every table follows Waasmaier-Kirfel", () => {
  for (const model of ["waasmaier", "lobato", "peng_0_4", "doyle", "weickenmeier", "kirkland"]) {
    for (let z = 2; z <= 92; z++) {
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
      // Doyle-Turner (four Gaussians) and Weickenmeier-Kohl are the least
      // accurate fits here (5% for Si at Q = 8, 4.5% for He at Q = 12)
      const ratio = S.evalElectronNeutral(z, q, num) / mottBethe(z, q);
      assert.ok(Math.abs(ratio - 1) < (name === "doyle" || name === "weickenmeier" ? 0.06 : 0.03), `${name} Z=${z} Q=${q}: ratio ${ratio}`);
    }
  }
});

test("Weickenmeier-Kohl has no hydrogen", () => {
  assert.throws(() => S.evalElectronNeutral(1, 2, 4), /no parameters for Z=1/);
  assert.throws(() => S.evalXray(1, 2, "weickenmeier"), /no parameters for Z=1/);
});
