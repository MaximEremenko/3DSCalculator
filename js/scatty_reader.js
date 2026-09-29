"use strict";

// Reads configuration files of J. A. M. Paddison's Scatty:
// [stem]_spins_NN.txt (SPIN lines: magnetic moments on the ideal sites) and
// [stem]_atoms_NN.txt (ATOM lines: displaced atoms). Keywords:
//   TITLE text                CELL a b c alpha beta gamma   BOX Na Nb Nc
//   SITE x y z                (one per site of the unit cell)
//   OCC El occ [El occ ...]   (one line for all sites, or one per site)
//   ATOM site R1 R2 R3 u1 u2 u3 El      (u in unit-cell fractions)
//   SPIN site R1 R2 R3 M1 M2 M3 [El]    (M along a/|a|, b/|b|, c/|c|)
//   FORM_FACTOR_J0 / FORM_FACTOR_J2 A a B b C c D   (one line or one per site)
//   C2 value
// Sites are 1-based and cells R run from 0 to BOX-1.
(function (global) {
  function detect(text) {
    const head = String(text || "").slice(0, 200000);
    return /^\s*BOX\s/m.test(head) && /^\s*SITE\s/m.test(head) && /^\s*(SPIN|ATOM)\s/m.test(head);
  }

  function numbers(tokens, count, label) {
    const out = tokens.slice(0, count).map(Number);
    if (out.length < count || out.some((v) => !Number.isFinite(v))) {
      throw new Error(`Scatty file: invalid ${label} line`);
    }
    return out;
  }

  // Element of the largest occupancy on an OCC line ("El occ El occ ...").
  function majorElement(tokens) {
    let best = null;
    let occ = -Infinity;
    for (let i = 0; i + 1 < tokens.length; i += 2) {
      const v = Number(tokens[i + 1]);
      if (v > occ) { occ = v; best = tokens[i]; }
    }
    return best;
  }

  function parse(name, text) {
    const reader = global.RMC6fReader;
    const ff = global.MagneticFormFactors;
    if (!reader || !ff) throw new Error("Scatty reader needs RMC6fReader and MagneticFormFactors.");
    let cell = null;
    let box = null;
    let title = "";
    let c2 = 0;
    const sites = [];
    const occ = [];
    const entries = [];
    const j0 = [];
    const j2 = [];
    let kind = null;
    for (const raw of String(text || "").split(/\r?\n/)) {
      const t = raw.trim().split(/\s+/);
      if (!t[0]) continue;
      const key = t[0].toUpperCase();
      const rest = t.slice(1);
      if (key === "TITLE") title = raw.trim().slice(5).trim();
      else if (key === "CELL") cell = numbers(rest, 6, "CELL");
      else if (key === "BOX") box = numbers(rest, 3, "BOX").map(Math.round);
      else if (key === "SITE") sites.push(numbers(rest, 3, "SITE"));
      else if (key === "OCC") occ.push(rest);
      else if (key === "FORM_FACTOR_J0") j0.push(numbers(rest, 7, "FORM_FACTOR_J0"));
      else if (key === "FORM_FACTOR_J2") j2.push(numbers(rest, 7, "FORM_FACTOR_J2"));
      else if (key === "C2") c2 = numbers(rest, 1, "C2")[0];
      else if (key === "ATOM" || key === "SPIN") {
        if (kind && kind !== key) throw new Error("Scatty file: ATOM and SPIN lines cannot be mixed.");
        kind = key;
        const v = numbers(rest, 7, key);
        entries.push({ site: Math.round(v[0]), cell: v.slice(1, 4).map(Math.round), vec: v.slice(4, 7), element: rest[7] || null });
      }
    }
    if (!cell) throw new Error("Scatty file: missing CELL");
    if (!box || box.some((v) => !(v > 0))) throw new Error("Scatty file: missing or invalid BOX");
    if (!sites.length) throw new Error("Scatty file: missing SITE lines");
    if (!entries.length) throw new Error("Scatty file: no ATOM or SPIN lines");
    const perSite = (list) => (list.length === 1 ? () => list[0] : (s) => list[s]);
    const occOf = occ.length ? perSite(occ) : () => null;
    const j0Of = j0.length ? perSite(j0) : () => null;
    const j2Of = j2.length ? perSite(j2) : () => null;

    // Form factor and ion of each site (identified from its <j0> coefficients).
    const siteFF = sites.map((_, s) => {
      const c0 = j0Of(s);
      const c2Coeffs = j2Of(s);
      const ion = c0 ? ff.identify(c0) : null;
      return {
        ion,
        formFactor: c0
          ? { j0: ff.coefficients(c0), j2: c2Coeffs ? ff.coefficients(c2Coeffs) : null, c2: c2Coeffs ? c2 : 0 }
          : { j0: [0, 0, 0, 0, 0, 0, 0, 0, 1], j2: null, c2: 0 },
      };
    });

    const n = entries.length;
    const fractional = new Float64Array(n * 3);
    const cellIndex = new Int32Array(n * 3);
    const elements = new Array(n);
    entries.forEach((e, i) => {
      const s = e.site - 1;
      if (!(s >= 0 && s < sites.length)) throw new Error(`Scatty file: site ${e.site} is not defined`);
      const u = kind === "ATOM" ? e.vec : [0, 0, 0];
      for (let d = 0; d < 3; d++) {
        cellIndex[i * 3 + d] = e.cell[d];
        fractional[i * 3 + d] = (e.cell[d] + sites[s][d] + u[d]) / box[d];
      }
      const occLine = occOf(s);
      const ion = siteFF[s].ion;
      elements[i] = e.element || (occLine && majorElement(occLine)) || (ion ? ion.ion.replace(/\d+$/, "") : "X");
    });

    const structure = reader.buildSupercellStructure({ parentCellDeg: cell, supercell: box, fractional, cellIndex, elements });
    const parsed = {
      file: String(name || ""),
      ...structure,
      cellDeg: [cell[0] * box[0], cell[1] * box[1], cell[2] * box[2], cell[3], cell[4], cell[5]],
      sourceFormat: kind === "SPIN" ? "Scatty spin configuration" : "Scatty atom configuration",
      title,
    };

    if (kind === "SPIN") {
      // Unit vectors along a, b, c in the same Cartesian frame as the atoms.
      const axes = reader.cellGeometry(cell).direct.map((r) => {
        const len = Math.hypot(r[0], r[1], r[2]);
        return [r[0] / len, r[1] / len, r[2] / len];
      });
      const moments = new Float64Array(n * 3);
      entries.forEach((e, i) => {
        for (let d = 0; d < 3; d++) {
          moments[i * 3 + d] = e.vec[0] * axes[0][d] + e.vec[1] * axes[1][d] + e.vec[2] * axes[2][d];
        }
      });
      // One magnetic species per distinct form factor.
      const species = [];
      const keys = new Map();
      const siteSpecies = siteFF.map((sf) => {
        const key = JSON.stringify(sf.formFactor);
        if (!keys.has(key)) {
          keys.set(key, species.length);
          species.push({
            label: sf.ion ? `${sf.ion.ion.replace(/(\d+)$/, "$1+")}` : "custom",
            ion: sf.ion ? sf.ion.ion : null,
            formFactor: sf.formFactor,
          });
        }
        return keys.get(key);
      });
      const speciesOfAtom = new Int32Array(n);
      entries.forEach((e, i) => { speciesOfAtom[i] = siteSpecies[e.site - 1]; });
      parsed.magnetic = { species, speciesOfAtom, moments, source: "Scatty SPIN (components along a, b, c)" };
      // Without OCC lines Scatty treats the file as magnetic-only.
      parsed.magneticOnly = !occ.length && !entries.some((e) => e.element);
    }
    return parsed;
  }

  global.ScattyReader = Object.freeze({ detect, parse });
})(typeof window !== "undefined" ? window : globalThis);
