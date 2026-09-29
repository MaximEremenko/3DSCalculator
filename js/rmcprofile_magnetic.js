"use strict";

// RMCProfile magnetic input for an already loaded .rmc6f structure:
// spin configurations (<stem>.cfg) in the legacy form
//   rho  total_spins  n_spin_types
//   m_1 ... m_n
//   sx sy sz                 (one unit vector per spin, species by species)
// and the "RMCProfile magnetic configuration version 2" form (SPIN records
// naming the atom), plus the MAGNETISM block of the .dat control file
// (MAGNETIC_ATOMS, MAGNETIC_MOMENTS in muB per species, FORM_FACTOR
// records with 7 coefficients, or 9 with MAG_5D). Magnetic species t is atom
// type t in RMCProfile's type order; in legacy files spin j of species t
// belongs to the j-th atom of that type. Spins are Cartesian in RMCProfile's
// frame: the rmc6f "Lattice vectors" when present, otherwise c along z and b
// in the yz plane, which is also this calculator's frame.
(function (global) {
  const SIGNATURE = "RMCProfile magnetic configuration version 2";

  function nonEmptyLines(text) {
    return String(text || "").split(/\r?\n/).map((l) => l.trim()).filter(Boolean);
  }

  // "legacy", "v2" or null.
  function detectSpinConfiguration(text) {
    const lines = nonEmptyLines(String(text || "").slice(0, 4096));
    if (!lines.length) return null;
    if (lines[0] === SIGNATURE) return "v2";
    const head = lines[0].split(/\s+/).map(Number);
    if (head.length !== 3 || !head.every(Number.isFinite)) return null;
    const [, total, types] = head;
    if (!(Number.isInteger(total) && total > 0 && Number.isInteger(types) && types > 0)) return null;
    const counts = (lines[1] || "").split(/\s+/).map(Number);
    if (counts.length !== types || !counts.every((c) => Number.isInteger(c) && c >= 0)) return null;
    return counts.reduce((a, b) => a + b, 0) === total ? "legacy" : null;
  }

  function vector(tokens, where) {
    const v = tokens.slice(0, 3).map(Number);
    if (v.length < 3 || !v.every(Number.isFinite)) throw new Error(`Spin configuration: invalid spin ${where}`);
    return v;
  }

  // { format, counts, spins (3 per spin, species-major), atoms (1-based, v2 only), typeNames }
  function parseSpinConfiguration(text) {
    const format = detectSpinConfiguration(text);
    if (!format) throw new Error("Not an RMCProfile spin configuration.");
    const lines = nonEmptyLines(text);
    if (format === "legacy") {
      const total = Number(lines[0].split(/\s+/)[1]);
      const counts = lines[1].split(/\s+/).map(Number);
      if (lines.length < total + 2) throw new Error(`Spin configuration: ${total} spins declared, ${lines.length - 2} found`);
      const spins = new Float64Array(total * 3);
      for (let i = 0; i < total; i++) spins.set(vector(lines[i + 2].split(/\s+/), `on line ${i + 3}`), i * 3);
      return { format, counts, spins, atoms: null, typeNames: null };
    }
    let counts = null;
    let typeNames = null;
    let total = null;
    const records = [];
    for (const line of lines.slice(1)) {
      const t = line.split(/\s+/);
      if (t[0] === "SPIN_POPULATIONS") counts = t.slice(1).map(Number);
      else if (t[0] === "TYPE_NAMES") typeNames = t.slice(1);
      else if (t[0] === "SPINS") total = Number(t[1]);
      else if (t[0] === "SPIN") records.push(t.slice(1));
      else if (t[0] === "END_MAGNETIC_CONFIGURATION") break;
    }
    if (!counts || records.length !== total) {
      throw new Error("Spin configuration (version 2): the SPIN records do not match SPINS / SPIN_POPULATIONS");
    }
    const offsets = counts.map((_, t) => counts.slice(0, t).reduce((a, b) => a + b, 0));
    const spins = new Float64Array(total * 3);
    const atoms = new Int32Array(total);
    records.forEach((r, n) => {
      const [atom, type, local] = r.slice(0, 3).map(Number);
      if (!(type >= 1 && type <= counts.length && local >= 1 && local <= counts[type - 1])) {
        throw new Error(`Spin configuration (version 2): SPIN record ${n + 1} has an invalid species slot`);
      }
      const i = offsets[type - 1] + local - 1;
      atoms[i] = atom;
      spins.set(vector(r.slice(3), `in SPIN record ${n + 1}`), i * 3);
    });
    return { format, counts, spins, atoms, typeNames };
  }

  // The MAGNETISM block (and ATOMS, MAG_5D) of an RMCProfile .dat file, or
  // null when the file enables no magnetism.
  function parseMagnetismBlock(text) {
    const out = { types: null, magneticAtoms: null, moments: null, formFactors: [], mag5d: false, stem: null, formFactorFile: null };
    let block = null;
    let enabled = false;
    for (const raw of String(text || "").split(/\r?\n/)) {
      const line = raw.replace(/[!#].*$/, "").trim();
      if (!line) continue;
      if (line.startsWith(">")) {
        const m = /^>\s*([A-Za-z0-9_]+)\s*(?:::\s*(.*))?$/.exec(line);
        if (!m) continue;
        const key = m[1].toUpperCase();
        const value = (m[2] || "").trim();
        if (key === "MAG_5D" && (block === "MAGNETISM" || block === "FLAGS")) out.mag5d = true;
        if (block !== "MAGNETISM") continue;
        const tokens = value.split(/\s+/).filter(Boolean);
        if (key === "MAGNETIC_ATOMS") out.magneticAtoms = Number(tokens[0]);
        else if (key === "MAGNETIC_MOMENTS") out.moments = tokens.map(Number);
        else if (key === "MAGNETISM_FILE_STEM") out.stem = value;
        else if (key === "MAGNETIC_FORM_FACTOR_FILE" || key === "MAGNETIC_BRAGG_FORM_FACTOR_FILE") out.formFactorFile = value;
        else if (key === "FORM_FACTOR") {
          const numbered = /^\d+$/.test(tokens[0]);
          out.formFactors.push({
            index: numbered ? Number(tokens[0]) : out.formFactors.length + 1,
            coefficients: tokens.slice(1).map(Number),
          });
        }
        continue;
      }
      const m = /^([A-Za-z0-9_]+)\s*::\s*(.*)$/.exec(line);
      if (m) {
        block = m[1].toUpperCase();
        if (block === "ATOMS") out.types = m[2].trim().split(/\s+/).filter(Boolean);
        if (block === "MAGNETISM") enabled = !/^NO\b/i.test(m[2].trim());
      } else if (/^END\b/i.test(line)) {
        block = null;
      }
    }
    return enabled ? out : null;
  }

  function inv3(m) {
    const [a, b, c] = m[0];
    const [d, e, f] = m[1];
    const [g, h, i] = m[2];
    const A = e * i - f * h, B = -(d * i - f * g), C = d * h - e * g;
    const det = a * A + b * B + c * C;
    if (!Number.isFinite(det) || Math.abs(det) < 1e-20) return null;
    const s = 1 / det;
    return [
      [A * s, -(b * i - c * h) * s, (b * f - c * e) * s],
      [B * s, (a * i - c * g) * s, -(a * f - c * d) * s],
      [C * s, -(a * h - b * g) * s, (a * e - b * d) * s],
    ];
  }

  function matmul(x, y) {
    return x.map((row) => [0, 1, 2].map((j) => row[0] * y[0][j] + row[1] * y[1][j] + row[2] * y[2][j]));
  }

  // Rotation M (row vectors: s_ours = det(M) s_file M) from the file's
  // lattice-vector frame into this calculator's frame, or null when the file
  // gives no lattice vectors. `direct` holds the parent-cell rows in our frame.
  function frameRotation(fileLattice, direct, supercell) {
    if (!fileLattice) return null;
    const inverse = inv3(fileLattice);
    if (!inverse) throw new Error("rmc6f lattice vectors are singular");
    for (const scale of [supercell, [1, 1, 1]]) {
      const ours = direct.map((row, i) => row.map((v) => v * scale[i]));
      const M = matmul(inverse, ours);
      const MtM = [0, 1, 2].map((i) => [0, 1, 2].map((j) => M[0][i] * M[0][j] + M[1][i] * M[1][j] + M[2][i] * M[2][j]));
      if (MtM.every((row, i) => row.every((v, j) => Math.abs(v - (i === j ? 1 : 0)) < 1e-6))) return M;
    }
    throw new Error("rmc6f lattice vectors do not match the Cell line (the spin frame is ambiguous)");
  }

  function rotateSpin(s, M) {
    if (!M) return s;
    const det = M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
      - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
      + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
    const sign = det < 0 ? -1 : 1; // moments are axial vectors
    return [0, 1, 2].map((j) => sign * (s[0] * M[0][j] + s[1] * M[1][j] + s[2] * M[2][j]));
  }

  function formFactorFromRecord(record, mag5d) {
    const ff = global.MagneticFormFactors;
    const c = record.coefficients.filter(Number.isFinite);
    if (c.length !== (mag5d ? 9 : 7)) {
      throw new Error(`FORM_FACTOR ${record.index}: expected ${mag5d ? 9 : 7} coefficients, found ${c.length}`);
    }
    return { j0: ff.coefficients(c), j2: null, c2: 0 };
  }

  // Moments and form factors of the MAGNETISM block applied to the species of
  // an attached spin configuration.
  function applyMagnetism(magnetic, magnetism) {
    if (!magnetic || !magnetism) return [];
    const notes = [];
    magnetic.species.forEach((sp, t) => {
      if (magnetism.moments && Number.isFinite(magnetism.moments[t])) sp.scale = magnetism.moments[t];
      const record = magnetism.formFactors.find((r) => r.index === t + 1);
      if (record) {
        sp.formFactor = formFactorFromRecord(record, magnetism.mag5d);
        sp.fileFormFactor = sp.formFactor;
        sp.ion = null;
        sp.formFactorSource = "dat";
      }
    });
    if (magnetism.formFactorFile) notes.push("tabulated MAGNETIC_FORM_FACTOR_FILE is not read; the ion table is used instead");
    return notes;
  }

  // Attaches a parsed spin configuration to a structure parsed from .rmc6f.
  // typeOrder: RMCProfile's atom-type order (.dat ATOMS, else the spin file's
  // TYPE_NAMES, else the rmc6f "Atom types present", else first appearance).
  function attachSpins(parsed, cfg, magnetism) {
    const reader = global.RMC6fReader;
    const ffTable = global.MagneticFormFactors;
    const firstSeen = [];
    for (const e of parsed.elements) if (!firstSeen.includes(e)) firstSeen.push(e);
    const typeOrder = (magnetism && magnetism.types) || cfg.typeNames || parsed.typeOrder || firstSeen;
    const types = cfg.counts.length;
    if (typeOrder.length < types) throw new Error(`The spin file has ${types} magnetic species but only ${typeOrder.length} atom types are known`);
    const same = (a, b) => String(a).toLowerCase() === String(b).toLowerCase();
    const members = typeOrder.slice(0, types).map((el) => {
      const out = [];
      parsed.elements.forEach((e, i) => { if (same(e, el)) out.push(i); });
      return out;
    });
    const supercell = parsed.super.map((v) => Math.max(1, Math.round(Number(v) || 1)));
    const parent = [0, 1, 2].map((d) => parsed.cellDeg[d] / supercell[d]).concat(parsed.cellDeg.slice(3, 6));
    const found = frameRotation(parsed.fileLattice, reader.cellGeometry(parent).direct, supercell);
    const M = found && found.some((row, i) => row.some((v, j) => Math.abs(v - (i === j ? 1 : 0)) > 1e-9)) ? found : null;
    const speciesOfAtom = new Int32Array(parsed.atoms).fill(-1);
    const vectors = new Float64Array(parsed.atoms * 3);
    let worst = 0;
    let k = 0;
    for (let t = 0; t < types; t++) {
      if (cfg.counts[t] > members[t].length) {
        throw new Error(`${cfg.counts[t]} spins for ${typeOrder[t]} but only ${members[t].length} ${typeOrder[t]} atoms`);
      }
      for (let j = 0; j < cfg.counts[t]; j++, k++) {
        const atom = cfg.atoms ? cfg.atoms[k] - 1 : members[t][j];
        if (!(atom >= 0 && atom < parsed.atoms)) throw new Error(`Spin ${k + 1} names atom ${atom + 1}, outside the configuration`);
        if (!same(parsed.elements[atom], typeOrder[t])) {
          throw new Error(`Spin ${k + 1} (species ${typeOrder[t]}) is on a ${parsed.elements[atom]} atom`);
        }
        if (speciesOfAtom[atom] >= 0) throw new Error(`Atom ${atom + 1} carries two spins`);
        const s = rotateSpin(Array.from(cfg.spins.subarray(k * 3, k * 3 + 3)), M);
        worst = Math.max(worst, Math.abs(Math.hypot(s[0], s[1], s[2]) - 1));
        speciesOfAtom[atom] = t;
        vectors.set(s, atom * 3);
      }
    }
    const species = typeOrder.slice(0, types).map((el) => {
      const ion = (/\d/.test(el) && ffTable.find(el)) || ffTable.defaultIon(el);
      return {
        label: el,
        element: el,
        ion: ion ? ion.ion : null,
        scale: 1,
        formFactor: ion ? { j0: ion.j0, j2: ion.j2, c2: 0 } : { j0: [0, 0, 0, 0, 0, 0, 0, 0, 1], j2: null, c2: 0 },
        formFactorSource: ion ? "table" : "none",
      };
    });
    const magnetic = {
      species,
      speciesOfAtom,
      vectors,
      source: `RMCProfile spin configuration (${cfg.format === "v2" ? "version 2" : "legacy"}${M ? ", rotated from the rmc6f lattice-vector frame" : ""})`,
      notes: worst > 1e-3 ? [`spin vectors deviate from unit length by up to ${worst.toExponential(1)}`] : [],
    };
    magnetic.notes.push(...applyMagnetism(magnetic, magnetism));
    return magnetic;
  }

  global.RmcMagnetic = Object.freeze({
    SIGNATURE,
    detectSpinConfiguration,
    parseSpinConfiguration,
    parseMagnetismBlock,
    frameRotation,
    applyMagnetism,
    attachSpins,
  });
})(typeof window !== "undefined" ? window : globalThis);
