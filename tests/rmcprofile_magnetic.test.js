"use strict";

const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const { loadScripts } = require("./helpers/browser_env");

const ctx = loadScripts([
  "js/rmc6f_reader.js", "js/lammps_data_reader.js", "js/magnetic_form_factors.js", "js/diffuse_core.js",
  "js/diffuse_amplitude.js", "js/scatty_reader.js", "js/rmcprofile_magnetic.js",
]);
const core = ctx.DiffuseCore;
const RM = ctx.RmcMagnetic;
const plain = (v) => JSON.parse(JSON.stringify(v));

// The MAGNETISM block of RMCProfile's GdRu2Si2 fixture
// (tests/fixtures/delta_pdf3d/gdru2si2_finite8/combined_current_configuration.dat).
const GDRU2SI2_DAT = `ATOMS :: Gd Ru Si
IGNORE_HISTORY_FILE ::

FLAGS ::
  > NO_MOVEOUT

MAGNETISM :: YES
  > MAGNETISM_FILE_STEM :: gdru2si2_weak_disorder_skyrmion_16000_spin
  > MAGNETIC_ATOMS :: 1
  > FORM_FACTOR :: 1 0.0186 25.3867 0.2895 11.1421 0.7135 3.752 -0.0217
  > MAGNETIC_MOMENTS :: 7.000
  > MAX_SPIN_MOVEMENT :: 0.0
  > NO_VARY_SPIN_MOVE_RATE

END ::
`;

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

function magneticIntensity(parsed, axes) {
  return ctx.DiffuseAmplitude.computeIntensity({
    parsed, h: axes, k: axes, l: axes, Bq: parsed.Bq, backend: "test", sub: false,
    scattering: { type: "neutron", model: "fast" }, runType3: directSum,
    magnetic: { mode: "magnetic", qZero: "average" },
  });
}

function maxRelDiff(a, b) {
  let peak = 0, diff = 0;
  for (let i = 0; i < a.length; i++) {
    peak = Math.max(peak, Math.abs(b[i]));
    diff = Math.max(diff, Math.abs(a[i] - b[i]));
  }
  return diff / peak;
}

// Ordered MnO from the Scatty fixture: 864 Mn spins of 5 muB in a 6^3 box.
const SCATTY = core.attachNeutronCoefficients(ctx.ScattyReader.parse("mno_order_spins_01.txt",
  fs.readFileSync(path.join(__dirname, "fixtures", "magnetic", "mno_order_spins_01.txt"), "utf8")), 10);
const A = 4.4344, BOX = 6;
const AXIS = Array.from({ length: 7 }, (_, i) => Number((-0.5 + i / 6).toFixed(12)));

// rocksalt MnO as .rmc6f (Mn first, then O), spins of the Scatty fixture as unit vectors
function mnoFiles(rotation) {
  const mn = [], o = [];
  for (let cx = 0; cx < BOX; cx++) for (let cy = 0; cy < BOX; cy++) for (let cz = 0; cz < BOX; cz++) {
    for (const s of [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]) mn.push([cx, cy, cz, s]);
    for (const s of [[0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5], [0.5, 0.5, 0.5]]) o.push([cx, cy, cz, s]);
  }
  const rows = [...mn.map((a) => ["Mn", a]), ...o.map((a) => ["O", a])].map(([el, [cx, cy, cz, s]], i) =>
    `${i + 1} ${el} [1] ${[cx + s[0], cy + s[1], cz + s[2]].map((v) => (v / BOX).toFixed(15)).join(" ")} 1 ${cx} ${cy} ${cz}`);
  const L = BOX * A;
  let lattice = [[L, 0, 0], [0, L, 0], [0, 0, L]];
  if (rotation) lattice = lattice.map((r) => [0, 1, 2].map((j) => r[0] * rotation[0][j] + r[1] * rotation[1][j] + r[2] * rotation[2][j]));
  const rmc6f = [
    "(Version 6f format configuration file)", "Atom types present:  Mn O", "Number of each atom type:  864 864",
    `Number of atoms: ${rows.length}`, `Supercell dimensions: ${BOX} ${BOX} ${BOX}`, `Cell (Ang/deg): ${L} ${L} ${L} 90 90 90`,
    ...(rotation ? ["Lattice vectors (Ang):", ...lattice.map((r) => r.map((v) => v.toFixed(12)).join(" "))] : []),
    "Atoms:", ...rows,
  ].join("\n");
  // spins in Scatty order (cell-major, 4 sites) = Mn order above; unit vectors in the file frame
  const spins = [];
  for (let i = 0; i < 864; i++) {
    const v = Array.from(SCATTY.magnetic.vectors.subarray(i * 3, i * 3 + 3)).map((x) => x / 5);
    spins.push(rotation ? [0, 1, 2].map((j) => v[0] * rotation[0][j] + v[1] * rotation[1][j] + v[2] * rotation[2][j]) : v);
  }
  const legacy = ["0.0917 864 1", "864", ...spins.map((s) => s.map((v) => v.toFixed(15)).join(" "))].join("\n");
  const v2 = [RM.SIGNATURE, "SPIN_MODEL HEISENBERG", "CONFIGURATION 1728 2", "TYPE_POPULATIONS 864 864", "TYPE_NAMES Mn O",
    "NUMBER_DENSITY 0.0917", "MAGNETIC_SPECIES 1 864", "SPIN_POPULATIONS 864", "SPINS 864",
    ...spins.map((s, j) => `SPIN ${j + 1} 1 ${j + 1} ${s.map((v) => v.toFixed(15)).join(" ")}`), "END_MAGNETIC_CONFIGURATION"].join("\n");
  const dat = "ATOMS :: Mn O\nMAGNETISM ::\n  > MAGNETIC_ATOMS :: 1\n  > FORM_FACTOR :: 1 0.4220 17.6840 0.5948 6.0050 0.0043 -0.6090 -0.0219\n  > MAGNETIC_MOMENTS :: 5.0\nEND ::\n";
  return { rmc6f, legacy, v2, dat };
}

function loadRmc(files, spinText) {
  const parsed = core.attachNeutronCoefficients(ctx.RMC6fReader.parse("mno.rmc6f", files.rmc6f), 10);
  parsed.magnetic = RM.attachSpins(parsed, RM.parseSpinConfiguration(spinText), RM.parseMagnetismBlock(files.dat));
  return parsed;
}

test("parses the MAGNETISM block of RMCProfile's GdRu2Si2 fixture", () => {
  const m = RM.parseMagnetismBlock(GDRU2SI2_DAT);
  assert.deepEqual(plain(m.types), ["Gd", "Ru", "Si"]);
  assert.equal(m.magneticAtoms, 1);
  assert.deepEqual(plain(m.moments), [7]);
  assert.deepEqual(plain(m.formFactors), [{ index: 1, coefficients: [0.0186, 25.3867, 0.2895, 11.1421, 0.7135, 3.752, -0.0217] }]);
  assert.equal(m.stem, "gdru2si2_weak_disorder_skyrmion_16000_spin");
  assert.equal(m.mag5d, false);
  assert.equal(RM.parseMagnetismBlock("ATOMS :: Mn O\nMAGNETISM :: NO\n"), null);
  assert.equal(RM.parseMagnetismBlock("ATOMS :: Mn O\nFLAGS ::\n  > MAG_5D\nMAGNETISM :: YES\n").mag5d, true);
});

test("detects legacy and version-2 spin configurations", () => {
  const f = mnoFiles();
  assert.equal(RM.detectSpinConfiguration(f.legacy), "legacy");
  assert.equal(RM.detectSpinConfiguration(f.v2), "v2");
  assert.equal(RM.detectSpinConfiguration(f.rmc6f), null);
  assert.equal(RM.detectSpinConfiguration(f.dat), null);
});

test("RMCProfile rmc6f + spin cfg + MAGNETISM reproduce the Scatty-validated MnO intensities", async () => {
  const files = mnoFiles();
  const parsed = loadRmc(files, files.legacy);
  assert.equal(parsed.magnetic.species[0].scale, 5);
  assert.equal(parsed.magnetic.species[0].formFactorSource, "dat");
  const ours = await magneticIntensity(parsed, AXIS);
  const reference = await magneticIntensity(SCATTY, AXIS);
  assert.ok(maxRelDiff(ours.I, reference.I) < 1e-9);
});

test("version-2 spin files and rotated lattice-vector frames give the same intensities", async () => {
  const reference = await magneticIntensity(SCATTY, AXIS);
  const plainFiles = mnoFiles();
  const v2 = await magneticIntensity(loadRmc(plainFiles, plainFiles.v2), AXIS);
  assert.ok(maxRelDiff(v2.I, reference.I) < 1e-9);
  // file frame rotated by 37 degrees about (1,2,3): spins are rotated back
  const n = [1, 2, 3].map((v) => v / Math.sqrt(14)), c = Math.cos(0.645772), s = Math.sin(0.645772);
  const R = [0, 1, 2].map((i) => [0, 1, 2].map((j) =>
    (i === j ? c : 0) + (1 - c) * n[i] * n[j] + s * [[0, n[2], -n[1]], [-n[2], 0, n[0]], [n[1], -n[0], 0]][i][j]));
  const rotated = mnoFiles(R);
  const parsed = loadRmc(rotated, rotated.legacy);
  assert.match(parsed.magnetic.source, /rotated/);
  const out = await magneticIntensity(parsed, AXIS);
  assert.ok(maxRelDiff(out.I, reference.I) < 1e-9);
});

test("spin files that do not fit the structure are refused", () => {
  const files = mnoFiles();
  const parsed = core.attachNeutronCoefficients(ctx.RMC6fReader.parse("mno.rmc6f", files.rmc6f), 10);
  const tooMany = ["0.09 900 1", "900", ...Array(900).fill("0 0 1")].join("\n");
  assert.throws(() => RM.attachSpins(parsed, RM.parseSpinConfiguration(tooMany), null), /only 864 Mn atoms/);
  const wrongAtom = files.v2.replace("SPIN 1 1 1 ", "SPIN 1000 1 1 ");
  assert.throws(() => RM.attachSpins(parsed, RM.parseSpinConfiguration(wrongAtom), null), /on a O atom/);
});

test("LAMMPS atom_style spin gives the same MnO intensities", async () => {
  const L = BOX * A;
  const rows = [];
  for (let i = 0; i < SCATTY.atoms; i++) {
    const m = Array.from(SCATTY.magnetic.vectors.subarray(i * 3, i * 3 + 3));
    const sp = Math.hypot(...m);
    rows.push(`${i + 1} 1 ${SCATTY.x[i]} ${SCATTY.y[i]} ${SCATTY.z[i]} ${m.map((v) => v / sp).join(" ")} ${sp}`);
  }
  const data = [
    "LAMMPS data file (MnO spins)", "", `${SCATTY.atoms} atoms`, "1 atom types", "",
    `0 ${L} xlo xhi`, `0 ${L} ylo yhi`, `0 ${L} zlo zhi`, "", "Masses", "", "1 54.938 # Mn", "", "Atoms # spin", "", ...rows, "",
  ].join("\n");
  const parsed = core.attachNeutronCoefficients(ctx.LammpsDataReader.parse("mno.data", data), 10);
  assert.equal(parsed.magnetic.species[0].ion, "Mn2");
  parsed.magnetic.species[0].formFactor = SCATTY.magnetic.species[0].formFactor;
  // LAMMPS hkl refer to the whole box: h_box = 6 h
  const boxAxis = AXIS.map((v) => Number((v * BOX).toFixed(12)));
  const ours = await magneticIntensity(parsed, boxAxis);
  const reference = await magneticIntensity(SCATTY, AXIS);
  assert.ok(maxRelDiff(ours.I, reference.I) < 1e-9);
});
