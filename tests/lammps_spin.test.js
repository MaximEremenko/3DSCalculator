"use strict";

const test = require("node:test");
const assert = require("node:assert/strict");
const { loadScripts } = require("./helpers/browser_env");

const ctx = loadScripts(["js/lammps_data_reader.js", "js/magnetic_form_factors.js"]);
const L = ctx.LammpsDataReader;

function dataFile(header, rows) {
  return [header, "", `${rows.length} atoms`, "1 atom types", "", "0 4 xlo xhi", "0 4 ylo yhi", "0 4 zlo zhi", "",
    "Masses", "", "1 55.845", "", "Atoms # spin", "", ...rows, ""].join("\n");
}
const moments = (p) => Array.from({ length: p.atoms }, (_, i) => Array.from(p.magnetic.vectors.subarray(i * 3, i * 3 + 3)).map((v) => +v.toFixed(12)));
const positions = (p) => Array.from({ length: p.atoms }, (_, i) => [p.x[i], p.y[i], p.z[i]].map((v) => +v.toFixed(12)));

test("reads the current LAMMPS spin layout: id type x y z spx spy spz sp", () => {
  const p = L.parse("new.data", dataFile("LAMMPS data file via write_data, version 4 Feb 2025", [
    "1 1 0.0 0.0 0.0 1 0 0 2.2 0 0 0",
    "2 1 1.5 1.5 1.5 0 0 2 1.5",
  ]));
  assert.deepEqual(positions(p), [[0, 0, 0], [1.5, 1.5, 1.5]]);
  assert.deepEqual(moments(p), [[2.2, 0, 0], [0, 0, 1.5]]); // direction normalized as LAMMPS does
  assert.equal(p.sourceWarnings.some((w) => /pre-2020/.test(w)), false);
});

test("reads the pre-2020 layout: id type sp x y z spx spy spz", () => {
  const p = L.parse("old.data", dataFile("LAMMPS data file via write_data, version 4 Jan 2019, timestep = 0", [
    "1 1 2.2 0.0 0.0 0.0 1 0 0 0 0 0",
    "2 1 2.2 1.43 1.43 1.43 1 0 0 0 0 0",
  ]));
  assert.deepEqual(positions(p), [[0, 0, 0], [1.43, 1.43, 1.43]]);
  assert.deepEqual(moments(p), [[2.2, 0, 0], [2.2, 0, 0]]);
  assert.ok(p.sourceWarnings.some((w) => /pre-2020 LAMMPS layout/.test(w)));
});

test("unit directions decide the layout; the write_data version breaks ties", () => {
  // both triplets are unit vectors in this row
  const row = ["1 1 1 0 0 0 1 0 0"];
  const old = L.parse("a", dataFile("LAMMPS data file via write_data, version 28 Feb 2019", row));
  assert.deepEqual(positions(old), [[0, 0, 0]]);
  assert.deepEqual(moments(old), [[1, 0, 0]]);
  const cur = L.parse("b", dataFile("LAMMPS data file via write_data, version 29 Aug 2024", row));
  assert.deepEqual(positions(cur), [[1, 0, 0]]);
  assert.equal(cur.magnetic, null); // sp = 0 in the current layout: no moments
  // no version: the rows decide (legacy direction in columns 7-9)
  const rows = ["1 1 2.5 0 0 0 0.6 0.8 0", "2 1 2.5 1 1 0.3 0 0.6 0.8"];
  assert.deepEqual(moments(L.parse("c", dataFile("hand-written", rows))), [[1.5, 2, 0], [0, 1.5, 2]]);
});
