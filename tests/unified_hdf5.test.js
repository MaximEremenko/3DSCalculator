"use strict";

const test = require("node:test");
const assert = require("node:assert/strict");
const { loadH5, writeEntryData } = require("./helpers/browser_env");

const plain = (value) => JSON.parse(JSON.stringify(value));
const str = (text) => ({ data: [text], shape: [1], dtype: "S" + Math.max(1, text.length) });
const scalarStr = (text) => ({ data: text, shape: [], dtype: "S" + Math.max(1, text.length) });
const i32 = (values, shape = [values.length]) => ({ data: Int32Array.from(values), shape, dtype: "<i" });
const f64 = (values, shape = [values.length]) => ({ data: Float64Array.from(values), shape, dtype: "<d" });
const u64 = (values) => ({ data: BigUint64Array.from(values.map(BigInt)), shape: [values.length], dtype: "<Q" });

// Layout written by RMCProfile 6.8.0-rc.1 rmc6f_to_unified.exe: parent cell,
// positions in unit-cell units, 1-based (N,3) cells, ';'-joined type names.
function rmcProfile68(overrides = {}) {
  return {
    audit_conform_dict_name: str("Disorder structure"),
    unit_cell_lengths: f64([4, 9, 10]),
    unit_cell_angles: f64([90, 100, 110]),
    unit_cells: i32([2, 1, 1]),
    number_of_atoms: i32([4]),
    number_of_types: i32([2]),
    types_names: str("Na;Cl"),
    atom_type: i32([1, 1, 2, 2]),
    atom_position: f64([0.25, 0.25, 0.375, 1.25, 0.75, 0.875, 0.2, 0.2, 0.3, 1.2, 0.7, 0.8], [4, 3]),
    atom_unit_cell: i32([1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1], [4, 3]),
    coordinates_are_supercell_fractional: i32([0]),
    ...overrides,
  };
}

// Layout written by DISCUS / the H5FORTRAN reference: (3,N) arrays, scalar
// strings, uint64 counts and optional groups such as magnetic_spins.
function discusLayout() {
  const x = [0.0, 0.13, 1.13, 1.5, 0.5], y = [0.1, 0.0, 0.0, 0.5, 0.5], z = [0.2, 0.0, 0.0, 0.5, 0.0];
  return {
    audit_conform_dict_name: scalarStr("Disorder structure"),
    unit_cell_lengths: f64([10, 10, 10]),
    unit_cell_angles: f64([90, 90, 90]),
    unit_cells: i32([2, 1, 1]),
    number_of_atoms: u64([5]),
    number_of_types: u64([4]),
    types_names: scalarStr("O;H;N;H"),
    atom_type: i32([1, 2, 2, 3, 4]),
    atom_position: f64([...x, ...y, ...z], [3, 5]),
    atom_unit_cell: i32([1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], [3, 5]),
    magnetic_spins: f64(new Array(15).fill(-1), [3, 5]),
  };
}

test("reads RMCProfile 6.8 structure with ';'-joined type names", async () => {
  const ctx = await loadH5(["js/unified_hdf5.js"]);
  const bytes = writeEntryData(ctx.h5wasm, rmcProfile68());
  const model = ctx.UnifiedH5.readStructure(ctx.h5wasm, bytes, "rmc68.h5");
  assert.deepEqual(plain(model.typeNames), ["Na", "Cl"]);
  assert.deepEqual(plain(model.elements), ["Na", "Na", "Cl", "Cl"]);
  assert.deepEqual(plain(model.unitCells), [2, 1, 1]);
  assert.deepEqual(plain(model.cellLengths), [4, 9, 10]);
  assert.deepEqual(plain(model.atomPosition[1]), [1.25, 0.75, 0.875]);
  assert.deepEqual(plain(model.atomUnitCell[1]), [2, 1, 1]);
  assert.equal(model.atomUnitCellPresent, true);
  assert.equal(model.coordinatesAreSupercellFractional, 0);
  assert.equal(model.legacyContract, false);
});

test("accepts space- and comma-separated type names", async () => {
  const ctx = await loadH5(["js/unified_hdf5.js"]);
  for (const names of ["Na Cl", "Na,Cl", "Na; Cl"]) {
    const bytes = writeEntryData(ctx.h5wasm, rmcProfile68({ types_names: str(names) }));
    const model = ctx.UnifiedH5.readStructure(ctx.h5wasm, bytes, "names.h5");
    assert.deepEqual(plain(model.elements), ["Na", "Na", "Cl", "Cl"], names);
  }
});

test("reads DISCUS (3,N) layout with scalar strings and uint64 counts", async () => {
  const ctx = await loadH5(["js/unified_hdf5.js"]);
  const bytes = writeEntryData(ctx.h5wasm, discusLayout(), { status_flag_is_super_structure: "is_true" });
  const model = ctx.UnifiedH5.readStructure(ctx.h5wasm, bytes, "discus.h5");
  assert.equal(model.atomCount, 5);
  assert.deepEqual(plain(model.elements), ["O", "H", "H", "N", "H"]);
  assert.deepEqual(plain(model.atomPosition[2]), [1.13, 0, 0]);
  assert.deepEqual(plain(model.atomUnitCell[3]), [2, 1, 1]);
  assert.equal(model.atomUnitCellPresent, true);
});

test("reports a missing atom_unit_cell instead of hiding it", async () => {
  const ctx = await loadH5(["js/unified_hdf5.js"]);
  const layout = rmcProfile68();
  delete layout.atom_unit_cell;
  const model = ctx.UnifiedH5.readStructure(ctx.h5wasm, writeEntryData(ctx.h5wasm, layout), "nocells.h5");
  assert.equal(model.atomUnitCellPresent, false);
});

test("treats an all-ones atom_unit_cell in a multi-cell supercell as missing", async () => {
  const ctx = await loadH5(["js/unified_hdf5.js"]);
  const ones = i32(new Array(12).fill(1), [4, 3]);
  const defaulted = ctx.UnifiedH5.readStructure(ctx.h5wasm,
    writeEntryData(ctx.h5wasm, rmcProfile68({ atom_unit_cell: ones })), "ones.h5");
  assert.equal(defaulted.atomUnitCellPresent, false);
  assert.equal(defaulted.atomUnitCellDefaulted, true);
  const singleCell = ctx.UnifiedH5.readStructure(ctx.h5wasm,
    writeEntryData(ctx.h5wasm, rmcProfile68({ atom_unit_cell: ones, unit_cells: i32([1, 1, 1]) })), "one_cell.h5");
  assert.equal(singleCell.atomUnitCellPresent, true);
  assert.equal(singleCell.atomUnitCellDefaulted, false);
});

test("rejects atom types outside types_names", async () => {
  const ctx = await loadH5(["js/unified_hdf5.js"]);
  const bytes = writeEntryData(ctx.h5wasm, rmcProfile68({ atom_type: i32([1, 1, 2, 3]) }));
  assert.throws(() => ctx.UnifiedH5.readStructure(ctx.h5wasm, bytes, "bad.h5"), /outside types_names/);
});

test("writes with_bragg and grid layout into unified data", async () => {
  const ctx = await loadH5(["js/unified_hdf5.js"]);
  const h5 = ctx.h5wasm;
  const values = Float64Array.from({ length: 24 }, (_, i) => i);
  const bytes = ctx.UnifiedH5.writeData(h5, {
    dims: [2, 3, 4], corner: [-1, 0, 0.5], vectors: [[0.5, 0, 0], [0, 0.25, 0], [0, 0, 1]],
    values, cell: [4, 5, 6, 90, 90, 90], withBragg: "bragg_subtracted", radiation: "neutron",
  });
  const file = "/tmp/written_data.h5";
  h5.FS.writeFile(file, bytes);
  const f = new h5.File(file, "r");
  try {
    assert.equal(String(f.get("entry/data/data_type_with_bragg").value), "bragg_subtracted");
    const ds = f.get("entry/data/data_values");
    assert.deepEqual(plain(ds.shape), [2, 3, 4]);
    assert.equal(ds.value[1 * 12 + 2 * 4 + 3], 23);
    assert.deepEqual(plain(Array.from(f.get("entry/data/data_corner").value)), [-1, 0, 0.5]);
  } finally {
    f.close();
    h5.FS.unlink(file);
  }
});
