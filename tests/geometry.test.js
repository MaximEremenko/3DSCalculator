"use strict";

const test = require("node:test");
const assert = require("node:assert/strict");
const { loadScripts } = require("./helpers/browser_env");

const ctx = loadScripts(["js/rmc6f_reader.js"]);
const reader = ctx.RMC6fReader;

const CELLS = {
  orthogonal: [4, 3.5, 3, 90, 90, 90],
  monoclinic: [5.1, 4.2, 6.3, 90, 101.5, 90],
  hexagonal: [3.2, 3.2, 5.2, 90, 90, 120],
  triclinic: [4.4, 3.9, 3.1, 84, 95, 103],
};

const norm = (v) => Math.hypot(v[0], v[1], v[2]);
const dotp = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const angle = (a, b) => (Math.acos(dotp(a, b) / (norm(a) * norm(b))) * 180) / Math.PI;
const hklToQ = (h, Bq) => [0, 1, 2].map((j) => h[0] * Bq[0][j] + h[1] * Bq[1][j] + h[2] * Bq[2][j]);

function rmc6f(cell, supercell, atoms) {
  const rows = atoms.map(([element, u, c], i) =>
    `${i + 1} ${element} [1] ${u.map((v, d) => v / supercell[d]).join(" ")} 1 ${c.join(" ")}`);
  return [
    "(Version 6f format configuration file)",
    `Number of atoms: ${atoms.length}`,
    `Supercell dimensions: ${supercell.join(" ")}`,
    `Cell (Ang/deg): ${cell.slice(0, 3).map((v, d) => v * supercell[d]).join(" ")} ${cell.slice(3).join(" ")}`,
    "Atoms:",
    ...rows,
  ].join("\n");
}

for (const [label, cell] of Object.entries(CELLS)) {
  test(`${label} cell: lattice rows reproduce a, b, c and angles`, () => {
    const g = reader.cellGeometry(cell);
    const [a, b, c] = g.direct;
    for (const [got, want] of [[norm(a), cell[0]], [norm(b), cell[1]], [norm(c), cell[2]],
      [angle(b, c), cell[3]], [angle(a, c), cell[4]], [angle(a, b), cell[5]]]) {
      assert.ok(Math.abs(got - want) < 1e-9, `${got} vs ${want}`);
    }
    for (let i = 0; i < 3; i++)
      for (let j = 0; j < 3; j++)
        assert.ok(Math.abs(dotp(g.Bp[i], g.direct[j]) - (i === j ? 1 : 0)) < 1e-12);
  });

  test(`${label} cell: q.r equals 2*pi*h.u for RMC6f atoms`, () => {
    const supercell = [3, 2, 4];
    const atoms = [["Fe", [1.37, 0.21, 2.64], [1, 0, 2]], ["O", [2.9, 1.45, 3.05], [2, 1, 3]]];
    const parsed = reader.parse("t.rmc6f", rmc6f(cell, supercell, atoms));
    for (let i = 0; i < atoms.length; i++) {
      const r = [parsed.x[i], parsed.y[i], parsed.z[i]];
      const ra = [parsed.xa[i], parsed.ya[i], parsed.za[i]];
      for (const h of [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1.5, -2, 0.7]]) {
        const q = hklToQ(h, parsed.Bq);
        const phase = 2 * Math.PI * dotp(h, atoms[i][1]);
        const idealPhase = 2 * Math.PI * dotp(h, atoms[i][2]);
        assert.ok(Math.abs(dotp(q, r) - phase) < 1e-9, `${label} atom ${i} h=${h}`);
        assert.ok(Math.abs(dotp(q, ra) - idealPhase) < 1e-9, `${label} ideal atom ${i} h=${h}`);
      }
    }
    assert.deepEqual(Array.from(parsed.cellIndex), [1, 0, 2, 2, 1, 3]);
  });
}

test("structure without cell references keeps wrapped positions and no cell indices", () => {
  const s = reader.buildSupercellStructure({
    parentCellDeg: CELLS.hexagonal, supercell: [2, 2, 1],
    fractional: Float64Array.from([0.9, 0.1, 0.5]), cellIndex: null, elements: ["Fe"],
  });
  assert.equal(s.cellIndex, null);
  assert.deepEqual([s.xa[0], s.ya[0], s.za[0]], [0, 0, 0]);
  const q = hklToQ([1, 1, 1], s.Bq);
  // u = (0.9 - 1, 0.1, 0.5) * supercell after wrapping to the nearest image
  const u = [-0.2, 0.2, 0.5];
  assert.ok(Math.abs(dotp(q, [s.x[0], s.y[0], s.z[0]]) - 2 * Math.PI * dotp([1, 1, 1], u)) < 1e-9);
});
