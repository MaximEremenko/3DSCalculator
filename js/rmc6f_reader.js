"use strict";

(function (global) {
  const PI2 = 2 * Math.PI;

  class RMC6fParser {
    constructor() {
      this.header = [];
      this.rows = [];
      this.super = [1, 1, 1];
      this.cellDeg = [NaN, NaN, NaN, NaN, NaN, NaN];
      this.cellRad = [NaN, NaN, NaN, NaN, NaN, NaN];
      this.cols = 0;
      this.names = [];
    }

    parse(txt) {
      const lines = String(txt || "")
        .split(/\r?\n/)
        .map((line) => line.replace(/\x00/g, ""));
      const start = this.findHeader(lines);
      this.header = lines.slice(0, start);
      this.rows = this.parseRows(lines.slice(start));
      this.getSuper();
      this.getCell();
    }

    findHeader(lines) {
      for (let i = 0; i < Math.min(lines.length, 151); i++) {
        if (lines[i].includes("Atoms:")) return i + 1;
      }
      return 0;
    }

    parseRows(lines) {
      const out = [];
      for (const raw of lines) {
        const line = raw.trim();
        if (!line) continue;
        const p = line.split(/\s+/);
        if (!this.cols) {
          this.cols = p.length;
          if (this.cols === 10) {
            this.names = [
              "atomNumber",
              "element",
              "id",
              "x",
              "y",
              "z",
              "refNumber",
              "cellRefNumX",
              "cellRefNumY",
              "cellRefNumZ",
            ];
          } else if (this.cols === 9) {
            this.names = [
              "atomNumber",
              "element",
              "x",
              "y",
              "z",
              "refNumber",
              "cellRefNumX",
              "cellRefNumY",
              "cellRefNumZ",
            ];
          } else {
            throw new Error(`Unsupported RMC6f: ${this.cols} columns`);
          }
        }
        const row = {};
        for (let i = 0; i < this.cols; i++) {
          const name = this.names[i] || `c${i}`;
          row[name] = [
            "atomNumber",
            "x",
            "y",
            "z",
            "refNumber",
            "cellRefNumX",
            "cellRefNumY",
            "cellRefNumZ",
          ].includes(name)
            ? Number(p[i])
            : p[i];
        }
        out.push(row);
      }
      return out;
    }

    getSuper() {
      const line = this.header.find((v) => v.includes("Supercell"));
      if (!line) {
        this.super = [1, 1, 1];
        return;
      }
      const nums = line
        .split(/\s+/)
        .map((v) => parseInt(v, 10))
        .filter(Number.isInteger);
      this.super = nums.length >= 3 ? nums.slice(0, 3) : [1, 1, 1];
    }

    getCell() {
      const line =
        this.header.find((v) => v.includes("Cell (Ang/deg)")) ||
        this.header.find((v) => /^\s*Cell\b/i.test(v));
      if (!line) return;
      const values = line
        .split(/\s+/)
        .map(Number)
        .filter(Number.isFinite);
      if (values.length >= 6) {
        this.cellDeg = values.slice(0, 6);
        this.cellRad = [
          values[0],
          values[1],
          values[2],
          (values[3] * Math.PI) / 180,
          (values[4] * Math.PI) / 180,
          (values[5] * Math.PI) / 180,
        ];
      }
    }
  }

  function inv3(m) {
    const a = m[0][0];
    const b = m[0][1];
    const c = m[0][2];
    const d = m[1][0];
    const e = m[1][1];
    const f = m[1][2];
    const g = m[2][0];
    const h = m[2][1];
    const i = m[2][2];
    const A = e * i - f * h;
    const B = -(d * i - f * g);
    const C = d * h - e * g;
    const D = -(b * i - c * h);
    const E = a * i - c * g;
    const F = -(a * h - b * g);
    const G = b * f - c * e;
    const H = -(a * f - c * d);
    const I = a * e - b * d;
    const det = a * A + b * B + c * C;
    if (!Number.isFinite(det) || Math.abs(det) < 1e-20) {
      throw new Error("Singular matrix");
    }
    const s = 1 / det;
    return [
      [A * s, D * s, G * s],
      [B * s, E * s, H * s],
      [C * s, F * s, I * s],
    ];
  }

  function mul(v, m) {
    return [
      v[0] * m[0][0] + v[1] * m[1][0] + v[2] * m[2][0],
      v[0] * m[0][1] + v[1] * m[1][1] + v[2] * m[2][1],
      v[0] * m[0][2] + v[1] * m[1][2] + v[2] * m[2][2],
    ];
  }

  function transpose(m) {
    return [
      [m[0][0], m[1][0], m[2][0]],
      [m[0][1], m[1][1], m[2][1]],
      [m[0][2], m[1][2], m[2][2]],
    ];
  }

  function scaleMat(m, s) {
    return [
      [m[0][0] * s, m[0][1] * s, m[0][2] * s],
      [m[1][0] * s, m[1][1] * s, m[1][2] * s],
      [m[2][0] * s, m[2][1] * s, m[2][2] * s],
    ];
  }

  // Lattice vectors as the ROWS of the returned matrix (c along z, b in the
  // yz plane); angles in radians.
  function cell2vec(a, b, c, al, be, ga) {
    const v = [
      [0, 0, 0],
      [0, 0, 0],
      [0, 0, 0],
    ];
    v[2][0] = a * Math.cos(be);
    v[1][1] = b * Math.sin(al);
    v[2][1] = b * Math.cos(al);
    v[2][2] = c;
    v[1][0] = (a * b * Math.cos(ga) - v[2][0] * v[2][1]) / v[1][1];
    v[0][0] = Math.sqrt(
      Math.max(0, a * a - v[1][0] * v[1][0] - v[2][0] * v[2][0])
    );
    return [
      [v[0][0], v[1][0], v[2][0]],
      [v[0][1], v[1][1], v[2][1]],
      [v[0][2], v[1][2], v[2][2]],
    ];
  }

  function wrapDelta(x) {
    if (x < -0.5) return x + 1;
    if (x > 0.5) return x - 1;
    return x;
  }

  // Cartesian geometry of the parent cell (lengths in Angstrom, angles in
  // degrees). Rows of `direct` are the lattice vectors and rows of `Bp` the
  // reciprocal vectors (a*.a = 1), so r = u * direct and q = [h,k,l] * Bq
  // give q.r = 2*pi*h.u for u in parent-cell units, for any cell angles.
  function cellGeometry(cellDeg) {
    const rad = (deg) => (Number(deg) * Math.PI) / 180;
    const direct = cell2vec(
      Number(cellDeg[0]),
      Number(cellDeg[1]),
      Number(cellDeg[2]),
      rad(cellDeg[3]),
      rad(cellDeg[4]),
      rad(cellDeg[5])
    );
    const Bp = transpose(inv3(direct));
    return { direct, Bp, Bq: scaleMat(Bp, PI2) };
  }

  // Builds the parsed structure from supercell-fractional positions (3 per
  // atom) and the 0-based cell each atom belongs to (null when the source
  // has no cell references). Ideal positions are the cell origins and the
  // displacements are wrapped to the nearest periodic image.
  function buildSupercellStructure(spec) {
    const supercell = spec.supercell.map((v) => Math.max(1, Math.round(Number(v) || 1)));
    const geometry = cellGeometry(spec.parentCellDeg);
    const frac = spec.fractional;
    const cells = spec.cellIndex || null;
    const n = spec.elements.length;
    const x = new Float64Array(n);
    const y = new Float64Array(n);
    const z = new Float64Array(n);
    const xa = new Float64Array(n);
    const ya = new Float64Array(n);
    const za = new Float64Array(n);
    const dx = new Float64Array(n);
    const dy = new Float64Array(n);
    const dz = new Float64Array(n);
    const ideal = [0, 0, 0];
    const delta = [0, 0, 0];

    for (let i = 0; i < n; i++) {
      for (let d = 0; d < 3; d++) {
        const cell = cells ? cells[i * 3 + d] : 0;
        ideal[d] = cell;
        delta[d] = wrapDelta(Number(frac[i * 3 + d]) - cell / supercell[d]) * supercell[d];
      }
      const pA = mul(ideal, geometry.direct);
      const pD = mul(delta, geometry.direct);
      xa[i] = pA[0];
      ya[i] = pA[1];
      za[i] = pA[2];
      dx[i] = pD[0];
      dy[i] = pD[1];
      dz[i] = pD[2];
      x[i] = pA[0] + pD[0];
      y[i] = pA[1] + pD[1];
      z[i] = pA[2] + pD[2];
    }

    return {
      atoms: n,
      super: supercell,
      Bp: geometry.Bp,
      Bq: geometry.Bq,
      x,
      y,
      z,
      xa,
      ya,
      za,
      dx,
      dy,
      dz,
      cellIndex: cells,
      elements: Array.from(spec.elements, (e) => String(e || "").trim()),
    };
  }

  function parseRmc6f(name, text) {
    const r = new RMC6fParser();
    r.parse(text);
    if (!r.rows.length) throw new Error("No atom rows found");
    if (!r.cellRad.every(Number.isFinite)) throw new Error("Cell parameters missing");

    const [sx, sy, sz] = r.super;
    const [a, b, c, al, be, ga] = r.cellDeg;
    const n = r.rows.length;
    const fractional = new Float64Array(n * 3);
    const cellIndex = new Int32Array(n * 3);
    const elements = new Array(n);
    for (let i = 0; i < n; i++) {
      const row = r.rows[i];
      fractional[i * 3] = Number(row.x);
      fractional[i * 3 + 1] = Number(row.y);
      fractional[i * 3 + 2] = Number(row.z);
      cellIndex[i * 3] = Math.round(Number(row.cellRefNumX));
      cellIndex[i * 3 + 1] = Math.round(Number(row.cellRefNumY));
      cellIndex[i * 3 + 2] = Math.round(Number(row.cellRefNumZ));
      elements[i] = row.element;
    }

    return {
      file: String(name || ""),
      ...buildSupercellStructure({
        parentCellDeg: [a / sx, b / sy, c / sz, al, be, ga],
        supercell: [sx, sy, sz],
        fractional,
        cellIndex,
        elements,
      }),
      cellDeg: r.cellDeg.slice(0, 6),
    };
  }

  global.RMC6fReader = Object.freeze({
    parse: parseRmc6f,
    Parser: RMC6fParser,
    cellGeometry,
    buildSupercellStructure,
  });
})(typeof window !== "undefined" ? window : globalThis);
