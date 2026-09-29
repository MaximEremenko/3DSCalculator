"use strict";

// Laue symmetrization of intensities on an h, k, l grid. The 12 Laue classes
// and their operations on (h, k, l) are those of Scatty (spinsym_int): the
// trigonal and hexagonal classes use hexagonal axes, with i = -h - k.
// The symmetrized intensity at a grid point is the mean over its equivalents
// g(h, k, l) that fall on the grid; masked (NaN) values are left out.
(function (global) {
  const BASIS = { h: [1, 0, 0], k: [0, 1, 0], i: [-1, -1, 0], l: [0, 0, 1] };

  // "k,i,-l" -> rows of the new h, k, l in terms of the old ones
  function parseOp(text) {
    return text.split(",").map((t) => {
      const s = t.trim(), neg = s[0] === "-", b = BASIS[neg ? s.slice(1) : s];
      if (!b) throw new Error(`Laue symmetry: bad operation "${text}"`);
      return b.map((v) => (neg ? -v : v));
    });
  }
  const invert = (m) => m.map((row) => row.map((v) => -v));
  const withInversion = (ops) => ops.concat(ops.map(invert));
  const fromText = (list) => withInversion(list.map(parseOp));

  // signed permutations of (h, k, l) for the given permutations, all sign combinations
  function signedPermutations(perms) {
    const out = [];
    for (const p of perms) {
      for (let s = 0; s < 8; s++) out.push([0, 1, 2].map((r) => { const row = [0, 0, 0]; row[p[r]] = (s >> r) & 1 ? -1 : 1; return row; }));
    }
    return out;
  }

  const ROT6 = ["h,k,l", "k,i,l", "i,h,l", "-h,-k,l", "-k,-i,l", "-i,-h,l"];
  const CLASSES = {
    "-1": () => fromText(["h,k,l"]),
    "2/m": () => fromText(["h,k,l", "-h,k,-l"]),                               // b unique
    "mmm": () => signedPermutations([[0, 1, 2]]),
    "4/m": () => fromText(["h,k,l", "-k,h,l", "-h,-k,l", "k,-h,l"]),
    "4/mmm": () => signedPermutations([[0, 1, 2], [1, 0, 2]]),
    "-3": () => fromText(["h,k,l", "k,i,l", "i,h,l"]),
    "-31m": () => fromText(["h,k,l", "k,i,l", "i,h,l", "-k,-h,-l", "-h,-i,-l", "-i,-k,-l"]),
    "-3m1": () => fromText(["h,k,l", "k,i,l", "i,h,l", "-k,-h,l", "-h,-i,l", "-i,-k,l"]),
    "6/m": () => fromText(ROT6),
    "6/mmm": () => fromText(ROT6.concat(["k,h,-l", "h,i,-l", "i,k,-l", "-k,-h,-l", "-h,-i,-l", "-i,-k,-l"])),
    "m-3": () => signedPermutations([[0, 1, 2], [1, 2, 0], [2, 0, 1]]),
    "m-3m": () => signedPermutations([[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]]),
  };

  function operations(laue) {
    const make = CLASSES[laue];
    if (!make) throw new Error(`Laue symmetry: unknown class "${laue}"`);
    return make();
  }

  // Largest relative change of the reciprocal metric G* = Bq Bq^T under the
  // operations: 0 when the cell has the class's symmetry (|Q| is kept).
  function metricDeviation(ops, Bq) {
    const G = [0, 1, 2].map((a) => [0, 1, 2].map((b) => Bq[a][0] * Bq[b][0] + Bq[a][1] * Bq[b][1] + Bq[a][2] * Bq[b][2]));
    const scale = Math.max(G[0][0], G[1][1], G[2][2]);
    let worst = 0;
    for (const M of ops) {
      for (let a = 0; a < 3; a++) {
        for (let b = 0; b < 3; b++) {
          let v = 0;
          for (let c = 0; c < 3; c++) for (let d = 0; d < 3; d++) v += M[c][a] * G[c][d] * M[d][b];
          worst = Math.max(worst, Math.abs(v - G[a][b]) / scale);
        }
      }
    }
    return worst;
  }

  // grid: { h, k, l, shape } with l fastest. Returns { I, min, max, meta }.
  function symmetrize(I, grid, laue, options) {
    const ops = operations(laue);
    const axes = [grid.h, grid.k, grid.l], shape = grid.shape;
    const lo = axes.map((a) => a[0]);
    const step = axes.map((a) => (a.length > 1 ? (a[a.length - 1] - a[0]) / (a.length - 1) : 1));
    const [nh, nk, nl] = shape;
    const out = new Float64Array(I.length);
    const index = (v, d) => {
      const a = Math.round((v - lo[d]) / step[d]);
      return a >= 0 && a < shape[d] && Math.abs(lo[d] + a * step[d] - v) <= 1e-3 * Math.abs(step[d]) + 1e-12 ? a : -1;
    };
    let used = 0, min = Infinity, max = -Infinity;
    for (let i = 0; i < nh; i++) {
      const h = axes[0][i];
      for (let j = 0; j < nk; j++) {
        const k = axes[1][j];
        for (let m = 0; m < nl; m++) {
          const l = axes[2][m];
          let sum = 0, count = 0;
          for (const M of ops) {
            const a = index(M[0][0] * h + M[0][1] * k + M[0][2] * l, 0);
            if (a < 0) continue;
            const b = index(M[1][0] * h + M[1][1] * k + M[1][2] * l, 1);
            if (b < 0) continue;
            const c = index(M[2][0] * h + M[2][1] * k + M[2][2] * l, 2);
            if (c < 0) continue;
            const v = I[(a * nk + b) * nl + c];
            if (!Number.isFinite(v)) continue;
            sum += v;
            count++;
          }
          const v = count ? sum / count : NaN;
          out[(i * nk + j) * nl + m] = v;
          used += count;
          if (v < min) min = v;
          if (v > max) max = v;
        }
      }
    }
    const Bq = options && options.Bq;
    const deviation = Bq ? metricDeviation(ops, Bq) : 0;
    return { I: out, min, max, meta: { laue, operations: ops.length, coverage: used / (I.length * ops.length), metricDeviation: deviation } };
  }

  global.LaueSymmetry = Object.freeze({ CLASSES: Object.keys(CLASSES), operations, metricDeviation, symmetrize });
})(typeof window !== "undefined" ? window : globalThis);
