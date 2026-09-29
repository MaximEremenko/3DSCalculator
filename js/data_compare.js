"use strict";

// Compares a calculated intensity grid with measured or reference data on
// the h, k, l points they share. The model
//   I_model = s I_calc + b + c |Q|
// is fitted by weighted least squares with w = 1 / sigma^2 (equal weights
// when the data give no sigma), each of s, b and c refined or held at
// s = 1, b = c = 0: Spinteract's single-crystal fit. A data point is left
// out when it is masked (NaN), its sigma is not positive, the calculation
// there is not finite, or Q = 0. With r = I_data - I_model:
//   chi2 = sum w r^2        R_wp = 100 sqrt(chi2 / sum w I_data^2)  (as Spinteract)
//   R = 100 sum |r| / sum |I_data|        chi2_nu = chi2 / (n - parameters)
// Grids are h, k, l axes with l fastest, as the calculator and
// IntensityListReader.buildGrid give them; Q = [h, k, l] * Bq.
(function (global) {
  // Index of each value of `from` in the regular axis `to`, or -1.
  function axisIndex(from, to) {
    const out = new Int32Array(from.length).fill(-1);
    const n = to.length;
    if (!n) return out;
    const step = n > 1 ? (to[n - 1] - to[0]) / (n - 1) : 1;
    const tol = n > 1 ? Math.min(2e-6, Math.abs(step) / 4) : 2e-6;
    for (let i = 0; i < from.length; i++) {
      const j = n > 1 ? Math.round((from[i] - to[0]) / step) : 0;
      if (j >= 0 && j < n && Math.abs(to[j] - from[i]) <= tol) out[i] = j;
    }
    return out;
  }

  // Where each data grid point sits in the calculated grid; `shared` counts
  // the data grid points the calculation has, `same` means identical grids.
  function matchGrids(data, calc) {
    const ih = axisIndex(data.h, calc.h);
    const ik = axisIndex(data.k, calc.k);
    const il = axisIndex(data.l, calc.l);
    const found = (a) => a.reduce((n, v) => n + (v >= 0 ? 1 : 0), 0);
    const shared = found(ih) * found(ik) * found(il);
    return { ih, ik, il, shared, same: shared === data.I.length && calc.I.length === data.I.length };
  }

  // |Q| from h, k, l with Q = [h, k, l] * Bq (rows of Bq are the reciprocal vectors).
  function qLength(Bq) {
    const b = Bq || [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    return function (h, k, l) {
      const x = h * b[0][0] + k * b[1][0] + l * b[2][0];
      const y = h * b[0][1] + k * b[1][1] + l * b[2][1];
      const z = h * b[0][2] + k * b[1][2] + l * b[2][2];
      return Math.sqrt(x * x + y * y + z * z);
    };
  }

  function noneSkipped() {
    return { masked: 0, outside: 0, sigma: 0, calc: 0, origin: 0 };
  }

  // Calls fn(at, y, w, x, q) for each data point that enters the fit: `at`
  // indexes the data grid, y is the data, w its weight, x the calculation
  // and q = |Q|. The others are counted in `skipped`.
  function eachPoint(data, calc, match, Bq, fn, skipped) {
    const [nh, nk, nl] = data.shape;
    const ck = calc.shape[1], cl = calc.shape[2];
    const sigma = data.sigma || null;
    const w0 = data.sigmaConstant > 0 ? 1 / (data.sigmaConstant * data.sigmaConstant) : 1;
    const qOf = qLength(Bq);
    for (let i = 0; i < nh; i++) {
      const ci = match.ih[i];
      for (let j = 0; j < nk; j++) {
        const cj = match.ik[j];
        for (let m = 0; m < nl; m++) {
          const at = (i * nk + j) * nl + m;
          const y = data.I[at];
          if (!Number.isFinite(y)) { skipped.masked++; continue; }
          const cm = match.il[m];
          if (ci < 0 || cj < 0 || cm < 0) { skipped.outside++; continue; }
          let w = w0;
          if (sigma) {
            const s = sigma[at];
            if (!(s > 0 && s < Infinity)) { skipped.sigma++; continue; }
            w = 1 / (s * s);
          }
          const x = calc.I[(ci * ck + cj) * cl + cm];
          if (!Number.isFinite(x)) { skipped.calc++; continue; }
          const q = qOf(data.h[i], data.k[j], data.l[m]);
          if (q < 1e-9) { skipped.origin++; continue; }
          fn(at, y, w, x, q);
        }
      }
    }
  }

  // Solves the symmetric system M p = v (n <= 3) by Gaussian elimination with
  // partial pivoting on the equilibrated matrix; null when it is singular.
  function solve(M, v) {
    const n = v.length;
    const d = M.map((row, i) => (row[i] > 0 ? 1 / Math.sqrt(row[i]) : 1));
    const a = M.map((row, i) => row.map((x, j) => x * d[i] * d[j]).concat([v[i] * d[i]]));
    for (let c = 0; c < n; c++) {
      let piv = c;
      for (let r = c + 1; r < n; r++) if (Math.abs(a[r][c]) > Math.abs(a[piv][c])) piv = r;
      if (!(Math.abs(a[piv][c]) > 1e-12)) return null;
      [a[c], a[piv]] = [a[piv], a[c]];
      for (let r = c + 1; r < n; r++) {
        const f = a[r][c] / a[c][c];
        for (let k = c; k <= n; k++) a[r][k] -= f * a[c][k];
      }
    }
    const p = new Array(n).fill(0);
    for (let r = n - 1; r >= 0; r--) {
      let s = a[r][n];
      for (let k = r + 1; k < n; k++) s -= a[r][k] * p[k];
      p[r] = s / a[r][r];
    }
    return p.map((x, i) => x * d[i]);
  }

  // Fits the model to the data. options: { scale: refine s (default true),
  // background: "none" | "flat" (default) | "linear" (b + c |Q|), Bq }.
  // Returns { scale, background, linear, n, parameters, chi2, chi2nu, rwp, r,
  // weights: "sigma" | "constant" | "equal", skipped } or { error, ... }.
  function fit(data, calc, match, options) {
    const o = options || {};
    const bg = o.background || "flat";
    const free = [o.scale !== false, bg === "flat" || bg === "linear", bg === "linear"];
    const skipped = noneSkipped();
    // Moments of the columns [x, 1, q] and of z = y (or y - x with s held at 1)
    let n = 0, Sxx = 0, Sx = 0, Sxq = 0, S1 = 0, Sq = 0, Sqq = 0, Sxz = 0, Sz = 0, Sqz = 0;
    eachPoint(data, calc, match, o.Bq, function (at, y, w, x, q) {
      const z = free[0] ? y : y - x;
      n++;
      Sxx += w * x * x; Sx += w * x; Sxq += w * x * q;
      S1 += w; Sq += w * q; Sqq += w * q * q;
      Sxz += w * x * z; Sz += w * z; Sqz += w * q * z;
    }, skipped);
    const weights = data.sigma ? "sigma" : data.sigmaConstant > 0 ? "constant" : "equal";
    const cols = [0, 1, 2].filter((c) => free[c]);
    const parameters = cols.length;
    const base = { n, parameters, weights, skipped };
    if (n <= parameters) return { ...base, error: n ? `only ${n} data point(s) for ${parameters} parameter(s)` : "no data point is shared with the calculation" };
    const M = [[Sxx, Sx, Sxq], [Sx, S1, Sq], [Sxq, Sq, Sqq]];
    const v = [Sxz, Sz, Sqz];
    const p = [free[0] ? 0 : 1, 0, 0];
    if (parameters) {
      const sol = solve(cols.map((a) => cols.map((b) => M[a][b])), cols.map((a) => v[a]));
      if (!sol) return { ...base, error: "the calculation does not vary on the data points, so the scale and the background cannot both be fitted" };
      cols.forEach((c, i) => { p[c] = sol[i]; });
    }
    const [scale, background, linear] = p;
    let chi2 = 0, absR = 0, absY = 0, wyy = 0;
    eachPoint(data, calc, match, o.Bq, function (at, y, w, x, q) {
      const r = y - (scale * x + background + linear * q);
      chi2 += w * r * r;
      absR += Math.abs(r);
      absY += Math.abs(y);
      wyy += w * y * y;
    }, noneSkipped());
    return {
      ...base,
      scale, background, linear, chi2,
      chi2nu: chi2 / (n - parameters),
      rwp: wyy > 0 ? 100 * Math.sqrt(chi2 / wyy) : NaN,
      r: absY > 0 ? (100 * absR) / absY : NaN,
    };
  }

  // s I_calc + b + c |Q| on the calculated grid.
  function modelGrid(calc, result, Bq) {
    const [nh, nk, nl] = calc.shape;
    const out = new Float64Array(calc.I.length);
    const s = result.scale, b = result.background, c = result.linear;
    const qOf = qLength(Bq);
    for (let i = 0; i < nh; i++) {
      for (let j = 0; j < nk; j++) {
        for (let m = 0; m < nl; m++) {
          const at = (i * nk + j) * nl + m;
          out[at] = s * calc.I[at] + b + (c ? c * qOf(calc.h[i], calc.k[j], calc.l[m]) : 0);
        }
      }
    }
    return out;
  }

  // I_data - I_model on the data grid, at the points of the fit; NaN elsewhere.
  function differenceGrid(data, calc, match, result, Bq) {
    const out = new Float64Array(data.I.length).fill(NaN);
    const s = result.scale, b = result.background, c = result.linear;
    eachPoint(data, calc, match, Bq, function (at, y, w, x, q) {
      out[at] = y - (s * x + b + c * q);
    }, noneSkipped());
    return out;
  }

  // Largest |v| and its `quantile` (default 0.995, from at most ~200k
  // samples), for color levels symmetric about zero.
  function symmetricLevels(values, quantile) {
    const qn = Number.isFinite(quantile) ? quantile : 0.995;
    const stride = Math.max(1, Math.floor(values.length / 200000));
    const sample = [];
    let max = 0;
    for (let i = 0; i < values.length; i++) {
      const a = Math.abs(values[i]);
      if (!Number.isFinite(a)) continue;
      if (a > max) max = a;
      if (i % stride === 0) sample.push(a);
    }
    if (!sample.length) return { max: 0, level: 0 };
    sample.sort((a, b) => a - b);
    const level = sample[Math.min(sample.length - 1, Math.floor(qn * (sample.length - 1)))];
    return { max, level: level > 0 ? level : max };
  }

  global.DataCompare = Object.freeze({ axisIndex, matchGrids, fit, modelGrid, differenceGrid, symmetricLevels });
})(typeof window !== "undefined" ? window : globalThis);
