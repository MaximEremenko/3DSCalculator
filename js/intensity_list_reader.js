"use strict";

// Reads reciprocal-space intensity lists, one point per line:
//   h k l I [sigma]
// as written by Scatty ([title]_[name]_sc_list.txt, sigma = 1) and used for
// Spinteract's single-crystal data (h k l intensity error). A fifth number is
// the uncertainty sigma when every line has one; lines with more numbers are
// read as h k l I. The points are placed on the regular h, k, l grid they
// span; grid points without a line (masked data) are NaN.
// Also reads the cell of a Spinteract configuration (CELL a b c alpha beta
// gamma) and recognises Scatty's settings files, which hold no structure.
(function (global) {
  const QUANT = 1e6; // hkl values are matched on a 1e-6 grid
  const MAX_ABS_HKL = 2000;
  const MAX_GRID_POINTS = 60e6;

  const POW10 = Array.from({ length: 23 }, (_, i) => Number(`1e${i}`));
  function pow10(e) {
    return e < POW10.length ? POW10[e] : Math.pow(10, e);
  }

  function isNumber(token) {
    return token !== "" && Number.isFinite(Number(token));
  }

  // True when the first data lines of `head` hold 4 to 8 numbers each.
  function detect(head) {
    let rows = 0;
    for (const raw of String(head || "").split(/\r?\n/, 80)) {
      const t = raw.trim();
      if (!t || t[0] === "#") continue;
      const tokens = t.split(/\s+/);
      if (tokens.length < 4 || tokens.length > 8 || !tokens.every(isNumber)) return false;
      if (++rows >= 6) break;
    }
    return rows >= 3;
  }

  // "spinteract" (a configuration with CELL), "scatty-settings"
  // (scatty_config.txt), "scatty-info" (..._scatty_info.txt), or null.
  function detectConfig(text) {
    const head = String(text || "").slice(0, 20000);
    const has = (key) => new RegExp(`^\\s*${key}\\b`, "im").test(head);
    if (has("CELL") && !has("SPIN") && !has("ATOM") && !has("BOX") && !/Atoms\s*:/.test(head)) return "spinteract";
    if (has("X_AXIS") && has("NAME") && !has("CELL")) return "scatty-settings";
    if (/^\s*Name\s*:/im.test(head) && /Window type\s*:/i.test(head)) return "scatty-info";
    return null;
  }

  // { title, cell: [a, b, c, alpha, beta, gamma] } or null.
  function parseConfig(text) {
    let cell = null;
    let title = "";
    for (const raw of String(text || "").split(/\r?\n/)) {
      const t = raw.trim().split(/\s+/);
      const key = (t[0] || "").toUpperCase();
      if (key === "TITLE") title = raw.trim().slice(5).trim();
      else if (key === "CELL") cell = t.slice(1, 7).map(Number);
    }
    if (!cell || cell.length < 6 || !cell.every(Number.isFinite) || !cell.every((v) => v > 0)) return null;
    return { title, cell };
  }

  // Incremental parser: push text chunks in order, then finish().
  function createParser() {
    let capacity = 1 << 16;
    let count = 0;
    let H = new Int32Array(capacity);
    let K = new Int32Array(capacity);
    let L = new Int32Array(capacity);
    let V = new Float64Array(capacity);
    let E = new Float32Array(capacity);
    let carry = "";
    let skipped = 0;
    let sigmaRows = 0;
    const cols = new Float64Array(6);

    function grow() {
      capacity *= 2;
      const h = new Int32Array(capacity); h.set(H); H = h;
      const k = new Int32Array(capacity); k.set(K); K = k;
      const l = new Int32Array(capacity); l.set(L); L = l;
      const v = new Float64Array(capacity); v.set(V); V = v;
      const e = new Float32Array(capacity); e.set(E); E = e;
    }

    // Decimal number in s[i, j): sign, digits, point, exponent (E or Fortran D).
    function number(s, i, j) {
      let k = i;
      let c = s.charCodeAt(k);
      const negative = c === 45;
      if (c === 45 || c === 43) k++;
      let mantissa = 0, digits = 0, scale = 0, dot = false, any = false;
      for (; k < j; k++) {
        c = s.charCodeAt(k);
        if (c >= 48 && c <= 57) {
          any = true;
          if (digits < 17) {
            mantissa = mantissa * 10 + (c - 48);
            if (mantissa) digits++;
            if (dot) scale--;
          } else if (!dot) {
            scale++;
          }
        } else if (c === 46 && !dot) {
          dot = true;
        } else {
          break;
        }
      }
      if (!any) return NaN;
      if (k < j) {
        if (c !== 69 && c !== 101 && c !== 68 && c !== 100) return NaN;
        k++;
        c = s.charCodeAt(k);
        const down = c === 45;
        if (c === 45 || c === 43) k++;
        let e = 0, exponent = false;
        for (; k < j; k++) {
          c = s.charCodeAt(k);
          if (c < 48 || c > 57) return NaN;
          e = e * 10 + (c - 48);
          exponent = true;
        }
        if (!exponent) return NaN;
        scale += down ? -e : e;
      }
      const v = scale === 0 ? mantissa : scale > 0 ? mantissa * pow10(scale) : mantissa / pow10(-scale);
      return negative ? -v : v;
    }

    // h k l I from the first four numbers of a line, and sigma from a fifth
    // when the line has exactly five.
    function line(s, a, b) {
      let n = 0;
      let i = a;
      while (i < b && n < 6) {
        while (i < b && s.charCodeAt(i) <= 32) i++;
        if (i >= b) break;
        let j = i;
        while (j < b && s.charCodeAt(j) > 32) j++;
        cols[n++] = number(s, i, j);
        i = j;
      }
      if (n === 0) return;
      if (n < 4 || !(Number.isFinite(cols[0]) && Number.isFinite(cols[1]) && Number.isFinite(cols[2]))) {
        skipped++;
        return;
      }
      for (let d = 0; d < 3; d++) {
        if (Math.abs(cols[d]) > MAX_ABS_HKL) throw new Error(`Intensity list: |hkl| above ${MAX_ABS_HKL} on data line ${count + skipped + 1}`);
      }
      if (count === capacity) grow();
      H[count] = Math.round(cols[0] * QUANT);
      K[count] = Math.round(cols[1] * QUANT);
      L[count] = Math.round(cols[2] * QUANT);
      V[count] = cols[3];
      E[count] = n === 5 ? cols[4] : NaN;
      if (n === 5) sigmaRows++;
      count++;
    }

    function push(chunk) {
      const s = carry + chunk;
      let start = 0;
      for (;;) {
        const nl = s.indexOf("\n", start);
        if (nl < 0) break;
        line(s, start, nl);
        start = nl + 1;
      }
      carry = s.slice(start);
    }

    function finish() {
      if (carry) line(carry, 0, carry.length);
      carry = "";
      return { count, skipped, h: H.subarray(0, count), k: K.subarray(0, count), l: L.subarray(0, count), value: V.subarray(0, count), sigma: E.subarray(0, count), sigmaRows };
    }

    return { push, finish, count: () => count };
  }

  function round9(v) {
    return Math.round(v * 1e9) / 1e9;
  }

  // Regular axis spanned by quantized values: every value on min + i * step.
  function regularAxis(q, name) {
    const unique = Array.from(new Set(q)).sort((a, b) => a - b);
    if (unique.length === 1) return { values: [round9(unique[0] / QUANT)], min: unique[0], step: 1 };
    let step = Infinity;
    for (let i = 1; i < unique.length; i++) step = Math.min(step, unique[i] - unique[i - 1]);
    const n = Math.round((unique[unique.length - 1] - unique[0]) / step) + 1;
    const exact = (unique[unique.length - 1] - unique[0]) / (n - 1);
    for (const v of unique) {
      const r = (v - unique[0]) / exact;
      if (Math.abs(r - Math.round(r)) > 1e-3) {
        throw new Error(`Intensity list: the ${name} values do not lie on a regular grid (e.g. ${v / QUANT}); only axis-aligned h, k, l grids are supported`);
      }
    }
    return { values: Array.from({ length: n }, (_, i) => round9((unique[0] + i * exact) / QUANT)), min: unique[0], step: exact };
  }

  // Places parsed points on their h, k, l grid (l fastest). Duplicate points
  // are averaged, with sigma sqrt(sum sigma^2) / n. When every line gives a
  // sigma, `sigma` is its grid (NaN where masked or where sigma <= 0), or
  // null with `sigmaConstant` set when all lines give the same positive one.
  function buildGrid(rows) {
    if (!rows || rows.count < 1) throw new Error("Intensity list: no data lines found");
    const ah = regularAxis(rows.h, "h");
    const ak = regularAxis(rows.k, "k");
    const al = regularAxis(rows.l, "l");
    const shape = [ah.values.length, ak.values.length, al.values.length];
    const total = shape[0] * shape[1] * shape[2];
    if (total > MAX_GRID_POINTS) throw new Error(`Intensity list: the grid ${shape.join(" x ")} has more than ${MAX_GRID_POINTS / 1e6} million points`);
    const I = new Float64Array(total).fill(NaN);
    const hits = new Uint16Array(total);
    const hasSigma = rows.sigma && rows.sigmaRows > 0 && rows.sigmaRows === rows.count;
    const S2 = hasSigma ? new Float64Array(total) : null;
    let duplicates = 0, sigmaInvalid = 0, sigmaMin = Infinity, sigmaMax = -Infinity;
    for (let p = 0; p < rows.count; p++) {
      const i = Math.round((rows.h[p] - ah.min) / ah.step);
      const j = Math.round((rows.k[p] - ak.min) / ak.step);
      const m = Math.round((rows.l[p] - al.min) / al.step);
      const at = (i * shape[1] + j) * shape[2] + m;
      const v = rows.value[p];
      if (hits[at] === 0) I[at] = v;
      else { duplicates++; I[at] += (v - I[at]) / (hits[at] + 1); }
      if (hits[at] < 65535) hits[at]++;
      if (S2) {
        const e = rows.sigma[p];
        if (e > 0 && e < Infinity) {
          S2[at] += e * e;
          if (e < sigmaMin) sigmaMin = e;
          if (e > sigmaMax) sigmaMax = e;
        } else {
          S2[at] = NaN;
          sigmaInvalid++;
        }
      }
    }
    let min = Infinity, max = -Infinity, filled = 0;
    for (let i = 0; i < total; i++) {
      const v = I[i];
      if (!Number.isFinite(v)) { I[i] = NaN; continue; }
      filled++;
      if (v < min) min = v;
      if (v > max) max = v;
    }
    if (!filled) throw new Error("Intensity list: no finite intensities");
    let sigma = null, sigmaConstant = null;
    if (S2 && !sigmaInvalid && sigmaMin === sigmaMax) {
      sigmaConstant = sigmaMin;
    } else if (S2) {
      for (let i = 0; i < total; i++) S2[i] = hits[i] && Number.isFinite(I[i]) ? Math.sqrt(S2[i]) / hits[i] : NaN;
      sigma = S2;
    }
    return { h: ah.values, k: ak.values, l: al.values, shape, I, min, max, points: rows.count, filled, masked: total - filled, duplicates, sigma, sigmaConstant, sigmaInvalid };
  }

  function parseText(text) {
    const parser = createParser();
    parser.push(String(text || ""));
    return buildGrid(parser.finish());
  }

  global.IntensityListReader = Object.freeze({ detect, detectConfig, parseConfig, createParser, buildGrid, parseText });
})(typeof window !== "undefined" ? window : globalThis);
