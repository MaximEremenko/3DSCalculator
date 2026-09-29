"use strict";

(function (global) {
  const NEUTRON = {
    va: 0,
    h: -0.3739,
    d: 0.6671,
    he: 0.326,
    li: -0.19,
    "7l": -0.222,
    be: 0.779,
    b: 0.53,
    c: 0.6646,
    n: 0.936,
    o: 0.5803,
    f: 0.5654,
    ne: 0.4566,
    na: 0.363,
    mg: 0.5375,
    al: 0.3449,
    si: 0.41491,
    p: 0.513,
    s: 0.2847,
    cl: 0.9577,
    ar: 0.1909,
    k: 0.367,
    ca: 0.47,
    sc: 1.229,
    ti: -0.3438,
    v: -0.03824,
    cr: 0.3635,
    mn: -0.373,
    fe: 0.945,
    co: 0.249,
    ni: 1.03,
    cu: 0.7718,
    zn: 0.568,
    ga: 0.7288,
    ge: 0.8185,
    as: 0.658,
    se: 0.797,
    br: 0.6795,
    kr: 0.781,
    rb: 0.709,
    sr: 0.702,
    y: 0.775,
    zr: 0.716,
    nb: 0.7054,
    mo: 0.6715,
    tc: 0.68,
    ru: 0.703,
    rh: 0.588,
    pd: 0.591,
    ag: 0.5922,
    cd: 0.487,
    in: 0.4065,
    sn: 0.6225,
    sb: 0.557,
    te: 0.58,
    i: 0.528,
    xe: 0.492,
    cs: 0.542,
    ba: 0.507,
    la: 0.824,
    ce: 0.484,
    pr: 0.458,
    nd: 0.769,
    pm: 1.26,
    sm: 0.08,
    eu: 0.722,
    gd: 0.65,
    tb: 0.738,
    dy: 1.69,
    ho: 0.801,
    er: 0.779,
    tm: 0.707,
    yb: 1.243,
    lu: 0.721,
    hf: 0.777,
    ta: 0.691,
    w: 0.486,
    re: 0.92,
    os: 1.07,
    ir: 1.06,
    pt: 0.96,
    au: 0.763,
    hg: 1.2692,
    tl: 0.8776,
    pb: 0.9405,
    bi: 0.8532,
    po: 0,
    at: 0,
    rn: 0,
    fr: 0,
    ra: 1,
    ac: 0,
    th: 1.031,
    pa: 0.91,
    u: 0.8417,
    np: 1.055,
    pu: 0,
    am: 0.83,
    cm: 0,
  };

  function mul(v, m) {
    return [
      v[0] * m[0][0] + v[1] * m[1][0] + v[2] * m[2][0],
      v[0] * m[0][1] + v[1] * m[1][1] + v[2] * m[2][1],
      v[0] * m[0][2] + v[1] * m[1][2] + v[2] * m[2][2],
    ];
  }

  function interleave3(x, y, z) {
    const n = x.length;
    const out = new Float64Array(n * 3);
    for (let i = 0; i < n; i++) {
      out[i * 3] = x[i];
      out[i * 3 + 1] = y[i];
      out[i * 3 + 2] = z[i];
    }
    return out;
  }

  function complexReal(a) {
    const out = new Float64Array(a.length * 2);
    for (let i = 0; i < a.length; i++) out[i * 2] = a[i];
    return out;
  }

  function complexOnes(atomCount) {
    const out = new Float64Array(Math.max(0, atomCount) * 2);
    for (let i = 0; i < atomCount; i++) out[i * 2] = 1;
    return out;
  }

  function targetsChunk(h, k, l, Bq, start, count) {
    const nh = h.length;
    const nk = k.length;
    const nl = l.length;
    const total = nh * nk * nl;
    if (start < 0 || count < 0 || start + count > total) {
      throw new Error(
        `Target chunk out of range: start=${start}, count=${count}, total=${total}`
      );
    }
    const t = new Float64Array(count * 3);
    const plane = nk * nl;
    for (let q = 0; q < count; q++) {
      const linear = start + q;
      const ih = Math.floor(linear / plane);
      const rem = linear - ih * plane;
      const ik = Math.floor(rem / nl);
      const il = rem - ik * nl;
      const v = mul([h[ih], k[ik], l[il]], Bq);
      t[q * 3] = v[0];
      t[q * 3 + 1] = v[1];
      t[q * 3 + 2] = v[2];
    }
    return t;
  }

  // |q| for a run of the C-ordered grid, without building the target array.
  function qMagnitudesChunk(h, k, l, Bq, start, count) {
    const nk = k.length;
    const nl = l.length;
    const plane = nk * nl;
    const out = new Float64Array(count);
    for (let q = 0; q < count; q++) {
      const linear = start + q;
      const ih = Math.floor(linear / plane);
      const rem = linear - ih * plane;
      const ik = Math.floor(rem / nl);
      const il = rem - ik * nl;
      const v = mul([h[ih], k[ik], l[il]], Bq);
      out[q] = Math.hypot(v[0], v[1], v[2]);
    }
    return out;
  }

  function attachNeutronCoefficients(parsed, fallback = 10) {
    if (!parsed || !Array.isArray(parsed.elements)) {
      throw new Error("Invalid parsed RMC structure");
    }
    const n = parsed.elements.length;
    const fca = new Float64Array(n);
    const present = new Set();
    const unknown = new Set();
    for (let i = 0; i < n; i++) {
      const elem = String(parsed.elements[i] || "").trim();
      const key = elem.toLowerCase();
      present.add(elem);
      if (Object.prototype.hasOwnProperty.call(NEUTRON, key)) {
        fca[i] = NEUTRON[key];
      } else {
        fca[i] = fallback;
        unknown.add(elem);
      }
    }
    return {
      ...parsed,
      fca,
      present: [...present].sort(),
      unknown: [...unknown].sort(),
    };
  }

  // Lattice sum over the supercell cells along one reciprocal axis,
  // sum_{c=c0}^{c0+n-1} exp(+2*pi*i*h*c). Only the fractional part of h
  // matters because c is an integer, which keeps integer h exact (= n).
  function laueAxis(values, n, c0) {
    const out = new Float64Array(values.length * 2);
    for (let i = 0; i < values.length; i++) {
      const e = Number(values[i]) - Math.round(Number(values[i]));
      const s = Math.sin(Math.PI * e);
      const magnitude = s === 0 ? n : Math.sin(Math.PI * n * e) / s;
      const phase = Math.PI * e * (n - 1 + 2 * c0);
      out[i * 2] = magnitude * Math.cos(phase);
      out[i * 2 + 1] = magnitude * Math.sin(phase);
    }
    return out;
  }

  // Product of the per-axis lattice sums for a run of the C-ordered h,k,l grid.
  function laueChunk(Lh, Lk, Ll, start, count) {
    const nk = Lk.length / 2;
    const nl = Ll.length / 2;
    const plane = nk * nl;
    const out = new Float64Array(count * 2);
    for (let q = 0; q < count; q++) {
      const linear = start + q;
      const ih = Math.floor(linear / plane);
      const rem = linear - ih * plane;
      const ik = Math.floor(rem / nl);
      const il = rem - ik * nl;
      const ar = Lh[ih * 2];
      const ai = Lh[ih * 2 + 1];
      const br = Lk[ik * 2];
      const bi = Lk[ik * 2 + 1];
      const cr = Ll[il * 2];
      const ci = Ll[il * 2 + 1];
      const abr = ar * br - ai * bi;
      const abi = ar * bi + ai * br;
      out[q * 2] = abr * cr - abi * ci;
      out[q * 2 + 1] = abr * ci + abi * cr;
    }
    return out;
  }

  // First cell index along each axis, after checking that the per-atom
  // indices (3 per atom) fit inside the supercell.
  function cellOrigin(cellIndex, atoms, supercell) {
    const lo = [Infinity, Infinity, Infinity];
    const hi = [-Infinity, -Infinity, -Infinity];
    for (let i = 0; i < atoms; i++) {
      for (let d = 0; d < 3; d++) {
        const c = cellIndex[i * 3 + d];
        if (c < lo[d]) lo[d] = c;
        if (c > hi[d]) hi[d] = c;
      }
    }
    for (let d = 0; d < 3; d++) {
      if (!(hi[d] - lo[d] < supercell[d])) {
        throw new Error(
          `Cell indices along axis ${d + 1} span ${hi[d] - lo[d] + 1} cells, more than the supercell (${supercell[d]})`
        );
      }
    }
    return lo;
  }

  // Neutron magnetic scattering length per Bohr magneton, gamma r_e / 2, in
  // 1e-12 cm (CODATA 2018: gamma = 1.91304273, r_e = 2.8179403262 fm); its
  // square is 0.0726529 barn, as in Scatty (RMCProfile uses 0.269536639).
  const MAGNETIC_LENGTH = (1.91304273 * 2.8179403262) / 2 / 10;

  // Empty accumulators for |dM|^2 and qhat.dM over a chunk.
  function magneticAccumulator(count) {
    return { t: new Float64Array(count), s: new Float64Array(count * 2), zero: new Uint8Array(count) };
  }

  // Adds Cartesian component `alpha` of a magnetic fluctuation amplitude dm
  // (interleaved complex) to |dM|^2 and to qhat.dM for a run of the grid.
  function accumulateMagneticComponent(acc, dm, alpha, h, k, l, Bq, start, count) {
    const nk = k.length;
    const nl = l.length;
    const plane = nk * nl;
    for (let i = 0; i < count; i++) {
      const linear = start + i;
      const ih = Math.floor(linear / plane);
      const rem = linear - ih * plane;
      const ik = Math.floor(rem / nl);
      const il = rem - ik * nl;
      const hv = h[ih];
      const kv = k[ik];
      const lv = l[il];
      const qx = hv * Bq[0][0] + kv * Bq[1][0] + lv * Bq[2][0];
      const qy = hv * Bq[0][1] + kv * Bq[1][1] + lv * Bq[2][1];
      const qz = hv * Bq[0][2] + kv * Bq[1][2] + lv * Bq[2][2];
      const qn = Math.hypot(qx, qy, qz);
      const w = qn > 1e-12 ? (alpha === 0 ? qx : alpha === 1 ? qy : qz) / qn : 0;
      if (!(qn > 1e-12)) acc.zero[i] = 1;
      const re = dm[i * 2];
      const im = dm[i * 2 + 1];
      acc.t[i] += re * re + im * im;
      acc.s[i * 2] += w * re;
      acc.s[i * 2 + 1] += w * im;
    }
  }

  // Writes, or with `add` adds, |dM_perp|^2 = |dM|^2 - |qhat.dM|^2 into out.
  // At Q = 0 the direction is undefined: 0 as in RMCProfile, or with
  // `qZeroAverage` the orientational average (2/3)|dM|^2 as in Scatty.
  // `offset` (optional, per point) is subtracted wherever the magnetic
  // intensity is defined, e.g. the ideal-paramagnet term.
  function finishMagneticChunk(out, start, count, acc, qZeroAverage, add, offset) {
    let min = Infinity;
    let max = -Infinity;
    for (let i = 0; i < count; i++) {
      const defined = !acc.zero[i] || qZeroAverage;
      let perp = acc.zero[i]
        ? (qZeroAverage ? (2 / 3) * acc.t[i] : 0)
        : Math.max(0, acc.t[i] - acc.s[i * 2] ** 2 - acc.s[i * 2 + 1] ** 2);
      if (offset && defined) perp -= offset[i];
      const v = add ? out[start + i] + perp : perp;
      out[start + i] = v;
      if (v < min) min = v;
      if (v > max) max = v;
    }
    return { min, max };
  }

  function accumulateIntensityChunk(
    out,
    start,
    count,
    q,
    qa,
    qd,
    mInv,
    subtractAverage
  ) {
    let min = Infinity;
    let max = -Infinity;
    for (let i = 0; i < count; i++) {
      const ar = q[i * 2];
      const ai = q[i * 2 + 1];
      let br = 0;
      let bi = 0;
      if (subtractAverage) {
        const avr = qa[i * 2];
        const avi = qa[i * 2 + 1];
        const dr = qd[i * 2];
        const di = qd[i * 2 + 1];
        br = (avr * dr - avi * di) * mInv;
        bi = (avr * di + avi * dr) * mInv;
      }
      const rr = ar - br;
      const ii = ai - bi;
      const v = rr * rr + ii * ii;
      out[start + i] = v;
      if (v < min) min = v;
      if (v > max) max = v;
    }
    return { min, max };
  }

  global.DiffuseCore = Object.freeze({
    NEUTRON,
    interleave3,
    complexReal,
    complexOnes,
    targetsChunk,
    qMagnitudesChunk,
    attachNeutronCoefficients,
    laueAxis,
    laueChunk,
    cellOrigin,
    MAGNETIC_LENGTH,
    magneticAccumulator,
    accumulateMagneticComponent,
    finishMagneticChunk,
    accumulateIntensityChunk,
  });
})(typeof window !== "undefined" ? window : globalThis);
