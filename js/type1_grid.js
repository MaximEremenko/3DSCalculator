"use strict";

// Maps uniform h,k,l grids onto type-1 NUFFTs. For grid points
// h = h_c + m*dh with centred mode indices m,
//   sum_j c_j exp(2*pi*i h.u_j) = sum_j c'_j exp(i m.x_j),
//   x_j = wrap(2*pi*dh*u_j),  c'_j = c_j exp(2*pi*i h_c.u_j),
// exactly, because m is an integer. u_j are parent-cell fractional
// coordinates. Every phase is formed in f64 here, so a GPU transform only
// sees |x| <= pi.
(function (global) {
  const TWO_PI = 2 * Math.PI;

  // {start, step, count} of a uniform axis, or null when it is not uniform.
  function uniformAxis(values) {
    const n = values.length;
    const start = Number(values[0]);
    if (n === 1) return { start, step: 1, count: 1 };
    const step = (Number(values[n - 1]) - start) / (n - 1);
    if (!(step > 0)) return null;
    const tol = 1e-9 * Math.max(1, Math.abs(step), Math.abs(start));
    for (let i = 1; i < n - 1; i++) {
      if (Math.abs(Number(values[i]) - (start + i * step)) > tol) return null;
    }
    return { start, step, count: n };
  }

  function gridFromAxes(h, k, l) {
    const axes = [uniformAxis(h), uniformAxis(k), uniformAxis(l)];
    return axes.every(Boolean) ? axes : null;
  }

  // u = Bq.r / (2*pi) for Cartesian points r, so that q.r = 2*pi*h.u.
  function fractionalCoordinates(packed, Bq) {
    const n = packed.length / 3;
    const u = new Float64Array(n * 3);
    for (let j = 0; j < n; j++) {
      const x = packed[j * 3];
      const y = packed[j * 3 + 1];
      const z = packed[j * 3 + 2];
      for (let d = 0; d < 3; d++) {
        u[j * 3 + d] = (Bq[d][0] * x + Bq[d][1] * y + Bq[d][2] * z) / TWO_PI;
      }
    }
    return u;
  }

  function wrap(t) {
    return t - TWO_PI * Math.round(t / TWO_PI);
  }

  // Point-major type-1 coordinates with the l axis first, so that the output
  // mode index il + nl*(ik + nk*ih) equals the calculator's C-ordered index.
  function type1Points(u, axes, output) {
    const n = u.length / 3;
    const out = output || new Float32Array(n * 3);
    const [ah, ak, al] = axes;
    for (let j = 0; j < n; j++) {
      out[j * 3] = wrap(TWO_PI * al.step * u[j * 3 + 2]);
      out[j * 3 + 1] = wrap(TWO_PI * ak.step * u[j * 3 + 1]);
      out[j * 3 + 2] = wrap(TWO_PI * ah.step * u[j * 3]);
    }
    return out;
  }

  // Grid point at mode 0 of a transform with `modes` = [nl, nk, planes]
  // modes (dimension 0 = l) whose first h-plane is ih0. Mode counts may
  // exceed the axis lengths: the extra modes extend each axis upwards, so
  // output index i still maps to axis point i.
  function slabCentre(axes, ih0, modes) {
    const [ah, ak, al] = axes;
    return [
      ah.start + (ih0 + Math.floor(modes[2] / 2)) * ah.step,
      ak.start + Math.floor(modes[1] / 2) * ak.step,
      al.start + Math.floor(modes[0] / 2) * al.step,
    ];
  }

  // Interleaved f32 strengths c_j exp(2*pi*i hc.u_j); the phase is reduced
  // to one turn in f64 before the trigonometry.
  function phasedStrengths(u, strengths, hc, out) {
    const n = u.length / 3;
    const res = out || new Float32Array(n * 2);
    for (let j = 0; j < n; j++) {
      const t = hc[0] * u[j * 3] + hc[1] * u[j * 3 + 1] + hc[2] * u[j * 3 + 2];
      const p = TWO_PI * (t - Math.round(t));
      const c = Math.cos(p);
      const s = Math.sin(p);
      const re = strengths[j * 2];
      const im = strengths[j * 2 + 1];
      res[j * 2] = re * c - im * s;
      res[j * 2 + 1] = re * s + im * c;
    }
    return res;
  }

  // Exact sum_j c_j exp(2*pi*i h.u_j) at one grid point, for spot checks.
  function directAmplitude(u, strengths, hkl) {
    let re = 0;
    let im = 0;
    for (let j = 0; j < u.length / 3; j++) {
      const t = hkl[0] * u[j * 3] + hkl[1] * u[j * 3 + 1] + hkl[2] * u[j * 3 + 2];
      const p = TWO_PI * (t - Math.round(t));
      const c = Math.cos(p);
      const s = Math.sin(p);
      re += strengths[j * 2] * c - strengths[j * 2 + 1] * s;
      im += strengths[j * 2] * s + strengths[j * 2 + 1] * c;
    }
    return [re, im];
  }

  // wgpuNUFFT fine-grid length for n modes: the smallest even length
  // >= max(ceil(sigma*n), 2w) whose factors are all in {2,3,5,7,11,13}.
  // wgpuNUFFT exponential-of-semicircle kernel width for eps and sigma.
  function kernelWidth(sigma, eps) {
    const w = sigma === 2
      ? Math.ceil(Math.log10(10 / eps))
      : Math.ceil(Math.log(1 / eps) / (Math.PI * Math.sqrt(1 - 1 / sigma)));
    return Math.max(2, Math.min(16, w));
  }

  function fineGridLength(n, sigma, eps) {
    let len = Math.max(Math.ceil(sigma * n), 2 * kernelWidth(sigma, eps));
    for (;; len++) {
      if (len % 2) continue;
      let r = len;
      for (const f of [2, 3, 5, 7, 11, 13]) while (r % f === 0) r /= f;
      if (r === 1) return len;
    }
  }

  global.DiffuseType1 = Object.freeze({
    uniformAxis,
    gridFromAxes,
    fractionalCoordinates,
    type1Points,
    slabCentre,
    phasedStrengths,
    directAmplitude,
    kernelWidth,
    fineGridLength,
  });
})(typeof window !== "undefined" ? window : globalThis);
