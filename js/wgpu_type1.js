"use strict";

// WebGPU type-1 NUFFT backend on wgpu-web (wgpuNUFFT, vendored as
// js/wgpu_web.js). It evaluates whole h-planes of a uniform h,k,l grid with
// the mapping in js/type1_grid.js; points are uploaded once per source set
// and plans are reused between calls.
(function (global) {
  const PLAN_CACHE = 3;
  const POINT_CACHE = 4;
  const FINE_GRID_BUDGET = 768 * 1024 * 1024; // bytes per fine-grid buffer (a plan holds two)
  const CHECK_POINTS = 4;
  const CHECK_LIMIT = 1e-3;
  // Short axes are padded to max(2B, B + 2w, B + w + 5) fine-grid cells,
  // B = [16, 16, 8] for the (l, k, h) dimensions, where wgpu-web 0.2 used
  // its block spreader (its per-cell gather stalled on clustered points).
  // wgpu-web 0.3 spreads every grid alike; the padding now gives small grids,
  // such as the lattice grids of A_delta, one shape whose shaders compile
  // once, instead of once per shape.
  const BLOCK = [16, 16, 8];

  let contextPromise = null;
  const plans = new Map();
  const pointSets = new Map();
  const stats = { lastCheckError: null, lastSlabs: 0, lastModes: null, adapter: null };

  function mapping() {
    if (!global.DiffuseType1) throw new Error("DiffuseType1 missing. Load ./js/type1_grid.js first.");
    return global.DiffuseType1;
  }

  async function context() {
    if (!contextPromise) {
      contextPromise = (async () => {
        if (!global.wgpuWeb) throw new Error("wgpu-web is not loaded (./js/wgpu_web.js).");
        const api = await global.wgpuWeb.load();
        const gpu = await api.WgpuFft.init();
        stats.adapter = gpu.adapterName || gpu.backend || "WebGPU";
        return { api, gpu };
      })().catch((error) => {
        contextPromise = null;
        throw error;
      });
    }
    return contextPromise;
  }

  function freeAll(cache) {
    for (const entry of cache.values()) entry.free();
    cache.clear();
  }

  function remember(cache, key, entry, limit) {
    cache.set(key, entry);
    while (cache.size > limit) {
      const [oldKey, oldEntry] = cache.entries().next().value;
      cache.delete(oldKey);
      oldEntry.free();
    }
    return entry;
  }

  // Drops every GPU object; the next call starts a fresh WebGPU context.
  function reset() {
    freeAll(plans);
    freeAll(pointSets);
    if (contextPromise) {
      contextPromise.then(({ gpu }) => { try { gpu.free(); } catch (_) {} }).catch(() => {});
    }
    contextPromise = null;
  }

  // Fractional coordinates and uploaded type-1 points of one source set,
  // cached by the array object; they depend only on the grid steps. Points
  // beyond the source count up to `capacity` are zero, so sources of
  // different sizes (element groups) can share one plan.
  function pointSet(ctx, packed, Bq, axes, capacity) {
    const key = `${capacity}|${axes[0].step}|${axes[1].step}|${axes[2].step}|${JSON.stringify(Bq)}`;
    const cached = pointSets.get(packed);
    if (cached && cached.key === key) return cached;
    if (cached) {
      pointSets.delete(packed);
      cached.free();
    }
    const T1 = mapping();
    const u = T1.fractionalCoordinates(packed, Bq);
    const points = new Float32Array(capacity * 3);
    points.set(T1.type1Points(u, axes));
    const buffer = ctx.gpu.upload(points);
    return remember(pointSets, packed, {
      key,
      u,
      buffer,
      count: u.length / 3,
      free() { buffer.free(); },
    }, POINT_CACHE);
  }

  async function plan(ctx, count, modes, eps, sigma) {
    const key = `${count}|${modes.join("x")}|${eps}|${sigma}`;
    const cached = plans.get(key);
    if (cached) {
      plans.delete(key);
      plans.set(key, cached);
      return cached;
    }
    const p = await ctx.gpu.createNufftType1Plan(
      new Uint32Array(modes), count, 1, eps, 1,
      ctx.api.WebNufftModeOrder.Centered, sigma, ctx.api.WebFftPrecision.F32
    );
    const output = ctx.gpu.createBuffer(Number(p.outputBytes));
    return remember(plans, key, {
      plan: p,
      output,
      free() { output.free(); p.free(); },
    }, PLAN_CACHE);
  }

  // Smallest mode count per (l, k, h) dimension; see BLOCK.
  function minimumModes(sigma, eps) {
    const T1 = mapping();
    const w = T1.kernelWidth(sigma, eps);
    return BLOCK.map((b) => {
      const fine = Math.max(2 * b, b + 2 * w, b + w + 5);
      let n = 1;
      while (T1.fineGridLength(n, sigma, eps) < fine) n++;
      return n;
    });
  }

  // Largest number of h-planes per transform whose fine grid fits a buffer.
  function planesPerSlab(ctx, nl, nk, planes, sigma, eps) {
    const T1 = mapping();
    const plane = 8 * T1.fineGridLength(nl, sigma, eps) * T1.fineGridLength(nk, sigma, eps);
    const limit = Math.min(
      FINE_GRID_BUDGET,
      Number(ctx.gpu.maxStorageBufferBindingSize) || 134217728,
      Number(ctx.gpu.maxBufferSize) || 268435456
    );
    let slab = planes;
    while (slab > 1 && plane * T1.fineGridLength(slab, sigma, eps) > limit) slab = Math.ceil(slab / 2);
    if (plane * T1.fineGridLength(slab, sigma, eps) > limit) {
      throw new Error("One h-plane of this grid needs more GPU memory than a single buffer allows; reduce the k or l range.");
    }
    return slab;
  }

  // Compares a few outputs with exact sums, relative to sum_j |c_j|.
  function spotCheck(u, strengths, grid, axes, out) {
    const T1 = mapping();
    let norm = 0;
    for (let j = 0; j < strengths.length / 2; j++) norm += Math.hypot(strengths[j * 2], strengths[j * 2 + 1]);
    if (!(norm > 0)) return 0;
    const plane = axes[1].count * axes[2].count;
    let worst = 0;
    for (let c = 0; c < CHECK_POINTS; c++) {
      const i = Math.floor(((c + 0.5) / CHECK_POINTS) * grid.count);
      const linear = grid.start + i;
      const ih = Math.floor(linear / plane);
      const rem = linear - ih * plane;
      const ik = Math.floor(rem / axes[2].count);
      const il = rem - ik * axes[2].count;
      const want = T1.directAmplitude(u, strengths, [grid.h[ih], grid.k[ik], grid.l[il]]);
      worst = Math.max(worst, Math.hypot(out[i * 2] - want[0], out[i * 2 + 1] - want[1]) / norm);
    }
    return worst;
  }

  // spec: { sourcesPacked (Cartesian, 3 per point), strengths (interleaved
  // complex), isign: 1, grid: { h, k, l, Bq, start, count }, pointCapacity }
  // where start and count cover whole h-planes and pointCapacity (optional)
  // is the largest source count that should share the plan. Returns
  // interleaved complex amplitudes.
  async function run(spec, opts, onStage) {
    const T1 = mapping();
    const grid = spec.grid;
    if (!grid) throw new Error("The wgpuNUFFT type-1 backend needs the h,k,l grid.");
    if (spec.isign !== 1) throw new Error("The wgpuNUFFT type-1 backend supports isign = +1 only.");
    const axes = T1.gridFromAxes(grid.h, grid.k, grid.l);
    if (!axes) throw new Error("The wgpuNUFFT type-1 backend needs uniform h, k and l axes.");
    const plane = axes[1].count * axes[2].count;
    if (grid.start % plane || grid.count % plane) {
      throw new Error("wgpuNUFFT type-1 chunks must hold whole h-planes.");
    }
    const tol = Number(opts && opts.tol);
    const eps = Number.isFinite(tol) ? Math.min(1e-2, Math.max(1e-6, tol)) : 1e-6;
    const up = Number(opts && opts.upsampfac);
    const sigma = Number.isFinite(up) && up > 1 ? up : 2;
    const stage = (name, details) => { if (typeof onStage === "function") onStage(name, details || {}); };
    const ih0 = grid.start / plane;
    const planes = grid.count / plane;

    try {
      stage("t1_setup", {});
      const ctx = await context();
      const count = spec.sourcesPacked.length / 3;
      const capacity = Math.max(count, Math.floor(Number(spec.pointCapacity) || 0));
      const source = pointSet(ctx, spec.sourcesPacked, grid.Bq, axes, capacity);
      // Pad short axes (see BLOCK); the extra modes extend each axis
      // upwards and are dropped when copying out.
      const minModes = minimumModes(sigma, eps);
      const nl = axes[2].count;
      const nk = axes[1].count;
      const ml = Math.max(nl, minModes[0]);
      const mk = Math.max(nk, minModes[1]);
      const slab = Math.max(planesPerSlab(ctx, ml, mk, planes, sigma, eps), minModes[2]);
      const modes = [ml, mk, slab];
      const p = await plan(ctx, capacity, modes, eps, sigma);
      const out = new Float64Array(grid.count * 2);
      const phased = new Float32Array(capacity * 2);
      const slabs = Math.ceil(planes / slab);
      stats.lastSlabs = slabs;
      stats.lastModes = modes;
      for (let s = 0, p0 = 0; p0 < planes; s++, p0 += slab) {
        stage("t1_phase", { slab: s + 1, slabs });
        T1.phasedStrengths(source.u, spec.strengths, T1.slabCentre(axes, ih0 + p0, modes), phased);
        const input = ctx.gpu.upload(phased);
        try {
          stage("t1_execute", { slab: s + 1, slabs, modes });
          await p.plan.execute(source.buffer, input, p.output, 1);
          stage("t1_download", { slab: s + 1, slabs });
          const bytes = await ctx.gpu.download(p.output);
          const words = bytes.byteOffset % 4 === 0
            ? new Float32Array(bytes.buffer, bytes.byteOffset, bytes.byteLength / 4)
            : new Float32Array(bytes.slice().buffer);
          const keep = Math.min(slab, planes - p0);
          for (let q = 0; q < keep; q++) {
            for (let ik = 0; ik < nk; ik++) {
              const from = 2 * ml * (ik + mk * q);
              const to = 2 * ((p0 + q) * plane + ik * nl);
              out.set(words.subarray(from, from + 2 * nl), to);
            }
          }
        } finally {
          input.free();
        }
      }
      stage("t1_check", {});
      const error = spotCheck(source.u, spec.strengths, grid, axes, out);
      stats.lastCheckError = error;
      if (!(error < CHECK_LIMIT)) {
        throw new Error(`wgpuNUFFT type-1 self-check failed: error ${error.toExponential(2)} relative to sum|c|.`);
      }
      stage("done", { outputs: grid.count });
      return out;
    } catch (error) {
      reset();
      throw error;
    }
  }

  global.WgpuType1 = Object.freeze({ run, reset, stats });
})(typeof window !== "undefined" ? window : globalThis);
