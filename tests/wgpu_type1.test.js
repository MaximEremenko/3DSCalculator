"use strict";

const test = require("node:test");
const assert = require("node:assert/strict");
const fs = require("node:fs");
const vm = require("node:vm");
const path = require("node:path");
const { loadScripts, ROOT } = require("./helpers/browser_env");

const MiB = 1024 * 1024;
const axis = (n) => Array.from({ length: n }, (_, i) => -1.25 + i * 0.125);

// Models the public GPU API, including explicit ownership. Small transforms
// use the type-1 definition, so slab/copy changes are checked numerically.
function fixture(nav = {}, bindingLimit = 128 * MiB, failure = null) {
  const ctx = loadScripts(["js/type1_grid.js", "js/wgpu_type1.js"]);
  ctx.navigator = nav;
  const log = { plans: [], uploads: 0, frees: 0, devices: 0, live: new Set() };
  function buffer(words) {
    const b = { words, free() { assert.ok(log.live.delete(b), "buffer freed exactly once"); } };
    log.live.add(b);
    return b;
  }
  ctx.wgpuWeb = { load: async () => ({
    WebNufftModeOrder: { Centered: 0 }, WebFftPrecision: { F32: 0 },
    WgpuFft: { init: async () => {
      log.devices++;
      return {
        maxStorageBufferBindingSize: bindingLimit, maxBufferSize: bindingLimit,
        upload(words) { log.uploads++; return buffer(words.slice()); },
        createBuffer(bytes) {
          if (failure === "output") throw new Error("output allocation failed");
          return buffer(new Float32Array(bytes / 4));
        },
        download: async (b) => new Uint8Array(b.words.buffer),
        createNufftType1Plan: async (modes) => {
          if (failure === "sizing") {
            log.plans.push({ shape: Array.from(modes) });
            throw new Error("sizing captured");
          }
          if (failure === "plan") throw new Error("plan allocation failed");
          const shape = Array.from(modes);
          const p = {
            outputBytes: shape.reduce((a, b) => a * b, 8), freed: false,
            free() { assert.equal(p.freed, false); p.freed = true; },
            execute: async (points, strengths, output) => {
              if (failure === "execute") throw new Error("device lost");
              const [nl, nk, nh] = shape;
              for (let ih = 0; ih < nh; ih++) for (let ik = 0; ik < nk; ik++) for (let il = 0; il < nl; il++) {
                const m = [il - Math.floor(nl / 2), ik - Math.floor(nk / 2), ih - Math.floor(nh / 2)];
                let re = 0, im = 0;
                for (let j = 0; j < strengths.words.length / 2; j++) {
                  const phase = m.reduce((s, v, d) => s + v * points.words[j * 3 + d], 0);
                  const c = Math.cos(phase), s = Math.sin(phase);
                  re += strengths.words[j * 2] * c - strengths.words[j * 2 + 1] * s;
                  im += strengths.words[j * 2] * s + strengths.words[j * 2 + 1] * c;
                }
                const i = 2 * (il + nl * (ik + nk * ih));
                output.words[i] = re; output.words[i + 1] = im;
              }
            },
          };
          log.plans.push({ shape, plan: p });
          return p;
        },
        free() { log.frees++; assert.equal(log.live.size, 0, "device freed after caller buffers"); },
      };
    } },
  }) };
  return { ctx, log };
}

function spec(shape, count = 3) {
  const [h, k, l] = shape.map(axis);
  const sourcesPacked = new Float64Array(count * 3), strengths = new Float64Array(count * 2);
  for (let i = 0; i < count; i++) {
    sourcesPacked.set([i * 0.23, i * 0.37, -i * 0.61], i * 3);
    strengths.set([1 + i / 4, -0.13 * i], i * 2);
  }
  return { isign: 1, sourcesPacked, strengths,
    grid: { h, k, l, Bq: [[2 * Math.PI, 0, 0], [0, 2 * Math.PI, 0], [0, 0, 2 * Math.PI]],
      start: 0, count: h.length * k.length * l.length } };
}

test("mobile policy includes desktop-mode iPads and devices reporting low RAM", () => {
  for (const nav of [{ userAgent: "iPhone" }, { userAgent: "Android" },
    { platform: "MacIntel", maxTouchPoints: 5 }, { deviceMemory: 4 }, { userAgentData: { mobile: true } }]) {
    const { ctx } = fixture(nav);
    assert.equal(ctx.WgpuType1.memoryPolicy().lowMemory, true);
    assert.equal(ctx.WgpuType1.memoryPolicy().maxChunkPoints, 120000);
  }
  assert.equal(fixture({ userAgent: "Chrome", deviceMemory: 8 }).ctx.WgpuType1.memoryPolicy().lowMemory, false);
});

test("mobile slabs fit the total budget even when the GPU permits a much larger buffer", async () => {
  const { ctx, log } = fixture({ userAgent: "iPhone" }, 768 * MiB, "sizing");
  const s = spec([81, 81, 81], 1);
  // Stop after sizing, before allocating or computing the large output.
  await assert.rejects(ctx.WgpuType1.run(s, { tol: 1e-6, upsampfac: 2 }), /sizing captured/);
  const modes = log.plans[0].shape;
  assert.ok(modes[2] < 81, "large hardware limit must not prevent slabbing on mobile");
  const fineBytes = modes.reduce((bytes, n) => bytes * ctx.DiffuseType1.fineGridLength(n, 2, 1e-6), 8);
  assert.ok(fineBytes <= 16 * MiB);
});

test("the standard PMN 161-cubed grid can run in mobile chunks", async () => {
  const { ctx, log } = fixture({ userAgent: "iPhone" }, 128 * MiB, "sizing");
  const s = spec([161, 161, 161], 1);
  s.pointCapacity = 320000;
  s.grid.count = 3 * 161 * 161;
  await assert.rejects(ctx.WgpuType1.run(s, { tol: 1e-6, upsampfac: 2 }), /sizing captured/);
  assert.equal(log.plans.length, 1, "the PMN chunk must pass preflight");
});

test("desktop sizing budgets the whole workspace rather than one fine-grid buffer", async () => {
  const { ctx, log } = fixture({}, 768 * MiB, "sizing");
  await assert.rejects(ctx.WgpuType1.run(spec([161, 161, 161], 1), {}), /sizing captured/);
  assert.ok(log.plans[0].shape[2] < 161);
});

test("oversized point sets fail before allocation even on high-limit GPUs", async () => {
  const { ctx, log } = fixture({ userAgent: "iPhone" }, 768 * MiB);
  const s = spec([2, 3, 4]);
  s.pointCapacity = 2000000;
  await assert.rejects(ctx.WgpuType1.run(s, {}), /memory budget/);
  assert.equal(log.plans.length, 0);
  assert.equal(log.uploads, 0);
});

test("checks the padded minimum before uploading points or allocating a plan", async () => {
  const { ctx, log } = fixture({}, 128 * 1024);
  await assert.rejects(ctx.WgpuType1.run(spec([1, 3, 3]), {}), /memory budget/);
  assert.equal(log.plans.length, 0);
  assert.equal(log.uploads, 0);
  assert.equal(log.frees, 1);
});

test("slab splitting preserves complex amplitudes and partial final slabs", async () => {
  const { ctx, log } = fixture({ userAgent: "iPhone" }, 30000);
  const s = spec([17, 7, 7]);
  // Multiple calls also exercise a chunk starting away from h=0.
  s.grid.start = 2 * 49; s.grid.count = 15 * 49;
  const out = await ctx.WgpuType1.run(s, { tol: 1e-6, upsampfac: 2 });
  assert.ok(ctx.WgpuType1.stats.lastSlabs > 1);
  const u = ctx.DiffuseType1.fractionalCoordinates(s.sourcesPacked, s.grid.Bq);
  for (let p = 0; p < 15; p++) for (let ik = 0; ik < 7; ik++) for (let il = 0; il < 7; il++) {
    const want = ctx.DiffuseType1.directAmplitude(u, s.strengths, [s.grid.h[p + 2], s.grid.k[ik], s.grid.l[il]]);
    const i = 2 * (p * 49 + ik * 7 + il);
    assert.ok(Math.hypot(out[i] - want[0], out[i + 1] - want[1]) < 1e-5);
  }
  await ctx.WgpuType1.reset();
  assert.equal(log.live.size, 0);
  assert.equal(log.frees, 1);
});

test("mobile changes of plan shape release the old device before creating the next", async () => {
  const { ctx, log } = fixture({ userAgent: "iPhone" });
  await ctx.WgpuType1.run(spec([2, 3, 4]), {});
  await ctx.WgpuType1.run(spec([3, 4, 5]), {});
  assert.equal(log.devices, 2);
  assert.equal(log.frees, 1);
  assert.equal(log.plans[0].plan.freed, true);
  await ctx.WgpuType1.reset();
  assert.equal(log.frees, 2);
});

test("allocation failures and device loss release all owned buffers and plans", async () => {
  for (const failure of ["plan", "output", "execute"]) {
    const { ctx, log } = fixture({ userAgent: "iPhone" }, 128 * MiB, failure);
    await assert.rejects(ctx.WgpuType1.run(spec([2, 3, 4]), {}), /failed|device lost/);
    assert.equal(log.live.size, 0);
    assert.equal(log.frees, 1);
    for (const p of log.plans) assert.equal(p.plan.freed, true);
  }
});

// Exercise the actual inline app functions without booting its UI or CDN
// dependencies. Mobile chunking and cleanup must be wired into compute(),
// not just implemented by the adapter.
function appFixture() {
  const { ctx, log } = fixture({ userAgent: "iPhone" });
  const text = fs.readFileSync(path.join(ROOT, "index.html"), "utf8");
  // Extract complete function declarations, using their next top-level
  // declaration as the boundary (these functions have no nested names).
  function source(name, next) {
    const start = text.indexOf(`function ${name}(`);
    return text.slice(start, text.indexOf(`\n${next}`, start));
  }
  const values = { backend: "wgpu_type1", gridLimit: "8000000", hMin: "-1", hMax: "1", hStep: ".1",
    kMin: "-1", kMax: "1", kStep: ".1", lMin: "-1", lMax: "1", lStep: ".1",
    tol: "1e-6", upsamp: "2", subAvg: "0", normMode: "none", smooth3d: "0", smoothType: "lanczos",
    chebDb: "100", smoothScale: "1", magMode: "nuclear", magQ0: "average", magPara: "0", laueSym: "none" };
  const elements = {};
  ctx.$ = (id) => elements[id] || (elements[id] = { value: values[id] || "0" });
  ctx.state = { rmc: { atoms: 3, Bq: [], file: "test", unknown: [] }, render3dEpoch: 0 };
  ctx.BACKEND_PRESETS = { wgpu_type1: { tol: "1e-6", upsamp: "2.0", gridLimit: "8000000" } };
  ctx.EXACT_SUM_WORK_LIMIT = 2e8;
  for (const name of ["status", "capture3dCameras", "getIncludedElements", "getScatteringConfigFromUi", "needDiffuseAmplitude",
    "syncResultControls", "updateSliceUI", "applyPendingPlane", "draw", "schedule3dRender", "renderCompareUi"]) ctx[name] = () => {};
  ctx.filterParsedByElements = (r) => r;
  ctx.getScatteringConfigFromUi = () => ({ type: "neutron" });
  ctx.axis = () => axis(81);
  ctx.ensembleConfigs = () => [ctx.state.rmc];
  ctx.fmtMs = () => "0";
  ctx.planeSig = () => "test";
  ctx.runType3 = () => {};
  vm.runInContext(source("applyBackendDefaults", "function "), ctx);
  vm.runInContext("async " + source("compute", "function slice("), ctx);
  return { ctx, log };
}

test("mobile app applies smaller defaults and caps chunks even if the limit was raised", async () => {
  const { ctx } = appFixture();
  ctx.applyBackendDefaults("wgpu_type1");
  assert.equal(ctx.$("gridLimit").value, "120000");
  ctx.$("gridLimit").value = "8000000";
  let chunkSize;
  ctx.DiffuseAmplitude = { computeIntensity: async (args) => {
    chunkSize = args.chunkSize;
    return { I: new Float64Array(81 ** 3), min: 0, max: 0, timings: {}, profile: {} };
  } };
  await ctx.compute();
  assert.equal(chunkSize, 18 * 81 * 81);
  assert.ok(chunkSize <= 120000);
});

test("mobile app releases type-1 resources on success and computation failure", async () => {
  for (const fail of [false, true]) {
    const { ctx } = appFixture();
    let resets = 0;
    ctx.WgpuType1 = { memoryPolicy: () => ({ lowMemory: true, maxChunkPoints: 120000, maxGridPoints: 8000000 }),
      reset: async () => { resets++; }, stats: {} };
    ctx.DiffuseAmplitude = { computeIntensity: async () => {
      if (fail) throw new Error("compute failed");
      return { I: new Float64Array(81 ** 3), min: 0, max: 0, timings: {}, profile: {} };
    } };
    if (fail) await assert.rejects(ctx.compute(), /compute failed/);
    else await ctx.compute();
    assert.equal(resets, 2);
  }
});

test("mobile app refuses a volume that cannot fit before starting any transform", async () => {
  const { ctx } = appFixture();
  ctx.axis = () => axis(201);
  let calls = 0;
  ctx.DiffuseAmplitude = { computeIntensity: async () => { calls++; } };
  await assert.rejects(ctx.compute(), /mobile memory limit/);
  assert.equal(calls, 0);
});
