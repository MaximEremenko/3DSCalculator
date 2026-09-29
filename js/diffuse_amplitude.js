"use strict";

(function (global) {
  function needCore() {
    if (!global.DiffuseCore) {
      throw new Error("DiffuseCore missing. Load ./js/diffuse_core.js first.");
    }
    return global.DiffuseCore;
  }

  function needCoeffSource() {
    if (!global.ScatteringCoeffSource) {
      throw new Error(
        "ScatteringCoeffSource missing. Load ./js/scattering_coeff_sources.js first."
      );
    }
    return global.ScatteringCoeffSource;
  }

  function nowMs() {
    if (typeof performance !== "undefined" && performance.now) {
      return performance.now();
    }
    return Date.now();
  }

  function normElem(sym) {
    return String(sym || "").trim().toLowerCase();
  }

  function normalizeScatteringConfig(cfg) {
    const src = cfg || {};
    const type = String(src.type || "neutron").trim().toLowerCase();
    const rawTableNum = Number.isFinite(Number(src.tableNum))
      ? Math.floor(Number(src.tableNum))
      : 1;
    const out = {
      type: type === "xray" || type === "electron" ? type : "neutron",
      model: "fast",
      tableNum: Math.max(1, Math.min(5, rawTableNum)),
      defaultValence: Number.isFinite(Number(src.defaultValence))
        ? Math.trunc(Number(src.defaultValence))
        : 0,
      valenceMapText: String(src.valenceMapText || ""),
    };
    if (out.type === "neutron") {
      const model = String(src.model || "fast").trim().toLowerCase();
      out.model = model === "grouped_exact" ? "grouped_exact" : "fast";
    } else if (out.type === "xray") {
      const model = String(src.model || "waasmaier").trim().toLowerCase();
      out.model = model === "table" ? "table" : "waasmaier";
      out.tableNum = Math.max(0, Math.min(5, rawTableNum));
    } else {
      const model = String(src.model || "neutral_table").trim().toLowerCase();
      out.model = model === "ion_peng" ? "ion_peng" : "neutral_table";
    }
    return out;
  }

  function parseValenceMap(text) {
    const map = Object.create(null);
    const raw = String(text || "").trim();
    if (!raw) return map;
    const parts = raw.split(/[;,]+/);
    for (const part of parts) {
      const s = part.trim();
      if (!s) continue;
      const m = /^([A-Za-z0-9]+)\s*:\s*([+-]?\d+)$/.exec(s);
      if (!m) continue;
      const key = normElem(m[1]);
      const val = Number.parseInt(m[2], 10);
      if (!Number.isInteger(val)) continue;
      map[key] = val;
    }
    return map;
  }

  function packTripletsByIndex(a, b, c, indices) {
    const n = indices.length;
    const out = new Float64Array(n * 3);
    for (let i = 0; i < n; i++) {
      const ix = indices[i];
      const j = i * 3;
      out[j] = a[ix];
      out[j + 1] = b[ix];
      out[j + 2] = c[ix];
    }
    return out;
  }

  function groupAtomsByElement(parsed) {
    const groups = new Map();
    for (let i = 0; i < parsed.atoms; i++) {
      const raw = String(parsed.elements[i] || "").trim();
      const key = normElem(raw);
      if (!groups.has(key)) {
        groups.set(key, {
          key,
          element: raw || key,
          indices: [],
        });
      }
      groups.get(key).indices.push(i);
    }
    return Array.from(groups.values());
  }

  function makeQMagnitudes(targetsPacked) {
    const n = Math.floor(targetsPacked.length / 3);
    const out = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const j = i * 3;
      const qx = targetsPacked[j];
      const qy = targetsPacked[j + 1];
      const qz = targetsPacked[j + 2];
      out[i] = Math.hypot(qx, qy, qz);
    }
    return out;
  }

  function fillConstant(out, v) {
    for (let i = 0; i < out.length; i++) out[i] = v;
    return out;
  }

  function accumulateScaledComplex(dst, src, scale) {
    const n = scale.length;
    for (let i = 0; i < n; i++) {
      const s = scale[i];
      const j = i * 2;
      dst[j] += src[j] * s;
      dst[j + 1] += src[j + 1] * s;
    }
  }

  function describeConfig(cfg) {
    if (cfg.type === "neutron") {
      return cfg.model === "grouped_exact"
        ? "Neutron grouped exact"
        : "Neutron fast";
    }
    if (cfg.type === "xray") {
      return cfg.model === "table"
        ? `X-ray (table ${cfg.tableNum})`
        : "X-ray (Waasmaier)";
    }
    return cfg.model === "ion_peng"
      ? "Electron ion (Peng)"
      : `Electron neutral (table ${cfg.tableNum})`;
  }

  function buildGroupedEvaluator(cfg) {
    if (cfg.type === "neutron") {
      return function evalGroup(group, qMag) {
        return fillConstant(new Float64Array(qMag.length), group.neutronB);
      };
    }
    const coeff = needCoeffSource();
    if (cfg.type === "xray") {
      const xrayModel =
        cfg.model === "table"
          ? { model: "table", tableNum: cfg.tableNum }
          : { model: "waasmaier" };
      return function evalGroup(group, qMag) {
        const out = new Float64Array(qMag.length);
        for (let i = 0; i < qMag.length; i++) {
          out[i] = coeff.evalXray(group.key, qMag[i], xrayModel);
        }
        return out;
      };
    }
    if (cfg.model === "ion_peng") {
      return function evalGroup(group, qMag) {
        const out = new Float64Array(qMag.length);
        for (let i = 0; i < qMag.length; i++) {
          out[i] = coeff.evalElectronIonPeng(group.key, group.valence, qMag[i]);
        }
        return out;
      };
    }
    return function evalGroup(group, qMag) {
      const out = new Float64Array(qMag.length);
      for (let i = 0; i < qMag.length; i++) {
        out[i] = coeff.evalElectronNeutral(group.key, qMag[i], cfg.tableNum);
      }
      return out;
    };
  }

  function validateGroups(groups, evalGroup) {
    const errs = [];
    const qProbe = new Float64Array([0]);
    for (const g of groups) {
      try {
        evalGroup(g, qProbe);
      } catch (e) {
        errs.push(`${g.element}: ${e && e.message ? e.message : String(e)}`);
      }
    }
    if (errs.length) {
      throw new Error(`Unsupported scattering coefficients: ${errs.join(" | ")}`);
    }
  }

  function buildProfile(parsed, scatteringCfg) {
    const core = needCore();
    const cfg = normalizeScatteringConfig(scatteringCfg);
    const profile = {
      config: cfg,
      description: describeConfig(cfg),
      atomCount: parsed.atoms,
      warnings: [],
    };

    const tPrep = nowMs();
    profile.src = core.interleave3(parsed.x, parsed.y, parsed.z);
    profile.srcD = core.interleave3(parsed.dx, parsed.dy, parsed.dz);
    const incSetRaw = parsed.includedElements || null;
    let incSet = null;
    if (incSetRaw && typeof incSetRaw.forEach === "function") {
      incSet = new Set();
      incSetRaw.forEach((v) => {
        const s = String(v || "").trim();
        if (!s) return;
        incSet.add(s);
        incSet.add(normElem(s));
      });
    }
    let includedIndices = [];
    if (!incSet) {
      includedIndices = Array.from({ length: parsed.atoms }, (_, i) => i);
    } else {
      for (let i = 0; i < parsed.atoms; i++) {
        const raw = String(parsed.elements[i] || "").trim();
        const key = normElem(raw);
        if (incSet.has(raw) || incSet.has(key)) includedIndices.push(i);
      }
      if (!includedIndices.length) {
        includedIndices = Array.from({ length: parsed.atoms }, (_, i) => i);
        profile.warnings.push("No included elements resolved; using all atoms.");
      }
    }
    profile.includedAtomCount = includedIndices.length;

    if (cfg.type === "neutron" && cfg.model === "fast") {
      profile.engine = "neutron_fast";
      profile.cAtoms = core.complexReal(parsed.fca);
      profile.prepareMs = nowMs() - tPrep;
      return profile;
    }

    const valenceMap = parseValenceMap(cfg.valenceMapText);
    const rawGroups = groupAtomsByElement(parsed)
      .filter((g) => !incSet || incSet.has(g.element) || incSet.has(g.key));
    const groups = rawGroups.map((g) => {
      const gi = g.indices[0];
      const valence = Number.isInteger(valenceMap[g.key])
        ? valenceMap[g.key]
        : cfg.defaultValence;
      return {
        key: g.key,
        element: g.element,
        count: g.indices.length,
        valence,
        neutronB: Number(parsed.fca[gi]),
        src: packTripletsByIndex(parsed.x, parsed.y, parsed.z, g.indices),
        srcD: packTripletsByIndex(parsed.dx, parsed.dy, parsed.dz, g.indices),
        ones: core.complexOnes(g.indices.length),
      };
    });

    const evalGroup = buildGroupedEvaluator(cfg);
    validateGroups(groups, evalGroup);

    profile.engine = "grouped_exact";
    profile.groups = groups;
    profile.evalGroupCoeffs = evalGroup;
    profile.prepareMs = nowMs() - tPrep;
    return profile;
  }

  // Indices where a per-axis lattice sum (interleaved complex, magnitude up
  // to n) is not zero. On steps commensurate with the supercell only integer
  // h survive, and they form an arithmetic progression.
  function latticeSupport(L, n) {
    const out = [];
    for (let i = 0; i < L.length / 2; i++) {
      if (Math.hypot(L[i * 2], L[i * 2 + 1]) > 1e-9 * n) out.push(i);
    }
    return out;
  }

  function isProgression(values) {
    if (values.length < 3) return true;
    const step = (values[values.length - 1] - values[0]) / (values.length - 1);
    const tol = 1e-9 * Math.max(1, Math.abs(step), Math.abs(values[0]));
    return values.every((v, i) => Math.abs(v - (values[0] + i * step)) <= tol);
  }

  // Copies sub-grid values into a run of the full C-ordered grid; points
  // outside the sub-grid stay zero.
  function gatherSubGrid(values, index, start, count, nk, nl) {
    const out = new Float64Array(count * 2);
    const plane = nk * nl;
    const nkSub = index.nk;
    const nlSub = index.nl;
    for (let q = 0; q < count; q++) {
      const linear = start + q;
      const ih = Math.floor(linear / plane);
      const rem = linear - ih * plane;
      const ik = Math.floor(rem / nl);
      const a = index.h[ih];
      const b = index.k[ik];
      const c = index.l[rem - ik * nl];
      if (a < 0 || b < 0 || c < 0) continue;
      const j = (a * nkSub + b) * nlSub + c;
      out[q * 2] = values[j * 2];
      out[q * 2 + 1] = values[j * 2 + 1];
    }
    return out;
  }

  async function computeIntensity(args) {
    const core = needCore();
    const parsed = args.parsed;
    const h = args.h;
    const k = args.k;
    const l = args.l;
    const Bq = args.Bq;
    const backend = args.backend;
    const opts = args.opts || {};
    const sub = !!args.sub;
    const runType3 = args.runType3;
    const onStatus = args.onStatus;
    const onStageTiming = args.onStageTiming;
    const gridAware = !!args.gridAware;
    // Only for backends that stay accurate on very small grids: the f32
    // type-3 paths are not, so they keep the full-grid A_delta.
    const deltaOnLattice = !!args.deltaOnLattice;
    if (typeof runType3 !== "function") {
      throw new Error("computeIntensity requires runType3 callback");
    }

    const profile = args.profile || buildProfile(parsed, args.scattering);
    const grid = h.length * k.length * l.length;
    const chunkSize = Math.max(
      1,
      Math.floor(Number.isFinite(args.chunkSize) ? args.chunkSize : grid)
    );
    const totalChunks = Math.max(1, Math.ceil(grid / chunkSize));
    const I = new Float64Array(grid);
    const timings = {
      prepare: profile.prepareMs || 0,
      a: 0,
      aavg: 0,
      adelta: 0,
      finalize: 0,
    };
    // Average structure amplitude: L(h) * A_delta(h) / N_cell, where L is the
    // lattice sum over the supercell cells and A_delta uses in-cell offsets.
    // Unlike a transform of the cell origins, this stays exact when cells
    // hold different numbers of atoms (vacancies, element filters).
    let laue = null;
    let cellCount = 1;
    if (sub) {
      if (!parsed.cellIndex || parsed.hasSupercell === false) {
        throw new Error(
          "Subtract avg needs per-atom unit-cell indices (RMC6f, or unified HDF5 with atom_unit_cell). Set Subtract avg OFF for this structure."
        );
      }
      const supercell = parsed.super.map((v) => Math.max(1, Math.round(Number(v) || 1)));
      const c0 = core.cellOrigin(parsed.cellIndex, parsed.atoms, supercell);
      laue = {
        h: core.laueAxis(h, supercell[0], c0[0]),
        k: core.laueAxis(k, supercell[1], c0[1]),
        l: core.laueAxis(l, supercell[2], c0[2]),
      };
      cellCount = supercell[0] * supercell[1] * supercell[2];
    }
    const mInv = 1 / cellCount;
    profile.averageCellCount = sub ? cellCount : null;
    // Largest source count, so grid-aware backends can share one plan.
    const pointCapacity = profile.engine === "neutron_fast"
      ? parsed.atoms
      : Math.max(1, ...profile.groups.map((g) => g.count));

    // L vanishes off integer h on axes whose steps are commensurate with the
    // supercell, so A_delta is only needed on that (small) sub-grid.
    let deltaSub = null;
    let deltaIndex = null;
    if (sub && deltaOnLattice) {
      const axes = [h, k, l];
      const supercell = parsed.super.map((v) => Math.max(1, Math.round(Number(v) || 1)));
      const keep = [laue.h, laue.k, laue.l].map((L, d) => {
        const idx = latticeSupport(L, supercell[d]);
        return isProgression(idx.map((i) => axes[d][i])) ? idx : axes[d].map((_, i) => i);
      });
      if (keep.some((idx, d) => idx.length < axes[d].length)) {
        const [hs, ks, ls] = keep.map((idx, d) => idx.map((i) => axes[d][i]));
        const n = hs.length * ks.length * ls.length;
        deltaIndex = { nk: ks.length, nl: ls.length };
        ["h", "k", "l"].forEach((name, d) => {
          const pos = new Int32Array(axes[d].length).fill(-1);
          keep[d].forEach((i, j) => { pos[i] = j; });
          deltaIndex[name] = pos;
        });
        deltaSub = new Float64Array(n * 2);
        if (n > 0) {
          if (typeof onStatus === "function") onStatus(`Computing Adelta(hkl) on ${n} lattice points ...`);
          const trgSub = gridAware ? null : core.targetsChunk(hs, ks, ls, Bq, 0, n);
          const gridSub = { h: hs, k: ks, l: ls, Bq, start: 0, count: n };
          const transformSub = async (sourcesPacked, strengths) =>
            (
              await runType3(
                { dim: 3, isign: 1, sourcesPacked, targetsPacked: trgSub, strengths, grid: gridSub, pointCapacity },
                opts,
                backend,
                onStageTiming
              )
            ).out;
          const tAd = nowMs();
          if (profile.engine === "neutron_fast") {
            deltaSub.set(await transformSub(profile.srcD, profile.cAtoms));
          } else {
            const qMagSub = core.qMagnitudesChunk(hs, ks, ls, Bq, 0, n);
            for (const g of profile.groups) {
              accumulateScaledComplex(deltaSub, await transformSub(g.srcD, g.ones), profile.evalGroupCoeffs(g, qMagSub));
            }
          }
          timings.adelta += nowMs() - tAd;
        }
      }
    }
    profile.deltaPoints = deltaSub ? deltaSub.length / 2 : null;
    let min = Infinity;
    let max = -Infinity;

    for (let chunk = 0, start = 0; start < grid; chunk++, start += chunkSize) {
      const count = Math.min(chunkSize, grid - start);
      // Grid-aware backends (type-1) take the h,k,l grid instead of targets.
      const trg = gridAware ? null : core.targetsChunk(h, k, l, Bq, start, count);
      const gridSpec = { h, k, l, Bq, start, count };
      const transform = async (sourcesPacked, strengths) =>
        (
          await runType3(
            { dim: 3, isign: 1, sourcesPacked, targetsPacked: trg, strengths, grid: gridSpec, pointCapacity },
            opts,
            backend,
            onStageTiming
          )
        ).out;
      const chunkTag = totalChunks > 1 ? ` (${chunk + 1}/${totalChunks})` : "";
      let q = null;
      let qa = null;
      let qd = null;

      if (sub) {
        const tAa = nowMs();
        qa = core.laueChunk(laue.h, laue.k, laue.l, start, count);
        if (deltaSub) qd = gatherSubGrid(deltaSub, deltaIndex, start, count, k.length, l.length);
        timings.aavg += nowMs() - tAa;
      }
      const fullDelta = sub && !deltaSub;

      if (profile.engine === "neutron_fast") {
        if (typeof onStatus === "function") onStatus(`Computing A(hkl)${chunkTag} ...`);
        const tA = nowMs();
        q = await transform(profile.src, profile.cAtoms);
        timings.a += nowMs() - tA;

        if (fullDelta) {
          if (typeof onStatus === "function") {
            onStatus(`Computing Adelta(hkl)${chunkTag} ...`);
          }
          const tAd = nowMs();
          qd = await transform(profile.srcD, profile.cAtoms);
          timings.adelta += nowMs() - tAd;
        }
      } else {
        const qMag = trg ? makeQMagnitudes(trg) : core.qMagnitudesChunk(h, k, l, Bq, start, count);
        q = new Float64Array(count * 2);
        if (fullDelta) qd = new Float64Array(count * 2);

        for (let gi = 0; gi < profile.groups.length; gi++) {
          const g = profile.groups[gi];
          const gTag =
            profile.groups.length > 1
              ? ` [${gi + 1}/${profile.groups.length} ${g.element}]`
              : "";
          const coeff = profile.evalGroupCoeffs(g, qMag);

          if (typeof onStatus === "function") {
            onStatus(`Computing A(hkl)${chunkTag}${gTag} ...`);
          }
          const tA = nowMs();
          accumulateScaledComplex(q, await transform(g.src, g.ones), coeff);
          timings.a += nowMs() - tA;

          if (fullDelta) {
            if (typeof onStatus === "function") {
              onStatus(`Computing Adelta(hkl)${chunkTag}${gTag} ...`);
            }
            const tAd = nowMs();
            accumulateScaledComplex(qd, await transform(g.srcD, g.ones), coeff);
            timings.adelta += nowMs() - tAd;
          }
        }
      }

      if (typeof onStatus === "function") {
        onStatus(`Finalizing intensity${chunkTag} ...`);
      }
      const tFin = nowMs();
      const mmChunk = core.accumulateIntensityChunk(
        I,
        start,
        count,
        q,
        qa,
        qd,
        mInv,
        sub
      );
      timings.finalize += nowMs() - tFin;
      if (mmChunk.min < min) min = mmChunk.min;
      if (mmChunk.max > max) max = mmChunk.max;
    }

    return {
      I,
      min,
      max,
      timings,
      profile,
      chunkSize,
      totalChunks,
    };
  }

  global.DiffuseAmplitude = Object.freeze({
    normalizeScatteringConfig,
    parseValenceMap,
    buildProfile,
    computeIntensity,
  });
})(typeof window !== "undefined" ? window : globalThis);
