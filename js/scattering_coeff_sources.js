"use strict";

(function (global) {
  function needDb() {
    if (!global.ScatteringCoeffDB) {
      throw new Error("ScatteringCoeffDB missing. Load ./js/scattering_coeff_db.js first.");
    }
    return global.ScatteringCoeffDB;
  }

  function normElem(sym) {
    return String(sym || "").trim().toLowerCase();
  }

  function resolveZ(symOrZ) {
    const db = needDb();
    if (Number.isInteger(symOrZ)) return symOrZ;
    const key = normElem(symOrZ);
    const z = db.elementToZ[key];
    if (!Number.isInteger(z)) {
      throw new Error(`Unsupported element symbol: "${symOrZ}"`);
    }
    return z;
  }

  function getElementByZ(z) {
    const db = needDb();
    return db.zToElement[String(z)] || "";
  }

  function gAndG2FromQ(qTemp) {
    const db = needDb();
    const g = Number(qTemp) / db.constants.pi4;
    const gsq = g * g;
    return { g, gsq };
  }

  // Kirkland and Lobato parameterize in q = 1/d = 2 sin(theta)/lambda, the
  // other tables in s = sin(theta)/lambda; qTemp is |Q| = 4 pi s.
  function qSqFromQ(qTemp) {
    const q = Number(qTemp) / (2 * Math.PI);
    return q * q;
  }

  // Weickenmeier & Kohl (1991) store B1..B6 and V per element (Z = 2-98);
  // the amplitudes follow from Z: A1..A3 = Z/(3(1+V)), A4..A6 = V Z/(3(1+V)),
  // so that sum A = Z. The table keeps V in each of the six "a" slots.
  function weickenmeierTerms(z, rec) {
    const b = rec.b || [];
    const v = Number((rec.a || [])[0]);
    if (!(v > 0) || !b.some((x) => Number(x) > 0)) {
      throw new Error(`Weickenmeier-Kohl has no parameters for Z=${z}; choose another table.`);
    }
    return b.slice(0, 6).map((bi, i) => ({ A: ((i < 3 ? 1 : v) * z) / (3 * (1 + v)), B: Number(bi) }));
  }

  function sumExpTerms(a, b, gsq, nTerms) {
    const n = Math.max(
      0,
      Math.min(
        Number.isFinite(nTerms) ? Math.floor(nTerms) : Math.max(a.length, b.length),
        Math.max(a.length, b.length)
      )
    );
    let s = 0;
    for (let i = 0; i < n; i++) {
      s += Number(a[i] || 0) * Math.exp(-Number(b[i] || 0) * gsq);
    }
    return s;
  }

  function evalExpSeriesModel(rec, gsq, model) {
    const cfg = model || {};
    const terms = Number.isFinite(cfg.terms)
      ? Math.floor(cfg.terms)
      : Math.max((rec.a || []).length, (rec.b || []).length);
    const sExp = sumExpTerms(rec.a || [], rec.b || [], gsq, terms);
    const factor = Number.isFinite(cfg.factor) ? Number(cfg.factor) : 1;
    const usesG2 = !!cfg.usesG2;
    const sign = Number.isFinite(cfg.sign) ? Number(cfg.sign) : 1;
    const constant = Number.isFinite(cfg.constant) ? Number(cfg.constant) : 0;
    const scaled = sExp * factor * (usesG2 ? gsq : 1);
    return constant + sign * scaled;
  }

  function evalXrayWaasmaier(symOrZ, qTemp) {
    const db = needDb();
    const z = resolveZ(symOrZ);
    const rec = db.xrayWk95ByZ[String(z)];
    if (!rec) throw new Error(`Waasmaier constants missing for Z=${z}`);
    const { gsq } = gAndG2FromQ(qTemp);
    return evalExpSeriesModel(rec, gsq, {
      terms: 5,
      factor: 1,
      usesG2: false,
      sign: 1,
      constant: Number(rec.c || 0),
    });
  }

  function evalXrayWk95(symOrZ, qTemp) {
    // Backward-compatible alias.
    return evalXrayWaasmaier(symOrZ, qTemp);
  }

  function tableNameByNum(tableNum) {
    switch (Number(tableNum) || 0) {
      case 1:
        return "lobato";
      case 2:
        return "peng_0_4";
      case 3:
        return "doyle";
      case 4:
        return "weickenmeier";
      case 5:
        return "kirkland";
      default:
        return "lobato";
    }
  }

  function xrayTableNameByNum(tableNum) {
    switch (Number(tableNum) || 0) {
      case 0:
        return "waasmaier";
      case 1:
        return "lobato";
      case 2:
        return "peng_0_4";
      case 3:
        return "doyle";
      case 4:
        return "weickenmeier";
      case 5:
        return "kirkland";
      default:
        return "waasmaier";
    }
  }

  function getNeutralTableRec(symOrZ, tableNum) {
    const db = needDb();
    const z = resolveZ(symOrZ);
    const table = tableNameByNum(tableNum);
    const byZ = db.scTables[table] || {};
    const rec = byZ[String(z)];
    if (!rec) throw new Error(`Scattering table "${table}" missing for Z=${z}`);
    return { z, rec, table };
  }

  function getXrayTableRec(symOrZ, tableNum) {
    const db = needDb();
    const z = resolveZ(symOrZ);
    const table = xrayTableNameByNum(tableNum);
    if (table === "waasmaier") {
      const rec = db.xrayWk95ByZ[String(z)];
      if (!rec) {
        throw new Error(`Waasmaier constants missing for Z=${z}`);
      }
      return { z, rec, table };
    }
    const byZ = db.scTables[table] || {};
    const rec = byZ[String(z)];
    if (!rec) throw new Error(`Scattering table "${table}" missing for Z=${z}`);
    return { z, rec, table };
  }

  function evalElectronNeutral(symOrZ, qTemp, tableNum) {
    const { z, rec, table } = getNeutralTableRec(symOrZ, tableNum);
    const { gsq } = gAndG2FromQ(qTemp);
    const a = rec.a || [];
    const b = rec.b || [];
    let f = 0;
    if (table === "lobato") {
      const qsq = qSqFromQ(qTemp);
      for (let i = 0; i < 5; i++) {
        const ai = Number(a[i] || 0);
        const bi = Number(b[i] || 0);
        const d = 1 + bi * qsq;
        f += ai * (2 + bi * qsq) / (d * d);
      }
      return f;
    }
    if (table === "peng_0_4") {
      return evalExpSeriesModel(rec, gsq, {
        terms: 5,
        factor: 1,
        usesG2: false,
        sign: 1,
        constant: 0,
      });
    }
    if (table === "doyle") {
      return evalExpSeriesModel(rec, gsq, {
        terms: 4,
        factor: 1,
        usesG2: false,
        sign: 1,
        constant: 0,
      });
    }
    if (table === "weickenmeier") {
      // Mott-Bethe of the X-ray form: f_e = sum A (1 - exp(-B s^2)) / (8 pi^2 a0 s^2)
      const eightPi2a0 = 4 * Number(needDb().constants.sqPi2a0);
      for (const t of weickenmeierTerms(z, rec)) {
        f += gsq > 0 ? (t.A * (1 - Math.exp(-t.B * gsq))) / gsq : t.A * t.B;
      }
      return f / eightPi2a0;
    }
    if (table === "kirkland") {
      const qsq = qSqFromQ(qTemp);
      for (let i = 0; i < 3; i++) {
        const ai = Number(a[i] || 0);
        const bi = Number(b[i] || 0);
        const ci = Number(a[3 + i] || 0);
        const di = Number(b[3 + i] || 0);
        if (bi + qsq !== 0) {
          f += ai / (bi + qsq) + ci * Math.exp(-di * qsq);
        } else {
          f += ci * Math.exp(-di * qsq);
        }
      }
      return f;
    }
    return f;
  }

  function evalXrayFromScTable(symOrZ, qTemp, tableNum) {
    const db = needDb();
    const { z, rec, table } = getXrayTableRec(symOrZ, tableNum);
    const { gsq } = gAndG2FromQ(qTemp);
    const a = rec.a || [];
    const b = rec.b || [];
    const sqPi2a0 = Number(db.constants.sqPi2a0);
    let f = 0;
    if (table === "waasmaier") {
      return evalExpSeriesModel(rec, gsq, {
        terms: 5,
        factor: 1,
        usesG2: false,
        sign: 1,
        constant: Number(rec.c || 0),
      });
    }
    if (table === "lobato") {
      const qsq = qSqFromQ(qTemp);
      for (let i = 0; i < 5; i++) {
        const ai = Number(a[i] || 0);
        const bi = Number(b[i] || 0);
        if (bi !== 0) f += (sqPi2a0 * ai) / (bi * (1 + bi * qsq) ** 2);
      }
      return f;
    }
    if (table === "peng_0_4") {
      // Mott-Bethe from the electron fit in s: f_x = Z - 8 pi^2 a0 s^2 f_e(s)
      return evalExpSeriesModel(rec, gsq, {
        terms: 5,
        factor: 4 * sqPi2a0,
        usesG2: true,
        sign: -1,
        constant: z,
      });
    }
    if (table === "doyle") {
      // Mott-Bethe from the electron fit in s: f_x = Z - 8 pi^2 a0 s^2 f_e(s)
      return evalExpSeriesModel(rec, gsq, {
        terms: 4,
        factor: 4 * sqPi2a0,
        usesG2: true,
        sign: -1,
        constant: z,
      });
    }
    if (table === "weickenmeier") {
      return weickenmeierTerms(z, rec).reduce((sum, t) => sum + t.A * Math.exp(-t.B * gsq), 0);
    }
    if (table === "kirkland") {
      const qsq = qSqFromQ(qTemp);
      f = z;
      for (let i = 0; i < 3; i++) {
        const ai = Number(a[i] || 0);
        const bi = Number(b[i] || 0);
        const ci = Number(a[3 + i] || 0);
        const di = Number(b[3 + i] || 0);
        if (bi + qsq !== 0) {
          f -= sqPi2a0 * qsq * ai / (bi + qsq);
          f -= sqPi2a0 * qsq * ci * Math.exp(-di * qsq);
        } else {
          f -= sqPi2a0 * qsq * ci * Math.exp(-di * qsq);
        }
      }
      return f;
    }
    return f;
  }

  function normalizeXrayModel(modelOrOptions) {
    if (modelOrOptions == null) return { model: "waasmaier", tableNum: 0 };
    if (Number.isFinite(Number(modelOrOptions))) {
      const n = Number(modelOrOptions);
      if (n === 0) {
        return { model: "waasmaier", tableNum: 0 };
      }
      return {
        model: "table",
        tableNum: Math.max(0, Math.min(5, Math.floor(Number(modelOrOptions)))),
      };
    }
    if (typeof modelOrOptions === "string") {
      const m = String(modelOrOptions).trim().toLowerCase();
      if (m === "waasmaier" || m === "wk95") {
        return { model: "waasmaier", tableNum: 0 };
      }
      if (m === "table") return { model: "table", tableNum: 0 };
      if (m === "peng_0_4") return { model: "table", tableNum: 2 };
      if (m === "doyle") return { model: "table", tableNum: 3 };
      if (m === "weickenmeier") return { model: "table", tableNum: 4 };
      if (m === "lobato") return { model: "table", tableNum: 1 };
      if (m === "kirkland") return { model: "table", tableNum: 5 };
      throw new Error(`Unsupported xray model: "${modelOrOptions}"`);
    }
    const rawModel = String(
      modelOrOptions.model || modelOrOptions.kind || "waasmaier"
    )
      .trim()
      .toLowerCase();
    const tableNum = Number.isFinite(Number(modelOrOptions.tableNum))
      ? Math.max(0, Math.min(5, Math.floor(Number(modelOrOptions.tableNum))))
      : 0;
    if (rawModel === "waasmaier" || rawModel === "wk95") {
      return { model: "waasmaier", tableNum: 0 };
    }
    if (rawModel === "table") return { model: "table", tableNum };
    if (rawModel === "peng_0_4") return { model: "table", tableNum: 2 };
    if (rawModel === "doyle") return { model: "table", tableNum: 3 };
    if (rawModel === "weickenmeier") return { model: "table", tableNum: 4 };
    if (rawModel === "lobato") return { model: "table", tableNum: 1 };
    if (rawModel === "kirkland") return { model: "table", tableNum: 5 };
    throw new Error(`Unsupported xray model: "${rawModel}"`);
  }

  function evalXray(symOrZ, qTemp, modelOrOptions) {
    const x = normalizeXrayModel(modelOrOptions);
    if (x.model === "waasmaier") return evalXrayWaasmaier(symOrZ, qTemp);
    return evalXrayFromScTable(symOrZ, qTemp, x.tableNum);
  }

  function evalElectronIonPeng(symOrZ, valence, qTemp) {
    const db = needDb();
    const z = resolveZ(symOrZ);
    const key = `${z}:${Number(valence)}`;
    const rec = db.pengIonByZValence[key];
    if (!rec) throw new Error(`Peng ion coefficients missing for key ${key}`);
    const { gsq } = gAndG2FromQ(qTemp);
    const a = rec.a || [];
    const b = rec.b || [];
    let f = evalExpSeriesModel(rec, gsq, {
      terms: 5,
      factor: 1,
      usesG2: false,
      sign: 1,
      constant: 0,
    });
    if (gsq !== 0) {
      f += Number(db.constants.ionCoulombFactor) * Number(valence) / gsq;
    }
    return f;
  }

  function evalArray(evalFn, symOrZ, qArray, opt) {
    const out = new Float64Array(qArray.length);
    for (let i = 0; i < qArray.length; i++) {
      out[i] = evalFn(symOrZ, qArray[i], opt);
    }
    return out;
  }

  global.ScatteringCoeffSource = Object.freeze({
    resolveZ,
    getElementByZ,
    tableNameByNum,
    xrayTableNameByNum,
    evalXray,
    evalXrayWaasmaier,
    evalXrayWk95,
    evalXrayFromScTable,
    evalElectronNeutral,
    evalElectronIonPeng,
    evalArray,
  });
})(typeof window !== "undefined" ? window : globalThis);
