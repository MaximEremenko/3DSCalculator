"use strict";

// Modulated structures, commensurate or incommensurate: a wave with
// wavevector q (reciprocal-lattice units of the parent cell; any real
// values) applied to a loaded structure. The phase at atom j is
//   theta_j = 2 pi q . r_j + phi,
// with r_j its unmodulated position in parent-cell fractions. Kinds:
//   displacive  u_j = A e sin(theta_j)                  (A in Angstrom)
//   sdw         mu_j = A e cos(theta_j)                 (spin-density wave, A in muB)
//   helix       mu_j = A (e1 cos(theta_j) + e2 sin(theta_j)), e1, e2 _|_ axis
// The scattering is computed for the finite box as it is, so no periodic
// approximant is needed: satellites at G +- q (and G +- m q for large
// displacements) come out with the box's finite-size width.
(function (global) {
  function cart(v, direct) {
    return [0, 1, 2].map((a) => v[0] * direct[0][a] + v[1] * direct[1][a] + v[2] * direct[2][a]);
  }

  function unit(v, what) {
    const n = Math.hypot(v[0], v[1], v[2]);
    if (!(n > 1e-12)) throw new Error(`Modulation: the ${what} must be a non-zero vector`);
    return v.map((x) => x / n);
  }

  function inverse3(m) {
    const [a, b, c] = m[0], [d, e, f] = m[1], [g, h, i] = m[2];
    const A = e * i - f * h, B = -(d * i - f * g), C = d * h - e * g;
    const det = a * A + b * B + c * C;
    const s = 1 / det;
    return [
      [A * s, -(b * i - c * h) * s, (b * f - c * e) * s],
      [B * s, (a * i - c * g) * s, -(a * f - c * d) * s],
      [C * s, -(a * h - b * g) * s, (a * e - b * d) * s],
    ];
  }

  // Parent-cell direct basis (rows a, b, c in Angstrom) from Bp = (direct^-1)^T.
  function directBasis(parsed) {
    const Bp = parsed.Bp;
    const t = [0, 1, 2].map((i) => [0, 1, 2].map((j) => Bp[j][i]));
    return inverse3(t);
  }

  // Whether q is commensurate with the box along each axis (N q integer).
  function commensurability(q, supercell) {
    return q.map((v, d) => {
      const n = Math.max(1, Math.round(Number(supercell && supercell[d]) || 1));
      const x = v * n;
      return Math.abs(x - Math.round(x)) < 1e-6;
    });
  }

  function applyModulation(parsed, spec) {
    const kind = String(spec.kind || "").toLowerCase();
    if (!["displacive", "sdw", "helix"].includes(kind)) throw new Error(`Modulation: unknown kind "${spec.kind}"`);
    const q = (spec.q || []).map(Number);
    if (q.length !== 3 || !q.every(Number.isFinite)) throw new Error("Modulation: q needs three numbers");
    const amplitude = Number(spec.amplitude);
    if (!Number.isFinite(amplitude)) throw new Error("Modulation: invalid amplitude");
    const phase = ((Number(spec.phaseDeg) || 0) * Math.PI) / 180;
    const direct = directBasis(parsed);
    const dir = unit(cart((spec.direction || [1, 0, 0]).map(Number), direct), "direction");
    const wanted = spec.elements && spec.elements.length ? new Set(spec.elements.map((e) => String(e).toLowerCase())) : null;
    const n = parsed.atoms;
    const selected = new Uint8Array(n);
    let count = 0;
    for (let j = 0; j < n; j++) {
      if (!wanted || wanted.has(String(parsed.elements[j]).toLowerCase())) { selected[j] = 1; count++; }
    }
    if (!count) throw new Error("Modulation: no atoms of the chosen elements");
    const Bp = parsed.Bp;
    const theta = new Float64Array(n);
    for (let j = 0; j < n; j++) {
      const r = [parsed.x[j], parsed.y[j], parsed.z[j]];
      let t = 0;
      for (let i = 0; i < 3; i++) t += q[i] * (r[0] * Bp[i][0] + r[1] * Bp[i][1] + r[2] * Bp[i][2]);
      theta[j] = 2 * Math.PI * t + phase;
    }
    const out = { ...parsed };
    const info = { kind, q, amplitude, phaseDeg: Number(spec.phaseDeg) || 0, atoms: count, commensurate: commensurability(q, parsed.super) };

    if (kind === "displacive") {
      for (const key of ["x", "y", "z", "dx", "dy", "dz"]) out[key] = Float64Array.from(parsed[key]);
      for (let j = 0; j < n; j++) {
        if (!selected[j]) continue;
        const s = amplitude * Math.sin(theta[j]);
        out.x[j] += s * dir[0]; out.y[j] += s * dir[1]; out.z[j] += s * dir[2];
        out.dx[j] += s * dir[0]; out.dy[j] += s * dir[1]; out.dz[j] += s * dir[2];
      }
      info.direction = dir;
    } else {
      let e1 = dir, e2 = null;
      if (kind === "helix") {
        // rotation plane _|_ axis; e1 is the axis-perpendicular part of a
        // reference direction (the lattice vector least aligned with the axis)
        const axis = dir;
        const refs = direct.map((row) => unit(row, "cell vector"));
        const ref = refs.reduce((best, r) => (Math.abs(r[0] * axis[0] + r[1] * axis[1] + r[2] * axis[2]) < Math.abs(best[0] * axis[0] + best[1] * axis[1] + best[2] * axis[2]) ? r : best));
        const d = ref[0] * axis[0] + ref[1] * axis[1] + ref[2] * axis[2];
        e1 = unit([ref[0] - d * axis[0], ref[1] - d * axis[1], ref[2] - d * axis[2]], "helix plane");
        e2 = [axis[1] * e1[2] - axis[2] * e1[1], axis[2] * e1[0] - axis[0] * e1[2], axis[0] * e1[1] - axis[1] * e1[0]];
        info.axis = axis;
      } else {
        info.direction = dir;
      }
      const ff = global.MagneticFormFactors;
      const species = [];
      const index = new Map();
      const speciesOfAtom = new Int32Array(n).fill(-1);
      const vectors = new Float64Array(n * 3);
      for (let j = 0; j < n; j++) {
        if (!selected[j]) continue;
        const el = String(parsed.elements[j]);
        if (!index.has(el)) {
          const ion = ff ? ((/\d/.test(el) && ff.find(el)) || ff.defaultIon(el)) : null;
          index.set(el, species.length);
          species.push({
            label: el, element: el, ion: ion ? ion.ion : null, scale: 1,
            formFactor: ion ? { j0: ion.j0, j2: ion.j2, c2: 0 } : { j0: [0, 0, 0, 0, 0, 0, 0, 0, 1], j2: null, c2: 0 },
            formFactorSource: ion ? "table" : "none",
          });
        }
        speciesOfAtom[j] = index.get(el);
        const c = Math.cos(theta[j]), s = Math.sin(theta[j]);
        for (let a = 0; a < 3; a++) vectors[j * 3 + a] = amplitude * (kind === "sdw" ? e1[a] * c : e1[a] * c + e2[a] * s);
      }
      out.magnetic = {
        species, speciesOfAtom, vectors,
        source: kind === "sdw" ? "spin-density wave (modulation)" : "helix (modulation)",
        notes: [],
      };
    }
    out.modulation = info;
    out.sourceFormat = `${parsed.sourceFormat || "structure"} + ${kind} modulation`;
    return out;
  }

  global.Modulation = Object.freeze({ applyModulation, commensurability });
})(typeof window !== "undefined" ? window : globalThis);
