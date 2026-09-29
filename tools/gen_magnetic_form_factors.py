"""Regenerates js/magnetic_form_factors.js from Sunny.jl's src/FormFactor.jl.

Usage (from the repository root):
    python tools/gen_magnetic_form_factors.py path/to/Sunny.jl/src/FormFactor.jl
"""
import re
import sys
from pathlib import Path

ROW = re.compile(
    r'^\s*\("([A-Za-z]+[0-9]+)",\s*"([^"]*)",\s*"([^"]*)",\s*\[([^\]]*)\],\s*\[([^\]]*)\]\),?\s*(#.*)?$'
)

TEMPLATE = r'''"use strict";

// Magnetic form factors in the dipole approximation:
//   f(Q) = <j0>(s) + C2 * <j2>(s),   s = Q / (4 pi),   Q in 1/Angstrom,
//   <j0>(s) = A e^{-a s^2} + B e^{-b s^2} + C e^{-c s^2} + D e^{-d s^2} + E,
//   <j2>(s) = s^2 * (the same form with the <j2> coefficients).
// C2 = (2 - gJ)/gJ for a J multiplet (4f, 5f) and about 0 for spin-only
// 3d ions. Coefficients: P. J. Brown (3d, 4d, 4f, 5f), Lisher and Forsyth
// (Ce3+), Kobayashi, Nagao and Ito (5d), as tabulated in Sunny.jl
// src/FormFactor.jl (MIT, see js/magnetic_form_factors-LICENSE.txt).
// Generated from that file; regenerate rather than edit the table.
(function (global) {
  // [ion (element + charge), electron configuration, ground term,
  //  <j0> A a B b C c D d E, <j2> A a B b C c D d E]
  const IONS = [
%TABLE%
  ];

  const ORBITALS = "SPDFGHIKLMNOQRTUV";

  // "Mn2+", "Mn+2", "mn2" -> "Mn2"; a bare element is the neutral atom.
  function ionKey(label) {
    const m = /^\s*([A-Za-z]{1,2})\s*\+?(\d*)\s*\+?\s*$/.exec(String(label || ""));
    if (!m) return null;
    return m[1][0].toUpperCase() + m[1].slice(1).toLowerCase() + (m[2] ? String(Number(m[2])) : "0");
  }

  function entry(row) {
    return { ion: row[0], config: row[1], term: row[2], j0: row[3], j2: row[4] };
  }

  // First tabulated ion for a label (optionally a specific configuration).
  function find(label, config) {
    const key = ionKey(label);
    const row = IONS.find((r) => r[0] === key && (!config || r[1] === config));
    return row ? entry(row) : null;
  }

  function ions() {
    return IONS.map(entry);
  }

  // Lande g of a Russell-Saunders term such as "6H15/2".
  function landeG(term) {
    const m = /^(\d+)([A-Z])(\d+)(?:\/(\d+))?$/.exec(String(term || ""));
    if (!m) return null;
    const S = (Number(m[1]) - 1) / 2;
    const L = ORBITALS.indexOf(m[2]);
    const J = m[4] ? Number(m[3]) / Number(m[4]) : Number(m[3]);
    if (L < 0 || !(J > 0)) return null;
    return 1 + (J * (J + 1) + S * (S + 1) - L * (L + 1)) / (2 * J * (J + 1));
  }

  // (2 - gJ)/gJ for 4f and 5f ions from their Hund's-rule term, else 0.
  function suggestedC2(ion) {
    if (!ion || !/[45]f/.test(ion.config)) return 0;
    const g = landeG(ion.term);
    return g ? (2 - g) / g : 0;
  }

  // 9 coefficients from 7 (A a B b C c D) or 9 (A a B b C c D d E) values.
  function coefficients(values) {
    const v = Array.from(values, Number);
    if (v.length === 7) return [v[0], v[1], v[2], v[3], v[4], v[5], 0, 0, v[6]];
    if (v.length === 9) return v;
    throw new Error("A magnetic form factor needs 7 or 9 coefficients.");
  }

  function gaussian(c, s2) {
    return c[0] * Math.exp(-c[1] * s2) + c[2] * Math.exp(-c[3] * s2) +
      c[4] * Math.exp(-c[5] * s2) + c[6] * Math.exp(-c[7] * s2) + c[8];
  }

  // f(Q) for a form factor { j0, j2, c2 } (9-coefficient arrays).
  function evaluate(ff, q) {
    const s2 = (q / (4 * Math.PI)) ** 2;
    const j0 = gaussian(ff.j0, s2);
    return ff.c2 && ff.j2 ? j0 + ff.c2 * s2 * gaussian(ff.j2, s2) : j0;
  }

  // Tabulated ion whose <j0> matches 7 or 9 given coefficients, if any.
  function identify(j0) {
    const c = coefficients(j0);
    const row = IONS.find((r) => r[3].every((v, i) => Math.abs(v - c[i]) <= 5e-4 * Math.max(1, Math.abs(v))));
    return row ? entry(row) : null;
  }

  global.MagneticFormFactors = Object.freeze({
    ionKey,
    find,
    ions,
    landeG,
    suggestedC2,
    coefficients,
    evaluate,
    identify,
  });
})(typeof window !== "undefined" ? window : globalThis);
'''


def number(v):
    s = repr(v)
    return s[:-2] if s.endswith(".0") else s


def main(source):
    rows = []
    for line in Path(source).read_text(encoding="utf-8").splitlines():
        m = ROW.match(line)
        if not m:
            continue
        j0 = [float(v) for v in m.group(4).split(",")]
        j2 = [float(v) for v in m.group(5).split(",")]
        if len(j0) != 9 or len(j2) != 9:
            raise SystemExit(f"unexpected coefficient count for {m.group(1)}")
        rows.append(
            f'    ["{m.group(1)}", "{m.group(2)}", "{m.group(3)}", '
            f'[{", ".join(number(v) for v in j0)}], [{", ".join(number(v) for v in j2)}]],'
        )
    Path("js/magnetic_form_factors.js").write_text(
        TEMPLATE.replace("%TABLE%", "\n".join(rows)), encoding="utf-8", newline="\n"
    )
    print(f"{len(rows)} ions")


if __name__ == "__main__":
    main(sys.argv[1])
