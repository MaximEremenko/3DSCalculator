"""Writes the LAMMPS spin examples in Examples/lammps (atom_style spin, the
current layout: id type x y z spx spy spz sp, units metal).

- fe_bcc_spins.data: bcc Fe, 10 x 10 x 10 cells (2000 atoms, a = 2.8665 A),
  2.2 muB spins tilted about +z with <cos theta> = 0.8 (exp(kappa cos theta)
  distribution) and Gaussian displacements of 0.07 A: magnetic and nuclear
  Bragg peaks plus diffuse scattering.
- skyrmion_film.data: a triangular monolayer of 60 x 60 Fe sites (a = 2.5 A)
  in a triclinic box, with a triple-q (Bloch) skyrmion lattice of period
  10 sites (cores down in a +z background) and 2.5 muB moments: six magnetic satellites around each
  (hk0) Bragg point, and rods along l for the single layer.

hkl of LAMMPS data refer to the whole box: bcc (1 1 0) is (10 10 0) and the
skyrmion wavevectors are at 6 box units.

Usage (from the repository root):
    python tools/make_lammps_examples.py
"""
import math
import random
from pathlib import Path

OUT = Path("Examples/lammps")


def write(path, title, box, tilt, rows):
    lines = [
        f"LAMMPS data file via 3DSCalculator tools/make_lammps_examples.py ({title}), units = metal",
        "",
        f"{len(rows)} atoms",
        "1 atom types",
        "",
        f"0.0 {box[0]:.6f} xlo xhi",
        f"0.0 {box[1]:.6f} ylo yhi",
        f"0.0 {box[2]:.6f} zlo zhi",
    ]
    if tilt:
        lines.append(f"{tilt[0]:.6f} {tilt[1]:.6f} {tilt[2]:.6f} xy xz yz")
    lines += ["", "Masses", "", "1 55.845 # Fe", "", "Atoms # spin", ""]
    for i, (r, s, mag) in enumerate(rows):
        lines.append(f"{i + 1} 1 {r[0]:.6f} {r[1]:.6f} {r[2]:.6f} {s[0]:.8f} {s[1]:.8f} {s[2]:.8f} {mag} 0 0 0")
    (OUT / path).write_text("\n".join(lines) + "\n", encoding="ascii", newline="\n")
    print(f"wrote {OUT / path}: {len(rows)} atoms")


def unit(v):
    n = math.sqrt(sum(x * x for x in v))
    return [x / n for x in v]


def fe_bcc():
    rng = random.Random(20260929)
    a, n, kappa = 2.8665, 10, 5.0
    rows = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for s in ((0, 0, 0), (0.5, 0.5, 0.5)):
                    r = [a * (c + d) + rng.gauss(0, 0.07) for c, d in zip((i, j, k), s)]
                    # cos(theta) from p ~ exp(kappa cos theta) about +z; phi uniform
                    u = rng.random()
                    ct = 1 + math.log(u + (1 - u) * math.exp(-2 * kappa)) / kappa
                    st = math.sqrt(max(0.0, 1 - ct * ct))
                    ph = 2 * math.pi * rng.random()
                    rows.append((r, [st * math.cos(ph), st * math.sin(ph), ct], 2.2))
    write("fe_bcc_spins.data", "bcc Fe, 2.2 muB spins with thermal disorder", [a * n] * 3, None, rows)


def skyrmion_film():
    a, n, period = 2.5, 60, 10
    a1 = (a, 0.0, 0.0)
    a2 = (a / 2, a * math.sqrt(3) / 2, 0.0)
    # reciprocal vectors b_i . a_j = 2 pi delta_ij (in-plane)
    area = a1[0] * a2[1] - a1[1] * a2[0]
    b1 = (2 * math.pi * a2[1] / area, -2 * math.pi * a2[0] / area, 0.0)
    b2 = (-2 * math.pi * a1[1] / area, 2 * math.pi * a1[0] / area, 0.0)
    qs = [[c / period for c in b1], [c / period for c in b2], [-(x + y) / period for x, y in zip(b1, b2)]]
    rows = []
    for i in range(n):
        for j in range(n):
            r = [i * a1[d] + j * a2[d] for d in range(3)]
            # cores pointing down (-z) on a triangular lattice in a +z background
            m = [0.0, 0.0, 0.5]
            for q in qs:
                qhat = unit(q)
                perp = (-qhat[1], qhat[0], 0.0)  # z x qhat: Bloch-type
                ph = q[0] * r[0] + q[1] * r[1]
                for d in range(3):
                    m[d] += -(1.0 if d == 2 else 0.0) * math.cos(ph) + perp[d] * math.sin(ph)
            rows.append((r, unit(m), 2.5))
    box = (n * a, n * a2[1], 5.0)
    write("skyrmion_film.data", "triple-q skyrmion lattice in a triangular Fe monolayer", box, (n * a2[0], 0.0, 0.0), rows)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    fe_bcc()
    skyrmion_film()
