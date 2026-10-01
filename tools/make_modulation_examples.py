"""Writes the structures in Examples/modulation.

KCP_Pt_chains.rmc6f: the Pt chains of KCP, K2Pt(CN)4Br0.3.3H2O, the example
of the Modulations row. Pt atoms 2.89 A apart form chains along c, 9.87 A
apart on a square lattice (17 x 17 chains of 40 Pt). Removing 0.3 electrons
per Pt leaves the Pt band 0.85 full, and each chain has a Peierls wave, a
longitudinal displacement with q = 2 kF = 0.85 per Pt spacing (amplitude
0.03 A). Above about 100 K the chains are not in step: each has its own
phase, so the satellites at l = n +- 0.85 spread into diffuse sheets in h and
k, the diffuse planes found by X-ray diffraction (Comes et al., 1973). The
phases across the chains are a Zadoff-Chu sequence, phi(m, n) = -pi (m (m + 1)
+ n (n + 1)) / 17. Its correlation is zero at every non-zero offset between
chains, so the chains are independent, and the sheets are smooth from this
one configuration; random phases would give the same sheets only on
average, with speckle that many configurations would have to average out.
With 17 chains, a prime, the multiples of the phases (the higher-order
sheets) are such sequences too. The waves are commensurate with the box
(40 x 0.85 = 34), so the sheets lie on grid planes, and on the grid points
of steps that are multiples of 1/17 in h and k.

Cr_bcc.rmc6f: bcc chromium, 21 x 21 x 21 cells, a plain lattice for trying
the Modulation tools. Chromium orders below 311 K in an incommensurate
spin-density wave with q = (1 - delta, 0, 0), delta about 0.048, transverse
(moments perpendicular to q) above 123 K. With 21 cells along a, q = 20/21
(delta = 1/21 = 0.0476) is commensurate with the box, so the satellites at
(1 +- delta, 0, 0) fall on grid points.

Usage (from the repository root):
    python tools/make_modulation_examples.py
"""
import cmath
import math
from pathlib import Path

OUT = Path("Examples/modulation")


def write_rmc6f(name, title, element, supercell, lengths, atoms):
    """atoms: (fractional x, y, z of the box, site, cell indices)."""
    rows = [f"{i}  {element}  {x:.8f}  {y:.8f}  {z:.8f}  {site}  {c[0]}  {c[1]}  {c[2]}"
            for i, (x, y, z, site, c) in enumerate(atoms, start=1)]
    volume = lengths[0] * lengths[1] * lengths[2]
    header = [
        "(Version 6f format configuration file)",
        f"Metadata title:     {title}",
        f"Atom types present:          {element}",
        f"Number of each atom type:    {len(rows)}",
        f"Number of atoms:             {len(rows)}",
        f"Number density (Ang^-3):     {len(rows) / volume:.6f}",
        f"Supercell dimensions:        {supercell[0]} {supercell[1]} {supercell[2]}",
        f"Cell (Ang/deg):    {lengths[0]:.6f} {lengths[1]:.6f} {lengths[2]:.6f} 90.000000 90.000000 90.000000",
        "Atoms:",
    ]
    path = OUT / name
    path.write_text("\n".join(header + rows) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {path}: {len(rows)} atoms")


def kcp_chains():
    chains, length = 17, 40          # chains along a and b (a prime); Pt per chain
    spacing, pt_pt = 9.87, 2.89      # A
    q, amplitude = 0.85, 0.03        # 2 kF per Pt spacing; A
    assert abs(q * length - round(q * length)) < 1e-9
    zc = [-math.pi * m * (m + 1) / chains for m in range(chains)]
    # Equal magnitude at every chain-lattice wavevector, for the first three
    # orders of the wave: no correlation between chains.
    for order in (1, 2, 3):
        for k in range(chains):
            s = sum(cmath.exp(1j * (order * zc[m] - 2 * math.pi * k * m / chains)) for m in range(chains))
            assert abs(abs(s) - math.sqrt(chains)) < 1e-9
    atoms = []
    for cx in range(chains):
        for cy in range(chains):
            phase = zc[cx] + zc[cy]
            for cz in range(length):
                u = amplitude * math.sin(2 * math.pi * q * cz + phase)
                atoms.append((cx / chains, cy / chains, (cz + u / pt_pt) / length, 1, (cx, cy, cz)))
    write_rmc6f("KCP_Pt_chains.rmc6f",
                f"Pt chains of KCP, {chains} x {chains} chains of {length} Pt, each with its own 2kF Peierls wave",
                "Pt", (chains, chains, length), (chains * spacing, chains * spacing, length * pt_pt), atoms)


def cr_bcc():
    n, a = 21, 2.884
    atoms = []
    for cx in range(n):
        for cy in range(n):
            for cz in range(n):
                for site, off in enumerate([(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)], start=1):
                    f = [(c + o) / n for c, o in zip((cx, cy, cz), off)]
                    atoms.append((f[0], f[1], f[2], site, (cx, cy, cz)))
    write_rmc6f("Cr_bcc.rmc6f", f"bcc Cr, {n} x {n} x {n} cells (structure of the spin-density-wave example)",
                "Cr", (n, n, n), (a * n, a * n, a * n), atoms)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    kcp_chains()
    cr_bcc()


if __name__ == "__main__":
    main()
