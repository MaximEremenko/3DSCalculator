"""Writes Examples/modulation/Cr_bcc.rmc6f: bcc chromium, 21 x 21 x 21 cells.

The structure for the spin-density-wave example. Chromium orders below 311 K
in an incommensurate SDW with q = (1 - delta, 0, 0), delta about 0.048,
transverse (moments perpendicular to q) above 123 K. With 21 cells along a,
q = 20/21 (delta = 1/21 = 0.0476) is commensurate with the box, so the
satellites at (1 +- delta, 0, 0) fall on grid points. The example button
applies the SDW with the calculator's Modulation tools; the file itself is a
plain bcc lattice (a = 2.884 A).

Usage (from the repository root):
    python tools/make_modulation_examples.py
"""
from pathlib import Path

N, A = 21, 2.884
OUT = Path("Examples/modulation")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sites = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)]
    rows = []
    for cx in range(N):
        for cy in range(N):
            for cz in range(N):
                for s, off in enumerate(sites, start=1):
                    f = [(c + o) / N for c, o in zip((cx, cy, cz), off)]
                    rows.append(f"{len(rows) + 1}  Cr  {f[0]:.8f}  {f[1]:.8f}  {f[2]:.8f}  {s}  {cx}  {cy}  {cz}")
    length = A * N
    header = [
        "(Version 6f format configuration file)",
        f"Metadata title:     bcc Cr, {N} x {N} x {N} cells (structure of the spin-density-wave example)",
        "Atom types present:          Cr",
        f"Number of each atom type:    {len(rows)}",
        f"Number of atoms:             {len(rows)}",
        f"Number density (Ang^-3):     {len(rows) / length ** 3:.6f}",
        f"Supercell dimensions:        {N} {N} {N}",
        f"Cell (Ang/deg):    {length:.6f} {length:.6f} {length:.6f} 90.000000 90.000000 90.000000",
        "Atoms:",
    ]
    path = OUT / "Cr_bcc.rmc6f"
    path.write_text("\n".join(header + rows) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {path}: {len(rows)} atoms")


if __name__ == "__main__":
    main()
