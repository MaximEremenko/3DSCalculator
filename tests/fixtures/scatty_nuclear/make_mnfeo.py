"""Rocksalt (Mn,Fe)O with random cation mixing and displacements, as a Scatty atoms file.

Writes mnfeo_atoms_01.txt and two scatty_config variants (total, and REMOVE_BRAGG F).
Mn and Fe are shifted by +-0.012 lattice units along x, so occupancy and displacement
are correlated; every atom also gets random displacements of 0.01 lattice units.
"""
import random
from pathlib import Path

random.seed(20260929)
A, BOX = 4.4344, 5
CATIONS = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
ANIONS = [(0.5, 0, 0), (0, 0.5, 0), (0, 0, 0.5), (0.5, 0.5, 0.5)]
sites = CATIONS + ANIONS
lines = ["TITLE mnfeo", f"CELL {A} {A} {A} 90 90 90"]
lines += [f"SITE {x} {y} {z}" for x, y, z in sites]
atoms = []
for cx in range(BOX):
    for cy in range(BOX):
        for cz in range(BOX):
            for s in range(8):
                el = ("Mn" if random.random() < 0.5 else "Fe") if s < 4 else "O"
                # Mn and Fe relax in opposite directions along x, plus random noise
                shift = (0.012 if el == "Mn" else -0.012 if el == "Fe" else 0.0)
                u = [random.gauss(0, 0.01) + (shift if d == 0 else 0) for d in range(3)]
                atoms.append(f"ATOM {s + 1} {cx} {cy} {cz} {u[0]:.6f} {u[1]:.6f} {u[2]:.6f} {el}")
# OCC: the supercell's overall cation composition on all four cation sites, so
# that Scatty's average structure is F-centred like the true average.
cat = [a.split()[-1] for a in atoms if int(a.split()[1]) <= 4]
mn = cat.count("Mn") / len(cat)
lines += [f"OCC Mn {mn:.10f} Fe {1 - mn:.10f}"] * 4 + ["OCC O 1.0000000000"] * 4
lines += [f"BOX {BOX} {BOX} {BOX}"] + atoms
Path("mnfeo_atoms_01.txt").write_text("\n".join(lines) + "\n")
base = ["CENTRE 0 0 0", "X_AXIS 2 0 0 10", "Y_AXIS 0 2 0 10", "Z_AXIS 0 0 2 10", "WINDOW 0", "RADIATION N", "SYMMETRY -1", "EXPANSION_MAX_ERROR 1e-12"]
Path("config_total.txt").write_text("\n".join(["NAME total"] + base) + "\n")
Path("config_diffuse.txt").write_text("\n".join(["NAME diffuse"] + base + ["REMOVE_BRAGG F"]) + "\n")
print("atoms", BOX ** 3 * 8)
