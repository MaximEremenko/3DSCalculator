"""Writes the MnO spin configurations used by tests/magnetic.test.js.

Scatty spin files of the fcc Mn sublattice of MnO (a = 4.4344 Angstrom and
the Mn2+ <j0> form factor of Spinteract's examples/MnO) in a 6x6x6 supercell
(864 spins of 5 muB):

  mno_random_spins_01.txt  uniformly random directions (ideal paramagnet)
  mno_order_spins_01.txt   type-II antiferromagnet, S = 5 [1 1 -2]/sqrt(6) (-1)^(x+y+z)

Usage: python make_mno_spins.py   (needs numpy; writes into the current directory)
"""
import numpy as np

A = 4.4344
BOX = 6
SITES = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
J0 = "0.4220 17.6840 0.5948 6.0050 0.0043 -0.6090 -0.0219"
MU = 5.0


def write(stem, spin_of):
    lines = [f"TITLE {stem}", f"CELL {A} {A} {A} 90 90 90"]
    lines += [f"SITE {x} {y} {z}" for x, y, z in SITES]
    lines += [f"BOX {BOX} {BOX} {BOX}", f"FORM_FACTOR_J0 {J0}"]
    for rx in range(BOX):
        for ry in range(BOX):
            for rz in range(BOX):
                for s, (x, y, z) in enumerate(SITES):
                    m = spin_of(rx + x, ry + y, rz + z)
                    lines.append(f"SPIN {s + 1} {rx} {ry} {rz} {m[0]:.12f} {m[1]:.12f} {m[2]:.12f}")
    with open(f"{stem}_spins_01.txt", "w") as f:
        f.write("\n".join(lines) + "\n")


rng = np.random.default_rng(20260929)


def random_spin(*_):
    v = rng.normal(size=3)
    return MU * v / np.linalg.norm(v)


AXIS = np.array([1.0, 1.0, -2.0]) / np.sqrt(6.0)


def ordered_spin(x, y, z):
    return MU * AXIS * (-1) ** int(round(x + y + z))


write("mno_random", random_spin)
write("mno_order", ordered_spin)
