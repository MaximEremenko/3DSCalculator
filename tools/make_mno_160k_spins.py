"""Writes Examples/magnetic/MnO_160K_spins_NN.txt: paramagnetic MnO at 160 K.

Classical Monte Carlo of the J1-J2 Heisenberg model on the fcc Mn lattice with
the couplings Paddison refined to single-crystal diffuse scattering at 160 K
(Spinteract, arXiv:2210.09016, Table I): J1 = 3.26 K, J2 = 4.45 K, with
antiferromagnetic couplings positive, spin pairs double counted and spins of
length sqrt(S(S+1)), S = 5/2. For unit vectors s the energy is the sum over
unordered pairs of K s_i.s_j with K = 2 J S(S+1).

A 16 x 16 x 16 box of conventional cells (16384 spins) is equilibrated with
heat-bath and over-relaxation sweeps over 32 sublattices that hold no nearest
or next-nearest pair. Snapshots follow every 20 sweeps. Moments are
g sqrt(S(S+1)) = 5.916 muB along the spins.

At 160 K the model is paramagnetic (T_N of MnO is 118 K) with short-range
type-II antiferromagnetic correlations: <s.s> about -0.05 for nearest and -0.15
for next-nearest neighbours. Averaged over snapshots, the calculator's
magnetic scattering matches Paddison's SXD data at 160 K with R_wp 9-10 %
(Spinteract's own fit: 8.76 %).

Usage (from the repository root):
    python tools/make_mno_160k_spins.py
"""
from pathlib import Path

import numpy as np

L, SNAPSHOTS, T = 16, 8, 160.0
J1, J2, S2 = 3.26, 4.45, 8.75
K1, K2 = 2 * J1 * S2, 2 * J2 * S2
MU = 2 * np.sqrt(S2)
A = 4.4344
OUT = Path("Examples/magnetic")
J0 = "0.4220 17.6840 0.5948 6.0050 0.0043 -0.6090 -0.0219"

n = 2 * L  # positions on the half-cell grid: x + y + z even
coords = np.array([(x, y, z) for x in range(n) for y in range(n) for z in range(n) if (x + y + z) % 2 == 0])
N = len(coords)
grid = np.full((n, n, n), -1, dtype=np.int64)
grid[coords[:, 0], coords[:, 1], coords[:, 2]] = np.arange(N)


def neighbours(vectors):
    out = np.empty((N, len(vectors)), dtype=np.int64)
    for j, v in enumerate(vectors):
        c = (coords + np.array(v)) % n
        out[:, j] = grid[c[:, 0], c[:, 1], c[:, 2]]
    return out


NN = neighbours([(a, b, 0) for a in (1, -1) for b in (1, -1)] + [(a, 0, b) for a in (1, -1) for b in (1, -1)] + [(0, a, b) for a in (1, -1) for b in (1, -1)])
NNN = neighbours([(2, 0, 0), (-2, 0, 0), (0, 2, 0), (0, -2, 0), (0, 0, 2), (0, 0, -2)])
key = (coords[:, 0] % 4) * 16 + (coords[:, 1] % 4) * 4 + (coords[:, 2] % 4)
SUBLATTICES = [np.where(key == k)[0] for k in np.unique(key)]

rng = np.random.default_rng(160)
s = rng.normal(size=(N, 3))
s /= np.linalg.norm(s, axis=1)[:, None]
beta = 1.0 / T


def field(idx):
    return -(K1 * s[NN[idx]].sum(1) + K2 * s[NNN[idx]].sum(1))


def heat_bath(idx):
    h = field(idx)
    hm = np.linalg.norm(h, axis=1)
    hh = h / hm[:, None]
    x = beta * hm
    u = rng.random(len(idx))
    c = np.clip(1 + np.log(u + (1 - u) * np.exp(-2 * x)) / x, -1, 1)
    ref = np.where(np.abs(hh[:, :1]) < 0.9, np.array([[1.0, 0, 0]]), np.array([[0, 1.0, 0]]))
    e1 = np.cross(hh, ref)
    e1 /= np.linalg.norm(e1, axis=1)[:, None]
    e2 = np.cross(hh, e1)
    phi = 2 * np.pi * rng.random(len(idx))
    st = np.sqrt(1 - c * c)
    s[idx] = c[:, None] * hh + st[:, None] * (np.cos(phi)[:, None] * e1 + np.sin(phi)[:, None] * e2)


def over_relax(idx):
    h = field(idx)
    hh = h / np.linalg.norm(h, axis=1)[:, None]
    s[idx] = 2 * (s[idx] * hh).sum(1)[:, None] * hh - s[idx]


def sweep():
    for sub in SUBLATTICES:
        heat_bath(sub)
    for _ in range(4):
        for sub in SUBLATTICES:
            over_relax(sub)


def write(path, k):
    site = {(0, 0, 0): 1, (1, 1, 0): 2, (1, 0, 1): 3, (0, 1, 1): 4}
    lines = [f"TITLE MnO_160K_J1J2_MC_{k:02d}", f"CELL {A} {A} {A} 90 90 90", "SITE 0 0 0", "SITE 0.5 0.5 0", "SITE 0.5 0 0.5",
             "SITE 0 0.5 0.5", f"BOX {L} {L} {L}", f"FORM_FACTOR_J0 {J0}"]
    for (x, y, z), m in zip(coords, s * MU):
        lines.append(f"SPIN {site[(x % 2, y % 2, z % 2)]} {x // 2} {y // 2} {z // 2} {m[0]:.3f} {m[1]:.3f} {m[2]:.3f}")
    path.write_text("\n".join(lines) + "\n", encoding="ascii", newline="\n")


def main():
    for _ in range(1000):
        sweep()
    for k in range(1, SNAPSHOTS + 1):
        for _ in range(20):
            sweep()
        write(OUT / f"MnO_160K_spins_{k:02d}.txt", k)
    c1 = (s[:, None, :] * s[NN]).sum(-1).mean()
    c2 = (s[:, None, :] * s[NNN]).sum(-1).mean()
    print(f"wrote {SNAPSHOTS} snapshots of {N} spins; <s.s> nearest {c1:+.3f}, next-nearest {c2:+.3f}")


if __name__ == "__main__":
    main()
