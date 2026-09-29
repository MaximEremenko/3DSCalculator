"""Writes the magnetic examples in Examples/magnetic from the MnO Scatty fixtures.

- MnO_order_spins.txt, MnO_random_spins.txt: the Scatty spin files of
  tests/fixtures/magnetic (made by make_mno_spins.py there).
- MnO.rmc6f, MnO_spins.cfg, MnO.dat: the ordered configuration in RMCProfile's
  magnetic format. MnO.rmc6f is rocksalt MnO (Mn first, then O), MnO_spins.cfg
  a legacy spin configuration with one unit vector per Mn atom in the same
  order, and MnO.dat holds the ATOMS line and the MAGNETISM block (5 muB and
  the Mn2+ <j0> form factor), the parts of a .dat that 3DSCalculator reads.

Usage (from the repository root):
    python tools/make_magnetic_examples.py
"""
import math
import shutil
from pathlib import Path

SOURCE = Path("tests/fixtures/magnetic")
TARGET = Path("Examples/magnetic")
MN_J0 = "0.4220 17.6840 0.5948 6.0050 0.0043 -0.6090 -0.0219"
O_SITES = [(0.5, 0, 0), (0, 0.5, 0), (0, 0, 0.5), (0.5, 0.5, 0.5)]


def read_scatty(path):
    cell, box, sites, spins = None, None, [], []
    for line in path.read_text(encoding="utf-8").splitlines():
        t = line.split()
        if not t:
            continue
        if t[0] == "CELL":
            cell = [float(v) for v in t[1:7]]
        elif t[0] == "BOX":
            box = [int(v) for v in t[1:4]]
        elif t[0] == "SITE":
            sites.append([float(v) for v in t[1:4]])
        elif t[0] == "SPIN":
            spins.append((int(t[1]), [int(v) for v in t[2:5]], [float(v) for v in t[5:8]]))
    return cell, box, sites, spins


def main():
    TARGET.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(SOURCE / "mno_order_spins_01.txt", TARGET / "MnO_order_spins.txt")
    shutil.copyfile(SOURCE / "mno_random_spins_01.txt", TARGET / "MnO_random_spins.txt")

    cell, box, sites, spins = read_scatty(SOURCE / "mno_order_spins_01.txt")
    assert cell[3:] == [90, 90, 90] and len(sites) == 4, "expects the cubic fcc MnO fixture"
    n = box[0] * box[1] * box[2]
    atoms = []
    for site, r, _ in spins:
        frac = [(r[d] + sites[site - 1][d]) / box[d] for d in range(3)]
        atoms.append(("Mn", frac, site, r))
    for cx in range(box[0]):
        for cy in range(box[1]):
            for cz in range(box[2]):
                for s, off in enumerate(O_SITES, start=5):
                    r = (cx, cy, cz)
                    atoms.append(("O", [(r[d] + off[d]) / box[d] for d in range(3)], s, r))
    n_mn = len(spins)
    lengths = [cell[d] * box[d] for d in range(3)]
    header = [
        "(Version 6f format configuration file)",
        "Metadata title:     MnO type-II antiferromagnet (from the Scatty example)",
        "Atom types present:          Mn O",
        f"Number of each atom type:    {n_mn} {len(atoms) - n_mn}",
        f"Number of atoms:             {len(atoms)}",
        f"Number density (Ang^-3):     {len(atoms) / (lengths[0] * lengths[1] * lengths[2]):.6f}",
        f"Supercell dimensions:        {box[0]} {box[1]} {box[2]}",
        f"Cell (Ang/deg):    {lengths[0]:.6f} {lengths[1]:.6f} {lengths[2]:.6f} 90.000000 90.000000 90.000000",
        "Atoms:",
    ]
    rows = [
        f"{i + 1}  {el}  {f[0]:.8f}  {f[1]:.8f}  {f[2]:.8f}  {site}  {r[0]}  {r[1]}  {r[2]}"
        for i, (el, f, site, r) in enumerate(atoms)
    ]
    (TARGET / "MnO.rmc6f").write_text("\n".join(header + rows) + "\n", encoding="utf-8", newline="\n")

    density = n_mn / (lengths[0] * lengths[1] * lengths[2])
    lines = [f"{density:.6f} {n_mn} 1", str(n_mn)]
    for _, _, m in spins:
        norm = math.sqrt(sum(v * v for v in m))
        lines.append(" ".join(f"{v / norm:.12f}" for v in m))
    (TARGET / "MnO_spins.cfg").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")

    dat = f"""TITLE :: MnO type-II antiferromagnet ({box[0]}x{box[1]}x{box[2]}, {n_mn} Mn spins)
ATOMS :: Mn O

MAGNETISM :: YES
  > MAGNETISM_FILE_STEM :: MnO_spins
  > MAGNETIC_ATOMS :: 1
  > FORM_FACTOR :: 1 {MN_J0}
  > MAGNETIC_MOMENTS :: 5.0

END ::
"""
    (TARGET / "MnO.dat").write_text(dat, encoding="utf-8", newline="\n")
    print(f"wrote {TARGET}: {len(atoms)} atoms, {n_mn} spins")


if __name__ == "__main__":
    main()
