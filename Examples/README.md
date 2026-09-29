# Diffuse Scattering Example Files

This directory contains example structures for the browser-based
diffuse-scattering calculator: `.rmc6f` profiles, magnetic configurations
and LAMMPS spin data.

Included files:

- `LiFeO2.rmc6f`
  - Chemical-order benchmark used to explore diffuse manifolds and `1/2(111)`-type condensation.
- `CaTiO3.rmc6f`
  - Displacement benchmark with rod-like and localized diffuse features.
- `PMN_300k.rmc6f`
  - PMN at 300 K benchmark for anisotropic diffuse scattering and advanced slice placement.

Magnetic examples of MnO (`magnetic/`). The first three are 6×6×6 cells with
864 Mn²⁺ spins of 5 μB:

- `MnO_order_spins.txt`
  - Scatty spin file of the type-II antiferromagnet; its magnetic Bragg
    peaks sit at (½ ½ ½)-type positions.
- `MnO_random_spins.txt`
  - Scatty spin file with random spin directions, an ideal paramagnet (as at
    infinite temperature).
- `MnO.rmc6f`, `MnO_spins.cfg`, `MnO.dat`
  - The antiferromagnet in RMCProfile's format: the structure, a legacy spin
    configuration, and the `ATOMS` line and `MAGNETISM` block (moments and
    ⟨j0⟩ form factor) of a `.dat` file. Load the three files together.
- `MnO_160K_spins_01.txt` to `MnO_160K_spins_08.txt`
  - Paramagnetic MnO at 160 K: Monte Carlo snapshots of 16×16×16 cells
    (16384 spins) of the J1–J2 Heisenberg model with the couplings Paddison
    refined to single-crystal data at 160 K (Spinteract, J1 = 3.26 K,
    J2 = 4.45 K). Selected together, they are averaged by a compute.
    `tools/make_mno_160k_spins.py` writes them.

The antiferromagnet and random Scatty files are copies of
`tests/fixtures/magnetic`, where the results are checked against `scatty.exe`;
`tools/make_magnetic_examples.py` writes them and the RMCProfile files.

LAMMPS spin examples (`lammps/`), data files with `atom_style spin` in the
current column order (`id type x y z spx spy spz sp`, units metal):

- `fe_bcc_spins.data`
  - bcc Fe, 10×10×10 cells (2000 atoms, a = 2.8665 Å), 2.2 μB spins tilted
    about +z with ⟨cos θ⟩ = 0.8, and Gaussian displacements (σ = 0.07 Å).
    Ferromagnetic Bragg peaks sit on the nuclear ones, on top of paramagnetic
    and thermal diffuse scattering.
- `skyrmion_film.data`
  - A triangular Fe monolayer of 60×60 sites (a = 2.5 Å) in a triclinic box,
    carrying a triple-q Bloch skyrmion lattice with a period of 10 sites
    (topological charge −1 per magnetic cell) and 2.5 μB moments: six
    magnetic satellites around every Bragg point (in the example's grid,
    around the origin) and rods along l, since the film is one layer thick.
    The example button computes magnetic scattering only, on a linear scale.

The calculator takes the whole LAMMPS box as the cell, so h, k, l count
reciprocal-box units: bcc (1 1 0) is (10 10 0) here, and the skyrmion
satellites lie 6 units from the origin, at (6 0 0), (0 6 0), (6 6 0) and
their opposites. A periodic box scatters only at integer box units, so the
examples, and LAMMPS data in general, use step 1. Between integer units, a
finite box shows only its size fringes. `tools/make_lammps_examples.py`
writes both files and checks the skyrmion texture's topological charge. They
are generated rather than copied from LAMMPS, whose own examples are GPL.

For workflow guidance and scientific context, see:

- [Diffuse-scattering examples guide](../docs/diffuse_scattering_examples.html)
