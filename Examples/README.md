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

Magnetic examples (`magnetic/`), 6×6×6 cells of MnO with 864 Mn²⁺ spins of 5 μB:

- `MnO_order_spins.txt`
  - Scatty spin file of the type-II antiferromagnet; its magnetic Bragg
    peaks sit at (½ ½ ½)-type positions.
- `MnO_random_spins.txt`
  - Scatty spin file with random spin directions, an ideal paramagnet.
- `MnO.rmc6f`, `MnO_spins.cfg`, `MnO.dat`
  - The antiferromagnet in RMCProfile's format: the structure, a legacy spin
    configuration, and the `ATOMS` line and `MAGNETISM` block (moments and
    ⟨j0⟩ form factor) of a `.dat` file. Load the three files together.

The Scatty files are copies of `tests/fixtures/magnetic`, where the results are
checked against `scatty.exe`. `tools/make_magnetic_examples.py` rewrites this
folder from those fixtures.

LAMMPS spin examples (`lammps/`), data files with `atom_style spin` in the
current column order (`id type x y z spx spy spz sp`, units metal):

- `fe_bcc_spins.data`
  - bcc Fe, 10×10×10 cells (2000 atoms, a = 2.8665 Å), 2.2 μB spins tilted
    about +z with ⟨cos θ⟩ = 0.8, and Gaussian displacements (σ = 0.07 Å).
    Ferromagnetic Bragg peaks sit on the nuclear ones, on top of paramagnetic
    and thermal diffuse scattering.
- `skyrmion_film.data`
  - A triangular Fe monolayer of 60×60 sites (a = 2.5 Å) in a triclinic box,
    carrying a triple-q Bloch skyrmion lattice with a period of 10 sites and
    2.5 μB moments: six magnetic satellites around every Bragg point (in the
    example's grid, around the origin) and rods along l, since the film is
    one layer thick.

The calculator takes the whole LAMMPS box as the cell, so h, k, l count
reciprocal-box units: bcc (1 1 0) is (10 10 0) here, and the skyrmion
satellites lie 6 units from the origin, at (6 0 0), (0 6 0), (6 6 0) and
their opposites. `tools/make_lammps_examples.py` writes both files. They are
generated rather than copied from LAMMPS, whose own examples are GPL.

For workflow guidance and scientific context, see:

- [Diffuse-scattering examples guide](../docs/diffuse_scattering_examples.html)
