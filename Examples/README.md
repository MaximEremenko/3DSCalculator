# Diffuse Scattering Example Files

This directory contains example `.rmc6f` profiles for the browser-based
RMC6f diffuse-scattering calculator.

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

For workflow guidance and scientific context, see:

- [Diffuse-scattering examples guide](../docs/diffuse_scattering_examples.html)
