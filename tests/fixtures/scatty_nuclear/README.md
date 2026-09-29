# Nuclear test fixture

Reference data for `tests/centring.test.js`, computed with J. A. M. Paddison's
[Scatty](https://doi.org/10.1107/S2053273318015632) rather than with
3DSCalculator itself.

## Configuration

`mnfeo_atoms_01.txt` is a Scatty atoms file written by `make_mnfeo.py`
(seed 20260929):

- **Structure:** rocksalt (Mn,Fe)O, a = 4.4344 Å, in a 5×5×5 box of 1000 atoms.
- **Disorder:** random Mn/Fe on the cation sites, with Mn and Fe shifted by
  ±0.012 lattice units along x, so occupancy and displacement are
  correlated. Every atom also has random displacements of 0.01 lattice
  units.
- **OCC lines:** the supercell's overall composition (Mn 0.502, Fe 0.498) on
  all four cation sites, so Scatty's average structure is F-centred.

## Scatty reference intensities

`scatty_mnfeo_reference.json` holds `scatty.exe` intensities (version
14 Oct 2022, from `spinteract_140923`'s `scatty_14Oct22exe.zip`) in barn sr⁻¹
per atom. The grid is h, k, l = −2…2 in steps of 0.2. The file stores all
125 lattice points and every 37th other point.

Both runs used this `scatty_config.txt`; the diffuse run added
`REMOVE_BRAGG F`:

```
CENTRE 0 0 0
X_AXIS 2 0 0 10
Y_AXIS 0 2 0 10
Z_AXIS 0 0 2 10
WINDOW 0
RADIATION N
SYMMETRY -1
EXPANSION_MAX_ERROR 1e-12
```

- **Why this setup is exact:** `WINDOW 0` and a step equal to the supercell
  resolution (1/5) make Scatty return exact intensities at the grid points.
- **Scattering lengths:** Scatty's own table has Fe 9.54 fm and O 5.805 fm,
  where Sears (1992), which 3DSCalculator uses, has 9.45 fm and 5.803 fm. The
  tests use Scatty's values.
- **Bragg removal:** with an F-centred average structure, Scatty removes the
  average only at F-allowed reflections. At F-forbidden lattice points, such
  as (0 0 1), the whole intensity is diffuse.
