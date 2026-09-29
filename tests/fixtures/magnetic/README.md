# Magnetic test fixtures

Reference data for `tests/magnetic.test.js`, computed with J. A. M. Paddison's
programs rather than with 3DSCalculator itself.

## Spin configurations

`mno_random_spins_01.txt` and `mno_order_spins_01.txt` are
[Scatty](https://doi.org/10.1107/S2053273318015632) spin files written by
`make_mno_spins.py` (seed 20260929).

- **Model:** the fcc Mn sublattice of MnO, using the parameters of
  Spinteract's `examples/MnO`: a = 4.4344 Å and the Mn²⁺ ⟨j0⟩ coefficients
  `0.4220 17.6840 0.5948 6.0050 0.0043 -0.6090 -0.0219`.
- **Supercell:** 6×6×6, which gives 864 spins of 5 μB.
- **Random file:** random spin directions, an ideal paramagnet.
- **Ordered file:** the type-II antiferromagnet, with spins
  ±5 μB [1 1 −2]/√6 and sign (−1)^(x+y+z).

## Scatty reference intensities

`scatty_mno_reference.json` holds `scatty.exe` intensities (version
14 Oct 2022, from `spinteract_140923`'s `scatty_14Oct22exe.zip`) in
barn sr⁻¹ per spin. The grid is h, k, l = −1…1 in steps of 1/6, with l
fastest.

The values are a subset of four runs on the ±2 volume. Each run used this
`scatty_config.txt`:

```
NAME vol
CENTRE 0 0 0
X_AXIS 2 0 0 12
Y_AXIS 0 2 0 12
Z_AXIS 0 0 2 12
WINDOW 0
RADIATION N
SYMMETRY -1
```

The `*_temp_subtract` series add `TEMP_SUBTRACT`, Scatty's subtraction of
the ideal-paramagnet term (2/3) Σ|μ|² f².

- **Why this setup is exact:** `WINDOW 0` and a step equal to the
  supercell resolution (1/6) make Scatty return exact intensities at the
  grid points.
- **Conventions:**
  - At Q = 0 Scatty uses the orientational average (2/3)|M|².
  - Its prefactor is (γr₀/2)² = 0.07265289 b.
  - It normalises per spin.
