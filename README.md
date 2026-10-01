# 3DSCalculator

**Live tool:** https://maximeremenko.github.io/3DSCalculator/

A browser-based forward diffuse-scattering calculator. It computes the
three-dimensional diffuse intensity `I(h,k,l)` in reciprocal space directly
from atomistic structure files — RMCProfile `.rmc6f` configurations, LAMMPS
data files, or unified structure HDF5 — and visualizes the result as
interactive 3D isosurfaces, 3D slice planes, and 2D slice maps. Everything
runs client-side in the browser: structure files never leave your machine.

## Documentation

The `docs/` directory contains the full documentation set
(live entry point: [https://maximeremenko.github.io/3DSCalculator/docs/](https://maximeremenko.github.io/3DSCalculator/docs/)):

- [User Guide](https://maximeremenko.github.io/3DSCalculator/docs/diffuse_scattering.html) — UI walkthrough and workflow.
- [Examples](https://maximeremenko.github.io/3DSCalculator/docs/diffuse_scattering_examples.html) — worked examples with the shipped benchmark files.
- [Theory](https://maximeremenko.github.io/3DSCalculator/docs/diffuse_scattering_theory.html) — the equations behind the calculator.
- [Supplementary](https://maximeremenko.github.io/3DSCalculator/docs/diffuse_scattering_supplementary.html) — backend provenance and scattering-table families.
- [Bibliography](https://maximeremenko.github.io/3DSCalculator/docs/diffuse_scattering_bibliography.html) — formal references.
- [Troubleshooting](https://maximeremenko.github.io/3DSCalculator/docs/diffuse_scattering_troubleshooting.html) — common issues and fixes.

## Features

- **Input formats**: RMCProfile `.rmc6f` (with the spin `.cfg` and `.dat`
  of magnetic RMC), LAMMPS data files (`.data` / `.lmp` / `.lammps`,
  including `atom_style spin`), unified structure HDF5 (`.h5` / `.hdf5`),
  and Scatty spin or atom files (`.txt`), loaded via drag-and-drop or the
  file picker (several files at once). Several configurations of one
  system selected together (RMC runs, Monte Carlo or MD snapshots) are
  averaged. Intensity data as `h k l I [σ]`
  lists (Scatty `_sc_list.txt`, Spinteract single-crystal data) load
  straight into the viewers.
- **Comparison with data**: a calculation is fitted to loaded intensity
  data (scale, flat or linear-in-|Q| background, weights 1/σ², as in
  Spinteract), with R_wp, R and χ² per degree of freedom. The views switch
  between the scaled calculation and the data, which share one color range,
  and a difference map in a diverging color map.
- **Compute backends**: the default, `wgpuNUFFT type-1 (WebGPU)`, evaluates
  the uniform `h,k,l` grid as a type-1 NUFFT on the GPU (bundled
  [wgpu-web](https://github.com/MaximEremenko/wgpuNUFFT)), with all phases
  prepared in fp64 and every run spot-checked against exact sums. The other
  engines are the type-3 WebGPU path (`Type-3 + webgpufft`), CPU FFT
  (`type3NufftCpu`), a CPU f64 NUFFT (`type3NufftCpuDirect`), and an exact
  CPU f64 direct sum for small validation grids. If WebGPU is unavailable
  the app falls back to the CPU FFT backend automatically; large grids on
  the WebGPU paths run in chunks.
- **Interface**: a Setup panel (structure, reciprocal grid, scattering,
  computation) with a pinned Compute button (`Ctrl+Enter`), three
  resizable and maximizable views, and a View panel for slice, color, 3D and export
  tools; light and dark themes; example structures one click away
  (including MnO spin configurations in Scatty and RMCProfile formats,
  paramagnetic MnO at 160 K as 8 Monte Carlo snapshots, and
  LAMMPS spin data for bcc Fe and a skyrmion film); a drawer layout for
  tablets and phones.
- **Magnetic scattering** (neutrons): magnetic-only or nuclear + magnetic
  intensity from Scatty spin files, RMCProfile magnetic configurations
  (moments and form factors from the `.dat` MAGNETISM block) and LAMMPS
  spins. It uses the perpendicular projection |M⊥|², Brown ⟨j0⟩/⟨j2⟩ form
  factors for 155 ions, and optional ideal-paramagnet subtraction. A species
  table sets the form factor, C2 and moment of each magnetic species. It is
  validated against J. A. M. Paddison's Scatty on MnO and spin ice.
- **Modulated and incommensurate structures** (optional tools, shown for
  the Cr spin-density-wave example or on request): a displacive wave,
  spin-density wave or helix with any wavevector q, applied to the loaded
  structure; satellites come out at G ± q with the box's finite-size width.
- **Radiation types**: neutron (fast and grouped-exact models), X-ray
  (Waasmaier table), and electron scattering (neutral-atom tables:
  Lobato, Peng, Doyle, Weickenmeier, Kirkland; ionic Peng model with
  user-specified valences).
- **Calculation options**: average-lattice subtraction, normalization
  controls, Laue symmetrization (Scatty's 12 Laue classes), optional 3D
  smoothing (Lanczos or Chebyshev filter), and an experimental per-element
  filter.
- **Visualization**: interactive Plotly 3D isosurface and 3D slice-plane
  views, plus a 2D slice heatmap with three slice modes — axis-aligned
  slices, arbitrary normal-plane slices, and volume-average slabs.
  Log/linear scaling, six colormaps, adjustable display levels, and a
  detachable slice window.
- **Exports**: `.dat` (`h k l intensity` columns), `.json`, unified
  diffuse-data HDF5 (`.h5`), ParaView `.vtk` (structured HKL grid),
  Gaussian `.cube` (ChimeraX/ParaView), 2D slice SVG and CSV, 3D plane/iso
  SVG snapshots, and standalone interactive Plotly HTML files of the
  3D views.

**Precision note**: the WebGPU paths and the CPU FFT path use `fp32`
arithmetic. On the shipped examples the type-1 backend stays within about
1e-6 of the peak intensity of exact sums. `type3NufftCpuDirect` runs the NUFFT in `fp64`, with an
accuracy set by the tolerance (1e-9 by default); it is not an exact sum.
For critical scientific checks, compare against the `Exact sum` backend on
a small grid.

## Getting Started

Open `index.html` in a modern browser. All local modules
load through plain script tags (no local `fetch` or workers), so opening the
file directly from disk (`file://`) works in practice. From disk, the browser
does not let the page read `Examples/` itself, so the example buttons load the
compressed copies in `Examples/embedded` (written by
`tools/make_embedded_examples.py`). The exception is PMN (30 MB), which opens
the file dialog. A local web server is the most reliable route:

```
python -m http.server
```

then browse to `http://localhost:8000/index.html`.

Notes:

- **Internet is required on first load** for the CDN scripts: the
  WebGPU-NUFFT type-3 kernels and WebGPU smoothing (jsDelivr) and Plotly
  (cdn.plot.ly, loaded on demand). HDF5 support (`h5wasm`) and the type-1
  engine (`wgpu-web`) are bundled locally in `js/`.
  The NUFFT and FFT libraries are pinned (WebGPU-NUFFT `v0.1.0`,
  WebGPU-FFT commit `fa45c93`), so results do not change when those
  repositories move on.
- The WebGPU backend requires a WebGPU-capable browser (e.g. current
  Chrome or Edge). Other browsers automatically use the CPU FFT backend.

Open **Examples** in the Structure card and pick one (or load your own structure), check the
`h,k,l` grid, and press **Compute diffuse** (`Ctrl+Enter`).

## Examples

Three benchmark `.rmc6f` configurations ship in [`Examples/`](Examples/):

| File | Description |
| --- | --- |
| `LiFeO2.rmc6f` | Chemical-order benchmark with diffuse manifold plus `1/2(111)` condensation. |
| `CaTiO3.rmc6f` | Displacement benchmark with overlapping rod-like and breathing-related diffuse features. |
| `PMN_300k.rmc6f` | Relaxor benchmark for anisotropic diffuse features and complex slice exploration. |

Magnetic examples sit in [`Examples/magnetic/`](Examples/magnetic/) (MnO
spin configurations as Scatty files and in RMCProfile's format) and
[`Examples/lammps/`](Examples/lammps/) (LAMMPS `atom_style spin` data for
bcc Fe and a skyrmion film, written by `tools/make_lammps_examples.py`).
[`Examples/README.md`](Examples/README.md) describes each file.

## Tests

The `tests/` directory holds Node tests for the shared JavaScript modules
(unified HDF5 I/O and the diffuse amplitude core). Run them from the
repository root with Node 22 or newer:

```
node --test
```

## Provenance

This repository was extracted (with full git history, via `git filter-repo`)
from the [MaximEremenko/Utilities](https://github.com/MaximEremenko/Utilities)
monorepo, where the tool lived under `RMCProfileUtilities/Diffuse_Scattering/`.
The shared modules `js/unified_hdf5.js` and `js/h5wasm.js` are vendored from
that monorepo's `Format_Converter` component; `h5wasm` is NIST-developed
software (see `js/h5wasm-LICENSE.txt`). `js/wgpu_web.js` is the standalone
build of wgpu-web 0.3.2 attached to the
[MaximEremenko/wgpuNUFFT](https://github.com/MaximEremenko/wgpuNUFFT)
[v0.3.2 release](https://github.com/MaximEremenko/wgpuNUFFT/releases/tag/v0.3.2)
(commit `69a3da9`, built by `wgpu-web/build_standalone.py` in its CI). Companion
tools from the Utilities collection remain in the monorepo and are linked
from the documentation.

## Scope

This tool is a forward diffuse-scattering calculator and visualization
surface. It does not perform the MOSAIC inverse reconstruction workflow;
the MOSAIC paper is referenced in the docs for scientific benchmark context
only.

## License

Apache License 2.0 — see [LICENSE](LICENSE). The vendored `h5wasm` bundle
carries its own NIST license notice in `js/h5wasm-LICENSE.txt`; the vendored
`wgpu-web` bundle is Apache-2.0 (`js/wgpu_web-LICENSE.txt`).
