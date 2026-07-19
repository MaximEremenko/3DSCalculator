# 3DSCalc

A browser-based forward diffuse-scattering calculator. It computes the
three-dimensional diffuse intensity `I(h,k,l)` in reciprocal space directly
from atomistic structure files — RMCProfile `.rmc6f` configurations, LAMMPS
data files, or unified structure HDF5 — and visualizes the result as
interactive 3D isosurfaces, 3D slice planes, and 2D slice maps. Everything
runs client-side in the browser: structure files never leave your machine.

## Features

- **Input formats**: RMCProfile `.rmc6f`, LAMMPS data files
  (`.data` / `.lmp` / `.lammps`), and unified structure HDF5
  (`.h5` / `.hdf5`), loaded via drag-and-drop or file picker.
- **Compute backends**: type-3 NUFFT on a user-defined `h,k,l` grid with
  three selectable engines — WebGPU (`Type-3 + webgpufft`), CPU FFT
  (`type3NufftCpu`), and direct CPU summation (`type3NufftCpuDirect`).
  If WebGPU is unavailable the app falls back to the CPU FFT backend
  automatically; large grids on the WebGPU path run as chunked NUFFT passes.
- **Radiation types**: neutron (fast and grouped-exact models), X-ray
  (Waasmaier table), and electron scattering (neutral-atom tables:
  Lobato, Peng, Doyle, Weickenmeier, Kirkland; ionic Peng model with
  user-specified valences).
- **Calculation options**: average-lattice subtraction, normalization
  controls, optional 3D smoothing (Lanczos or Chebyshev filter), and an
  experimental per-element filter.
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

**Precision note**: the WebGPU FFT/NUFFT stack and the CPU FFT path use
`fp32` arithmetic. For critical scientific checks, compare against the
`type3NufftCpuDirect` backend.

## Getting Started

Open `index.html` in a modern browser. All local modules
load through plain script tags (no local `fetch` or workers), so opening the
file directly from disk (`file://`) works in practice; a local web server is
the most reliable route:

```
python -m http.server
```

then browse to `http://localhost:8000/index.html`.

Notes:

- **Internet is required on first load** for the CDN scripts: the
  WebGPU-NUFFT compute kernels (jsDelivr) and Plotly (cdn.plot.ly, loaded
  on demand). HDF5 support (`h5wasm`) is bundled locally in `js/`.
- The WebGPU backend requires a WebGPU-capable browser (e.g. current
  Chrome or Edge). Other browsers automatically use the CPU FFT backend.

Load one of the files from `Examples/` (or your own structure), set the
`h,k,l` grid range, and press **Compute Diffuse**.

## Documentation

The `docs/` directory contains the full documentation set
(entry point: [`docs/index.html`](docs/index.html)):

- [User Guide](docs/diffuse_scattering.html) — UI walkthrough and workflow.
- [Examples](docs/diffuse_scattering_examples.html) — worked examples with the shipped benchmark files.
- [Theory](docs/diffuse_scattering_theory.html) — the equations behind the calculator.
- [Supplementary](docs/diffuse_scattering_supplementary.html) — backend provenance and scattering-table families.
- [Bibliography](docs/diffuse_scattering_bibliography.html) — formal references.
- [Troubleshooting](docs/diffuse_scattering_troubleshooting.html) — common issues and fixes.

## Examples

Three benchmark `.rmc6f` configurations ship in [`Examples/`](Examples/):

| File | Description |
| --- | --- |
| `LiFeO2.rmc6f` | Chemical-order benchmark with diffuse manifold plus `1/2(111)` condensation. |
| `CaTiO3.rmc6f` | Displacement benchmark with overlapping rod-like and breathing-related diffuse features. |
| `PMN_300k.rmc6f` | Relaxor benchmark for anisotropic diffuse features and complex slice exploration. |

## Provenance

This repository was extracted (with full git history, via `git filter-repo`)
from the [MaximEremenko/Utilities](https://github.com/MaximEremenko/Utilities)
monorepo, where the tool lived under `RMCProfileUtilities/Diffuse_Scattering/`.
The shared modules `js/unified_hdf5.js` and `js/h5wasm.js` are vendored from
that monorepo's `Format_Converter` component; `h5wasm` is NIST-developed
software (see `js/h5wasm-LICENSE.txt`). Companion tools from the Utilities
collection remain in the monorepo and are linked from the documentation.

## Scope

This tool is a forward diffuse-scattering calculator and visualization
surface. It does not perform the MOSAIC inverse reconstruction workflow;
the MOSAIC paper is referenced in the docs for scientific benchmark context
only.

## License

Apache License 2.0 — see [LICENSE](LICENSE). The vendored `h5wasm` bundle
carries its own NIST license notice in `js/h5wasm-LICENSE.txt`.
