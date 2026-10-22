# compps_conv_polarized

Restricted research extension of the convolutional compPS model for self-consistent, angle-resolved Stokes I and Q Comptonization spectra escaping both surfaces of a slab.

**This is a private derivative, not an independently relicensed compPS distribution.** Upstream contains no explicit license grant. See [LICENSING.md](LICENSING.md) before copying, sharing, or redistributing any source. Keeping a repository private does not establish redistribution permission.

## Origin and scope

The repository preserves the Git history of [mitsza/compps_conv](https://github.com/mitsza/compps_conv), based on compPS by Juri Poutanen and Roland Svensson and the convolution version by Michal Szanecki. The original source banner and original `README` remain. The modified Fortran includes thermal polarized redistribution, coupled I/Q source iteration, upper/lower escape spectra, stable short-characteristic transport, and positive linear incident-energy product integration (eight Gauss nodes per interval). This snapshot corresponds to the later linear-q8 parameter-scan implementation, not the earlier nodal-only prototype.

The standalone Python interface restricts geometry to a full plane-parallel slab with Lambert-law lower-boundary injection. It does not claim validation of every legacy geometry or parameter option in compPS. The internal polarized solver is Fortran; Python builds it and supplies spectra and analysis-friendly arrays.

## Requirements and quick start

Python ≥3.10, NumPy, GNU Fortran, and BLAS. On Debian/Ubuntu install `gfortran libblas-dev`; macOS requires a working gfortran toolchain and Apple command-line developer tools (Accelerate is used for BLAS).

```sh
python -m pip install -r requirements.txt
# If gfortran is not on PATH:
# export COMPPSC_GFORTRAN=/absolute/path/to/gfortran
python polarized.py --kte 20 --tau 1 --seed-kt 0.005 --mu 0.5 --output output/example
python -m unittest discover -s tests -v
PYTHONPATH=. python tests/check_transport.py
```

The first run compiles a shared library. Build filenames identify source, compiler, and `(depth, energy, angle)` mesh. `--resolution 51 200 7` selects the default starting mesh; it is not a universal convergence guarantee. Compiled files and calculated spectra are ignored by Git.

`spectra.csv` contains E in keV and I/Q photon-number densities per keV per sr per injected photon (for a unit-normalized seed). `audit.json` records parameters, input seed mass, and the scattering-order iteration diagnostic. Pi(E) is Q(E)/I(E); use masks where I is negligible.

```python
import numpy as np
from polarized import build, blackbody_bins, PolarizedComppsc
edges = np.geomspace(1e-8, 1e3, 1761)
model = PolarizedComppsc(build(), edges, blackbody_bins(edges, .005), 20., 1.)
spectra, audits = model.window(30, 60, order=24)
I_up, Q_up, I_down, Q_down = spectra
```

For an arbitrary seed, provide **bin-integrated photon numbers**, not E² times a differential spectrum. No seed renormalization is silently applied: output scales linearly with supplied seed mass. Ensure the input energy grid contains the seed and Comptonized spectrum; widen/refine it to test interpolation error. Maximum input bin count is 5000 due to the native output buffer.

## Conventions

- τ is the full vertical slab Thomson depth; kT_e and seed kT are in keV.
- The upper surface is opposite the injection boundary; the lower surface is the injection boundary.
- `spectrum(mu)` returns bin-integrated values **per unit mu**; divide by `2π ΔE` for per-sr spectral densities. `window()` already makes that conversion and integrates the angular window.
- The Python interface flips native compPS Q once: returned Q is positive parallel to the slab surface, matching the recursion/MCRT project convention. U vanishes by axisymmetry in this reference basis; this interface does not fit a sky polarization angle.
- Fortran COMMON state is not thread safe. This interface serializes native calls; use separate processes for parallel models. Do not bypass the lock with direct calls from another wrapper.

## Convergence and intended use

The interface rejects an unconverged source-order iteration. Passing that check does **not** imply energy/depth/angular grid convergence. Independently refine all three dimensions, input energy sampling, and escape-angle integration. Large optical depth and low temperature need special care. Native normalization, interpolation, and zero high-energy tails can affect comparisons; this repository does not assert exact photon conservation or universal benchmark accuracy.

The original `compile.xcm` / `lmodel_comppsc.dat` remain for the upstream scalar XSPEC entry point. They do **not** automatically register an I/Q/U XSPEC fitting model. The polarized API documented here is the standalone `xscomppspol2f` interface.

See [physical/numerical notes](docs/PHYSICS.md) and [provenance](docs/PROVENANCE.md). Source organization and documentation preparation were assisted by GPT 6.1 sol. Requests for access to the modified compPSc source may be directed to the corresponding author; any sharing remains subject to upstream rights.
