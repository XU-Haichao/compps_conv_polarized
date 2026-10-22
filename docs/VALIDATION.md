# Initial packaging checks

The four native source files were checked byte-for-byte against the selected frozen linear-q8 research snapshot; hashes are in `source_manifest.json`. Upstream history is preserved.

Local macOS arm64 / GNU Fortran 15.1 tests on 2026-10-22:

- Standalone shared-library compilation succeeded.
- A deliberately coarse smoke calculation (21 depth, 40 energy, 3 angle nodes; kT_e = 20 keV, τ = 1, 5 eV blackbody, μ = 0.5) produced finite I/Q on both surfaces and reached its source-iteration criterion. This is an interface test, not a converged benchmark spectrum.
- The unit-normalized blackbody input integral and basic physical Stokes bounds passed the interface tests.
- `tests/check_transport.py` separately checks exact constant-source attenuation, unscattered Lambert injection, spatial passivity, and signed BLAS source products against independent analytic/matrix expressions.

Always perform mesh refinement for the physical case of interest. No new all-temperature/all-depth convergence certification is implied by this release.
