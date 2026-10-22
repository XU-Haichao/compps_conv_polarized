# Source provenance

The initial upstream history is retained through `23bda36eb88ad68659d0dbb15ed423c1f611ed8d` from `https://github.com/mitsza/compps_conv`.

The four native build inputs (`compps.f`, `mscomppsq.inc`, `mscomppsj.inc`, `standalone_logging.f`) are copied unchanged from the research project's frozen positive-linear, eight-point incident-energy integration snapshot `9cc8f25add5c`. This was selected by the October 2026 preliminary parameter-scan runtime. Earlier nodal and rejected alternative interpolation variants are not substituted for it.

Compared with upstream, the research work adds coupled polarized source iteration and both escaping surfaces, exact short-characteristic formal transport, BLAS source products, configurable mesh sizes, source-iteration diagnostics, positive linear incident-energy product integration, and related array/safety fixes. The new `polarized.py` removes project-specific paths and supplies a standalone build/input/output interface.

The retained upstream `README`, `compile.xcm`, and scalar local-model definition are historical upstream interface files. The new README describes the polarized standalone workflow. No observational data, fitted spectra, MCRT results, or personal server settings are included.

See `source_manifest.json` for hashes of the exact imported native snapshot. A hash identifies provenance, not a claim of global numerical convergence or a license grant.
