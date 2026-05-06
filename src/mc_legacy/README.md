# Monte Carlo Legacy Helpers

This directory holds helper routines that are still required by the Monte Carlo
sampler and CreatorAll integration. They live under `src/` because they are
active source-level dependencies, but they are kept separate from the modern MC
modules because they still depend on Monte Carlo globals and, in some cases,
differ from similarly named transfer-matrix helpers.

- `sections/`: section-numbering routines used by the sampler.
- `topology/`: Monte Carlo versions of the left-flag tests.
- `transfer_matrix/`: legacy spectral helper still included by the sampler
  dependency wrapper.
- `utils/`: small helpers that are not yet shared elsewhere.

When a helper becomes independent of Monte Carlo globals, prefer moving it to a
shared `deps/` or modern `src/` module and deleting the legacy copy.
