# Monte Carlo Compatibility Helpers

This directory holds the legacy helper routines that are still required by the
Monte Carlo sampler and CreatorAll integration. They are kept separate from the
shared `deps/` helpers because they depend on the Monte Carlo globals and, in
some cases, differ from similarly named transfer-matrix helpers.

The subdirectory names intentionally start with `legacy_` so these files do not
look like duplicate copies of the modern shared utilities:

- `legacy_sections/`: section-numbering routines used by the sampler.
- `legacy_topology/`: Monte Carlo versions of the left-flag tests.
- `legacy_transfer_matrix/`: legacy spectral helper still included by the
  sampler dependency wrapper.
- `legacy_utils/`: small helpers that are not yet shared elsewhere.

When a helper becomes independent of Monte Carlo globals, prefer moving it to a
shared `deps/` or `src/` module and deleting the compatibility copy.
