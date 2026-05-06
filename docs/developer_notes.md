# Developer Notes

This repo still carries some archival structure because exact parity matters
more than cosmetic rewrites. Prefer small, parity-checked refactors.

## Current Architecture

- `src/MASTER_TMcalc.c` owns transfer-matrix graph construction.
- `src/MASTER_MCsample.c` owns the modern Monte Carlo entry point for modes 0
  and 1 and dispatches modes 2 and 3 to integrated source-level 2SAP engines.
- `src/MASTER_CreatorAll.c` owns exhaustive exact-span generation.
- `src/mc_legacy/` contains source files that are still active dependencies,
  but are isolated because they retain archival global-state assumptions.
- `src/mc_2sap_common.c` holds shared helpers for 2SAP/Ham2SAP mode support.

## Refactor Rules

- Keep physics-preserving changes small and run `make parity-audit` after each
  meaningful step.
- Avoid binary-level namespacing. All active logic should be readable C or
  Python source in the repo.
- Prefer explicit structs over new global variables.
- Prefer contiguous allocation for large tables, while preserving indexing
  conventions until call sites can be modernized safely.
- Add benchmarks before changing any transition, section, or rejection-sampling
  logic.

## Next Good Targets

- Move 2SAP and Ham2SAP engine state into a shared `Mc2SapEngine` struct.
- Extract CreatorAll output logic shared by standard and Hamiltonian 2SAP.
- Replace fixed 2SAP capacities with validated runtime capacities where the
  archival math allows it.
- Continue replacing repeated linear key lookups with exact hash tables.
