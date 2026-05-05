# SAP and 2SAP Transfer Matrix / Monte Carlo Engines

This repository contains research code for self-avoiding polygons (SAPs) and
two self-avoiding polygon systems (2SAPs) in `L x M` lattice tubes. It provides
two unified executables:

- `bin/tm_master`: transfer-matrix state generation, CSR export, and spectral
  calculations for modes 0-3.
- `bin/mc_master`: Monte Carlo sampling for modes 0-3, including integrated
  2SAP and 2SAP-Hamiltonian logic.
- `bin/creator_all`: exhaustive exact-span generation for modes 0-3.

The current code is source-level C. It does not rely on on-the-fly compilation,
`-Dmain=...`, `objcopy`, or symbol prefixing to access 2SAP modes.

[Ph.D. dissertation](https://harvest.usask.ca/items/021d9d39-cc85-4584-a7ca-2d594f462496)

## Requirements

- `gcc`
- `make`
- OpenMP support (`-fopenmp`)
- Python 3 for audit scripts
- Bash for `scripts/verify_all.sh`

## Quick Start

Build everything:

```bash
make clean && make
```

Run a transfer-matrix calculation:

```bash
bin/tm_master -L 2 -M 1 -m 1
```

Generate Monte Carlo samples:

```bash
bin/mc_master -L 2 -M 1 -m 0 -s 5 -n 50 -r 1 -S 12345
```

Run the automated parity audit:

```bash
make parity-audit
```

## Modes

| Mode | Name | Meaning |
| :--- | :--- | :------ |
| `0` | Standard SAP | One self-avoiding polygon |
| `1` | Hamiltonian SAP | One polygon visiting every vertex in each lattice section |
| `2` | 2SAP | Two polygons sharing the same physical lattice |
| `3` | 2SAP-Hamiltonian | Two polygons whose union satisfies the Hamiltonian condition |

## Build Targets

```bash
make              # Build bin/tm_master, bin/mc_master, and bin/creator_all
make tm           # Build only bin/tm_master
make sampler      # Build only bin/mc_master
make creator      # Build only bin/creator_all
make test         # Run scripts/audit_engine.py for a small TM case
make verify       # Run transfer-matrix benchmark table
make parity-audit # Run TM and MC parity checks
make postprocess-test # Run smoke tests for revived post-processing tools
make clean        # Remove build artifacts and generated data subtrees
```

`make clean` removes `build/`, `bin/`, `data/TransferMatrix`,
`data/MonteCarlo`, `data/CreatorAll`, and loose `data/*.txt` / `data/*.bin`
files.

## Transfer Matrix CLI

Executable: `bin/tm_master`

| Flag | Description |
| :--- | :---------- |
| `-L <int>` | Lattice width. Must be non-negative; at least one of `L`, `M` must be positive. |
| `-M <int>` | Lattice height. Must be non-negative; at least one of `L`, `M` must be positive. |
| `-m <0..3>` | Simulation mode. |
| `-S <int>` | Override maximum number of one-polygon sections. |
| `-K <int>` | Override maximum number of composite 2SAP states. |
| `-c <float>` | Power-method convergence threshold. Default: `1e-8`. |
| `-d` | Enable damping in the spectral solver. Automatically enabled for Hamiltonian modes. |
| `-x <float>` | Manually set fugacity and skip root finding. Useful for sampling exports. |
| `-E` | Export transition-indexed right eigenvectors for Monte Carlo workflows. |
| `-h` | Print usage. |

Examples:

Calculate the standard SAP connective constant for a `2 x 1` tube:

```bash
bin/tm_master -L 2 -M 1 -m 0
```

Calculate the Hamiltonian SAP case for the same tube. Damping is enabled
automatically for Hamiltonian modes:

```bash
bin/tm_master -L 2 -M 1 -m 1
```

Calculate the 2SAP transfer matrix for a `2 x 1` tube:

```bash
bin/tm_master -L 2 -M 1 -m 2
```

Calculate the 2SAP-Hamiltonian case on a `2 x 2` tube:

```bash
bin/tm_master -L 2 -M 2 -m 3
```

Solve at a fixed fugacity instead of root-finding, then export the
transition-indexed right eigenvector used by Monte Carlo tooling:

```bash
bin/tm_master -L 2 -M 1 -m 2 -x 1.0 -E
```

Run a larger exploratory case with explicit capacity overrides:

```bash
bin/tm_master -L 3 -M 1 -m 2 -S 200000 -K 800000
```

## Monte Carlo CLI

Executable: `bin/mc_master`

| Flag | Description |
| :--- | :---------- |
| `-L <int>` | Lattice width. Must be non-negative; at least one of `L`, `M` must be positive. |
| `-M <int>` | Lattice height. Must be non-negative; at least one of `L`, `M` must be positive. |
| `-m <0..3>` | Simulation mode. Modes 2 and 3 dispatch to the integrated 2SAP samplers. |
| `-s <int>` | Target span. Must be at least `2`. |
| `-n <int>` | Number of samples. Must be positive. |
| `-r <int>` | Run number used in the output filename. |
| `-S <uint>` | RNG seed. |

Examples:

Generate one standard SAP sample on a `1 x 1` tube at span 2. This is the
small deterministic case used in the parity audit:

```bash
bin/mc_master -L 1 -M 1 -m 0 -s 2 -n 1 -r 301 -S 227001
```

Generate one Hamiltonian SAP sample for the same lattice and span:

```bash
bin/mc_master -L 1 -M 1 -m 1 -s 2 -n 1 -r 302 -S 227002
```

Generate one 2SAP sample on a `2 x 1` tube. This checks the integrated mode 2
path and the `9.455960990693537` dominant eigenvalue benchmark:

```bash
bin/mc_master -L 2 -M 1 -m 2 -s 2 -n 1 -r 303 -S 822100
```

Generate one 2SAP-Hamiltonian sample:

```bash
bin/mc_master -L 2 -M 1 -m 3 -s 2 -n 1 -r 304 -S 931250
```

Generate a larger batch of standard SAP samples at span 5. The run number and
seed make the output filename and RNG stream reproducible:

```bash
bin/mc_master -L 2 -M 1 -m 0 -s 5 -n 50 -r 1 -S 12345
```

## Exhaustive Creator CLI

Executable: `bin/creator_all`

`creator_all` generates every SAP or unordered 2SAP system of exactly span `s`
for the requested tube and mode. It uses the same graph/endhinge discovery and
deduplication tables as the Monte Carlo engine, but walks the transition graph
exhaustively instead of sampling.

| Flag | Description |
| :--- | :---------- |
| `-L <int>` | Lattice width. Must be non-negative; at least one of `L`, `M` must be positive. |
| `-M <int>` | Lattice height. Must be non-negative; at least one of `L`, `M` must be positive. |
| `-m <0..3>` | Simulation mode. |
| `-s <int>` | Exact target span. Must be at least `2`. |
| `-C <uint>` | Maximum number of objects to write before refusing. Default: `100000`. |
| `-f` | Force output even when the count exceeds `-C`. |

Generate all standard SAPs of exact span 2 in a `1 x 1` tube:

```bash
bin/creator_all -L 1 -M 1 -m 0 -s 2
```

Generate all Hamiltonian SAPs for the same tube and span:

```bash
bin/creator_all -L 1 -M 1 -m 1 -s 2
```

Generate all unordered 2SAP systems of exact span 2 in a `2 x 1` tube:

```bash
bin/creator_all -L 2 -M 1 -m 2 -s 2
```

Generate all unordered 2SAP-Hamiltonian systems, allowing a larger output:

```bash
bin/creator_all -L 2 -M 1 -m 3 -s 2 -C 1000000
```

Intentionally refuse output if the enumeration is larger than a small cap:

```bash
bin/creator_all -L 1 -M 1 -m 0 -s 2 -C 10
```

## Output Files

Generated files live under `data/`. Directories are created on demand.

### Transfer Matrix Output

`tm_master` writes transfer-matrix results to `data/TransferMatrix/TMresults/`:

- `L_Evector_L<L>M<M>_<mode>.txt`
- `R_Evector_L<L>M<M>_<mode>.txt`
- `CSR_L<L>M<M>_<mode>.bin`

Mode suffixes are:

- `std`
- `ham`
- `2sap`
- `2sap_ham`

With `-E`, `tm_master` instead writes the transition-indexed right eigenvector
format expected by Monte Carlo workflows:

- `R_Evector_TS_L<L>M<M>.txt`
- `R_EvectorHam_TS_L<L>M<M>.txt`
- `2SAP_R_Evector_TS_L<L>M<M>.txt`
- `2SAP_R_EvectorHam_TS_L<L>M<M>.txt`

### Monte Carlo Output

`mc_master` writes sampled polygons in `UofS` format:

- Mode 0: `data/MonteCarlo/SAPs/MCpolysL<L>M<M>span<S>run<R>num<N>.txt`
- Mode 1: `data/MonteCarlo/HamSAPs/MCpolysHamL<L>M<M>span<S>run<R>num<N>.txt`
- Mode 2: `data/MonteCarlo/2SAPs/MC2SAPsL<L>M<M>span<S>run<R>num<N>.txt`
- Mode 3: `data/MonteCarlo/Ham2SAPs/MC2SAPsHamL<L>M<M>span<S>run<R>num<N>.txt`

The sampler also writes calculated transition-indexed eigenvectors to
`data/MonteCarlo/MC_Evectors/` for auditing and reuse.

### Exhaustive Creator Output

`creator_all` writes `UofS` files under `data/CreatorAll/<MODE>/`:

- Mode 0: `data/CreatorAll/All_SAPs/AllSAPsL<L>M<M>span<S>num<N>.txt`
- Mode 1: `data/CreatorAll/All_HamSAPs/AllHamSAPsL<L>M<M>span<S>num<N>.txt`
- Mode 2: `data/CreatorAll/All_2SAPs/All2SAPsL<L>M<M>span<S>num<N>.txt`
- Mode 3: `data/CreatorAll/All_Ham2SAPs/AllHam2SAPsL<L>M<M>span<S>num<N>.txt`
- Each file contains at most 10,000 objects and uses a `num<N>` suffix:
  for example, `AllSAPsL1M1span2num1.txt`, `...num2.txt`, and so on.

The terminal output includes the total number of generated SAPs or unordered
2SAP systems.

## Verification

The most important regression check is:

```bash
make parity-audit
```

This runs:

- transfer-matrix benchmarks across all four modes for `2x1`, `2x2`, and `3x1`
  lattices;
- Monte Carlo benchmarks across all four modes;
- eigenvalue, rejection-count, output-path, and output-file SHA checks for the
  MC cases.

The critical 2SAP sampler benchmark is mode 2 on a `2 x 1` grid:

```text
dominant eigenvalue = 9.455960990693537
```

Additional checks:

```bash
make verify
make postprocess-test
python3 scripts/audit_engine.py -L 2 -M 1 -m 0
```

## Post-Processing Tools

Archive post-processing workflows have been revived as importable Python modules
under `postprocess/` and command-line tools under `scripts/`.

### UofS File Tools

`scripts/uofs_tool.py` handles file conversion, validation, counting,
canonicalization, deduplication, and contact maps.

Examples:

```bash
python3 scripts/uofs_tool.py summary data/CreatorAll/All_SAPs/AllSAPsL1M1span2num1.txt
python3 scripts/uofs_tool.py validate data/CreatorAll/All_SAPs/AllSAPsL1M1span2num1.txt -L 1 -M 1 -s 2
python3 scripts/uofs_tool.py to-points data/CreatorAll/All_SAPs/AllSAPsL1M1span2num1.txt /tmp/saps_points.txt
python3 scripts/uofs_tool.py to-lp data/MonteCarlo/2SAPs/MC2SAPsL2M1span2run303num1.txt /tmp/links_lp.txt -p 2
python3 scripts/uofs_tool.py unique data/CreatorAll/All_2SAPs/All2SAPsL2M1span2num1.txt -p 2 --unordered
python3 scripts/uofs_tool.py contacts data/CreatorAll/All_SAPs/AllSAPsL1M1span2num1.txt -r 1.0
```

For 2SAP files, pass `-p 2` so two consecutive UofS polygons are treated as one
object. The revived tools default to this repo's generated UofS convention:
span on `x`, lattice width `L` on `y`, and lattice height `M` on `z`. Use
`--span-axis z` for older files that use `z` as the tube direction.

### Spectral Tools

`scripts/spectral_tool.py` provides standalone audits over exported CSR matrices
and eigenvectors:

```bash
bin/tm_master -L 2 -M 1 -m 2
python3 scripts/spectral_tool.py audit -L 2 -M 1 -m 2 -x 0.653914
python3 scripts/spectral_tool.py transition-check -L 2 -M 1 -m 2 -x 0.653914 --limit 10
```

### Topology Tools

`scripts/topology_tool.py` revives the archive's BFACF-shrink-style
knot/unknot workflow and label-file utilities:

```bash
python3 scripts/topology_tool.py shrink-id data/CreatorAll/All_SAPs/AllSAPsL1M1span2num1.txt /tmp/shrink_ids.txt
python3 scripts/topology_tool.py shrink-lengths data/CreatorAll/All_SAPs/AllSAPsL1M1span2num1.txt
python3 scripts/topology_tool.py tally-labels /tmp/shrink_ids.txt
python3 scripts/topology_tool.py split-by-label data/CreatorAll/All_SAPs/AllSAPsL1M1span2num1.txt /tmp/shrink_ids.txt /tmp/split_by_label
python3 scripts/topology_tool.py linking-number data/CreatorAll/All_2SAPs/All2SAPsL2M1span2num1.txt --keep-going
```

The shrink classifier is deterministic by default and intentionally conservative:
objects that do not shrink below the archive threshold are labeled
`knot_or_unresolved`.

## Source Layout

```text
src/
  MASTER_TMcalc.c          Unified transfer-matrix driver
  MASTER_MCsample.c        Unified Monte Carlo entry point
  MASTER_CreatorAll.c      Exhaustive exact-span generator
  tm_spectral.c            Shared OpenMP CSR spectral solver for TM
  mc_spectral.c            Transition-indexed spectral solver for 2SAP MC
  mc_sampler_weights.c     Precomputed rejection-sampling weights
  mc_builder.c             SAP graph/endhinge discovery for MC
  mc_2sap_integrated.c     Integrated mode 2 sampler
  mc_2sap_ham_integrated.c Integrated mode 3 sampler

include/
  tm_runtime.h             TM runtime helpers and exact section hash table
  tm_spectral.h            TM spectral problem interface
  mc_globals.h             MC globals and PolygonState definition
  mc_spectral.h            2SAP MC spectral interface
  mc_sampler_weights.h     Sampler weight interface

deps/
  topology/, utils/        Topology and allocation helpers used by the engines
  mc_compat/               Legacy-compatible helper code included by MC

scripts/
  parity_audit.py          Automated TM and MC parity suite
  audit_engine.py          CSR/eigenvector consistency audit
  uofs_tool.py             UofS conversion, validation, counting, and contacts
  spectral_tool.py         CSR/eigenvector post-processing audits
  topology_tool.py         BFACF shrink labels and label-file utilities
  verify_all.sh            TM benchmark table

postprocess/
  uofs.py                  Shared UofS/LP parser, writer, validator
  analysis.py              Edge/span/contact analyses
  spectral.py              CSR/eigenvector audit helpers
  bfacf.py                 BFACF-style shrink heuristic
  topology.py              Linking-number helpers
```

## Architecture Notes

- Section ids are assigned from canonical boundary templates. Hashes accelerate
  lookup, but full templates are retained so hash collisions cannot merge
  states.
- 2SAP transfer-matrix states are composite pairs of one-polygon sections.
- Sparse matrices are stored in 1-based CSR form to remain compatible with the
  archival Numerical Recipes-style vector helpers.
- The TM spectral solver prepares transpose CSR once so left-vector
  multiplication can run in parallel without atomic updates.
- Monte Carlo sampling uses transition-indexed right eigenvectors and
  precomputed left/right rejection envelopes.
- Modes 2 and 3 are source-level integrated C implementations linked into
  `bin/mc_master`.
