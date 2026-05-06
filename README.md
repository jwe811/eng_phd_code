# SAP and 2SAP Research Toolkit

This repository studies self-avoiding polygons (SAPs) and two-polygon systems
(2SAPs) in finite `L x M` lattice tubes. It can build transfer matrices,
calculate connective constants and eigenvectors, generate random Monte Carlo
samples, exhaustively enumerate all polygons of a fixed span, and post-process
the resulting polygon files.

The main tools are:

- `bin/tm_master`: transfer-matrix construction and spectral calculations.
- `bin/mc_master`: Monte Carlo sampling for SAP, Hamiltonian SAP, 2SAP, and
  Hamiltonian 2SAP modes.
- `bin/creator_all`: exhaustive exact-span generation for all four modes.
- `scripts/uofs_tool.py`: conversion, validation, counting, canonicalization,
  deduplication, and contact-map analysis for generated files.
- `scripts/spectral_tool.py`: audits over exported transfer-matrix CSR data and
  eigenvectors.
- `scripts/topology_tool.py`: BFACF-style shrink labels, label tallies,
  split-by-label output, and basic linking-number calculations.

The current engines are source-level C. They do not use on-the-fly compilation,
`-Dmain=...`, `objcopy`, or symbol prefixing for the 2SAP modes.

[Ph.D. dissertation](https://harvest.usask.ca/items/021d9d39-cc85-4584-a7ca-2d594f462496)

## First Five Minutes

Build everything:

```bash
make
```

Run the full parity check:

```bash
make parity-audit
```

Generate a tiny exhaustive SAP dataset:

```bash
bin/creator_all -L 1 -M 1 -m 0 -s 2
```

Summarize and validate that file:

```bash
python3 scripts/uofs_tool.py summary data/CreatorAll/All_SAPs/AllSAPsL1M1span2num1.txt
python3 scripts/uofs_tool.py validate data/CreatorAll/All_SAPs/AllSAPsL1M1span2num1.txt -L 1 -M 1 -s 2
```

Generate one random 2SAP sample:

```bash
bin/mc_master -L 2 -M 1 -m 2 -s 2 -n 1 -r 303 -S 822100
```

Compute a transfer-matrix result:

```bash
bin/tm_master -L 2 -M 1 -m 0
```

## Requirements

- `gcc`
- `make`
- OpenMP support for `gcc` via `-fopenmp`
- Python 3
- Bash, for `scripts/verify_all.sh`

No third-party Python packages are required.

## Concepts

### Lattice Tubes

The tools use an `L x M` tube cross-section and a span direction. In generated
UofS output from this repo, the span direction is `x`, width `L` is stored on
`y`, and height `M` is stored on `z`.

The post-processing tools default to this convention. Use `--span-axis z` only
for older files that use `z` as the tube direction.

### Modes

| Mode | Name | Meaning |
| :--- | :--- | :------ |
| `0` | Standard SAP | One self-avoiding polygon |
| `1` | Hamiltonian SAP | One polygon visiting every vertex in each lattice section |
| `2` | 2SAP | Two polygons sharing the same physical lattice |
| `3` | 2SAP-Hamiltonian | Two polygons whose union satisfies the Hamiltonian condition |

### UofS Files

Most generated polygon files use `UofS` direction format:

- The first line is `UofS`.
- Each polygon starts with `x y z`.
- Direction codes follow, one per line.
- `-111` ends one polygon.
- `-999` ends the file.
- In 2SAP files, two consecutive polygons form one 2SAP object.

Direction codes:

| Code | Move |
| :--- | :--- |
| `1` | `+x` |
| `2` | `-x` |
| `3` | `+y` |
| `4` | `-y` |
| `5` | `+z` |
| `6` | `-z` |

## What Can This Repo Do?

### 1. Transfer-Matrix Calculations

Use `tm_master` when you want connective constants, amplitudes, eigenvectors,
or exported sparse transfer matrices.

Examples:

```bash
# Standard SAP connective constant for a 2 x 1 tube
bin/tm_master -L 2 -M 1 -m 0

# Hamiltonian SAP case
bin/tm_master -L 2 -M 1 -m 1

# 2SAP transfer matrix
bin/tm_master -L 2 -M 1 -m 2

# 2SAP-Hamiltonian transfer matrix
bin/tm_master -L 2 -M 2 -m 3

# Fixed fugacity solve plus Monte Carlo eigenvector export
bin/tm_master -L 2 -M 1 -m 2 -x 1.0 -E

# Larger exploratory run with explicit capacity overrides
bin/tm_master -L 3 -M 1 -m 2 -S 200000 -K 800000
```

Important flags:

| Flag | Description |
| :--- | :---------- |
| `-L <int>` | Lattice width. Must be non-negative. |
| `-M <int>` | Lattice height. Must be non-negative. |
| `-m <0..3>` | Simulation mode. |
| `-S <int>` | Override maximum one-polygon section capacity. |
| `-K <int>` | Override maximum composite 2SAP state capacity. |
| `-c <float>` | Power-method convergence threshold. Default: `1e-8`. |
| `-d` | Enable damping. Hamiltonian modes enable it automatically. |
| `-x <float>` | Manually set fugacity and skip root finding. |
| `-E` | Export transition-indexed right eigenvectors for Monte Carlo workflows. |
| `-h` | Print usage. |

Transfer-matrix output goes to:

```text
data/TransferMatrix/TMresults/
```

Files include:

- `L_Evector_L<L>M<M>_<mode>.txt`
- `R_Evector_L<L>M<M>_<mode>.txt`
- `CSR_L<L>M<M>_<mode>.bin`

Mode suffixes are `std`, `ham`, `2sap`, and `2sap_ham`.

Every generated TM result also gets a small `.meta` sidecar with the command,
mode, lattice, git commit, timestamp, and OpenMP thread count. These metadata
files are for reproducibility only; they do not change the numeric result files.

### 2. Monte Carlo Sampling

Use `mc_master` when you want random samples at a target span.

Examples:

```bash
# One small standard SAP sample
bin/mc_master -L 1 -M 1 -m 0 -s 2 -n 1 -r 301 -S 227001

# One Hamiltonian SAP sample
bin/mc_master -L 1 -M 1 -m 1 -s 2 -n 1 -r 302 -S 227002

# One 2SAP sample
bin/mc_master -L 2 -M 1 -m 2 -s 2 -n 1 -r 303 -S 822100

# One 2SAP-Hamiltonian sample
bin/mc_master -L 2 -M 1 -m 3 -s 2 -n 1 -r 304 -S 931250

# Larger standard SAP batch
bin/mc_master -L 2 -M 1 -m 0 -s 5 -n 50 -r 1 -S 12345
```

Important flags:

| Flag | Description |
| :--- | :---------- |
| `-L <int>` | Lattice width. |
| `-M <int>` | Lattice height. |
| `-m <0..3>` | Simulation mode. |
| `-s <int>` | Target span. Must be at least `2`. |
| `-n <int>` | Number of samples. |
| `-r <int>` | Run number used in output filenames. |
| `-S <uint>` | RNG seed. |

Monte Carlo samples go to:

```text
data/MonteCarlo/SAPs/
data/MonteCarlo/HamSAPs/
data/MonteCarlo/2SAPs/
data/MonteCarlo/Ham2SAPs/
```

Calculated transition-indexed eigenvectors go to:

```text
data/MonteCarlo/MC_Evectors/
```

Monte Carlo sample and eigenvector files also get `.meta` sidecars with run
parameters such as command, seed, run number, mode, span, and dominant
eigenvalue. The UofS coordinate files themselves are unchanged.

### 3. Exhaustive Exact-Span Generation

Use `creator_all` when you want every SAP or unordered 2SAP object of exactly
span `s`.

Examples:

```bash
# All standard SAPs of exact span 2 in a 1 x 1 tube
bin/creator_all -L 1 -M 1 -m 0 -s 2

# All Hamiltonian SAPs
bin/creator_all -L 1 -M 1 -m 1 -s 2

# All unordered 2SAP systems
bin/creator_all -L 2 -M 1 -m 2 -s 2

# All unordered 2SAP-Hamiltonian systems with a larger count cap
bin/creator_all -L 2 -M 1 -m 3 -s 2 -C 1000000

# Refuse if output would exceed a small cap
bin/creator_all -L 1 -M 1 -m 0 -s 2 -C 10
```

Important flags:

| Flag | Description |
| :--- | :---------- |
| `-L <int>` | Lattice width. |
| `-M <int>` | Lattice height. |
| `-m <0..3>` | Simulation mode. |
| `-s <int>` | Exact target span. Must be at least `2`. |
| `-C <uint>` | Maximum number of objects to write before refusing. Default: `100000`. |
| `-f` | Force output even when the count exceeds `-C`. |

CreatorAll output goes to:

```text
data/CreatorAll/All_SAPs/
data/CreatorAll/All_HamSAPs/
data/CreatorAll/All_2SAPs/
data/CreatorAll/All_Ham2SAPs/
```

Each file contains at most 10,000 objects and uses a `num<N>` suffix.

### 4. UofS Post-Processing

Use `scripts/uofs_tool.py` for generated polygon files.

```bash
# Summarize a file
python3 scripts/uofs_tool.py summary data/CreatorAll/All_SAPs/AllSAPsL1M1span2num1.txt

# Validate SAP output
python3 scripts/uofs_tool.py validate data/CreatorAll/All_SAPs/AllSAPsL1M1span2num1.txt -L 1 -M 1 -s 2

# Validate 2SAP output. The -p 2 means two polygons per object.
python3 scripts/uofs_tool.py validate data/CreatorAll/All_2SAPs/All2SAPsL2M1span2num1.txt -p 2 -L 2 -M 1 -s 2

# Convert directions to point coordinates
python3 scripts/uofs_tool.py to-points data/CreatorAll/All_SAPs/AllSAPsL1M1span2num1.txt /tmp/saps_points.txt

# Convert 2SAP output to LP format for topology tooling
python3 scripts/uofs_tool.py to-lp data/MonteCarlo/2SAPs/MC2SAPsL2M1span2run303num1.txt /tmp/links_lp.txt -p 2

# Count objects by total edge length
python3 scripts/uofs_tool.py count-edges data/CreatorAll/All_SAPs/AllSAPsL1M1span2num1.txt

# Count objects by span
python3 scripts/uofs_tool.py count-spans data/CreatorAll/All_SAPs/AllSAPsL1M1span2num1.txt

# Count unique 2SAP objects, treating swapped polygon order as equivalent
python3 scripts/uofs_tool.py unique data/CreatorAll/All_2SAPs/All2SAPsL2M1span2num1.txt -p 2 --unordered

# Contact counts by cyclic distance
python3 scripts/uofs_tool.py contacts data/CreatorAll/All_SAPs/AllSAPsL1M1span2num1.txt -r 1.0
```

### 5. Spectral Post-Processing

Use `scripts/spectral_tool.py` for CSR/eigenvector consistency checks.

```bash
bin/tm_master -L 2 -M 1 -m 2
python3 scripts/spectral_tool.py audit -L 2 -M 1 -m 2 -x 0.653914
python3 scripts/spectral_tool.py transition-check -L 2 -M 1 -m 2 -x 0.653914 --limit 10
```

### 6. Topology Post-Processing

Use `scripts/topology_tool.py` for BFACF-style shrink labels, label-file
workflows, and simple linking-number checks.

```bash
# Label each SAP with the conservative shrink heuristic
python3 scripts/topology_tool.py shrink-id data/CreatorAll/All_SAPs/AllSAPsL1M1span2num1.txt /tmp/shrink_ids.txt

# Print shrunk lengths
python3 scripts/topology_tool.py shrink-lengths data/CreatorAll/All_SAPs/AllSAPsL1M1span2num1.txt

# Tally any label file
python3 scripts/topology_tool.py tally-labels /tmp/shrink_ids.txt

# Split a UofS file into one output file per label
python3 scripts/topology_tool.py split-by-label data/CreatorAll/All_SAPs/AllSAPsL1M1span2num1.txt /tmp/shrink_ids.txt /tmp/split_by_label

# Compute projection-based linking numbers for 2SAP objects
python3 scripts/topology_tool.py linking-number data/CreatorAll/All_2SAPs/All2SAPsL2M1span2num1.txt --keep-going
```

The shrink classifier is deterministic by default and intentionally
conservative. It writes one single-token label per polygon:

- `unknot`: the polygon shrank below the chosen threshold, so the heuristic
  treats it as unknotted.
- `knot_or_unresolved`: the polygon did not shrink below the threshold. This
  may be a true knot, or it may be an unknot that this simple shrink pass did
  not simplify far enough.

## Build And Test Commands

```bash
make                  # Build bin/tm_master, bin/mc_master, and bin/creator_all
make tm               # Build only bin/tm_master
make sampler          # Build only bin/mc_master
make creator          # Build only bin/creator_all
make test             # Run a small CSR/eigenvector audit
make verify           # Run transfer-matrix benchmark table
make parity-audit     # Run TM and MC parity checks
make postprocess-test # Run smoke tests for Python post-processing tools
make quick-check      # Build everything and run fast smoke tests
make check            # Build everything and run the full local audit suite
make clean            # Remove build artifacts only: build/ and bin/
make clean-data       # Remove generated data output trees
make distclean        # Remove both build artifacts and generated data
```

`make clean` intentionally leaves `data/` alone so research outputs are not
deleted by accident. Use `make clean-data` when you deliberately want to remove
`data/TransferMatrix`, `data/MonteCarlo`, `data/CreatorAll`, and loose
`data/*.txt` / `data/*.bin` files.

The most important regression check is:

```bash
make parity-audit
```

It runs transfer-matrix benchmarks across all four modes for `2x1`, `2x2`, and
`3x1`, plus Monte Carlo benchmarks across all four modes.

The critical 2SAP sampler benchmark is:

```text
mode 2, 2 x 1 grid dominant eigenvalue = 9.455960990693537
```

## Output Directory Map

```text
data/
  TransferMatrix/
    TMresults/       CSR matrices and left/right eigenvectors
  MonteCarlo/
    MC_Evectors/     Transition-indexed eigenvectors from sampling runs
    SAPs/            Mode 0 samples
    HamSAPs/         Mode 1 samples
    2SAPs/           Mode 2 samples
    Ham2SAPs/        Mode 3 samples
  CreatorAll/
    All_SAPs/        Mode 0 exhaustive output
    All_HamSAPs/     Mode 1 exhaustive output
    All_2SAPs/       Mode 2 exhaustive output
    All_Ham2SAPs/    Mode 3 exhaustive output
```

## Source Layout

```text
src/
  MASTER_TMcalc.c          Unified transfer-matrix driver
  MASTER_MCsample.c        Unified Monte Carlo entry point
  MASTER_CreatorAll.c      Exhaustive exact-span generator
  tm_spectral.c            Shared OpenMP CSR spectral solver for TM
  mc_spectral.c            Transition-indexed spectral solver for 2SAP MC
  mc_sampler_weights.c     Precomputed rejection-sampling weights
  mc_paths.c               Canonical mode labels, prefixes, and output paths
  run_metadata.c           Reproducibility sidecars for generated outputs
  mc_builder.c             SAP graph/endhinge discovery for MC
  mc_2sap_integrated.c     Integrated mode 2 sampler
  mc_2sap_ham_integrated.c Integrated mode 3 sampler
  mc_legacy/               Source-level legacy helpers still included by MC

include/
  tm_runtime.h             TM runtime helpers and exact section hash table
  tm_spectral.h            TM spectral problem interface
  mc_globals.h             MC globals and PolygonState definition
  mc_spectral.h            2SAP MC spectral interface
  mc_sampler_weights.h     Sampler weight interface
  mc_paths.h               Mode path/name interface
  mc_2sap_common.h         Shared helpers for 2SAP/Ham2SAP samplers
  run_metadata.h           Output metadata interface

deps/
  topology/, utils/        Topology and allocation helpers used by the engines
  numerics/                Small numerical routines shared by engine code

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
  historical vector conventions.
- The TM spectral solver prepares transpose CSR once so left-vector
  multiplication can run in parallel without atomic updates.
- Monte Carlo sampling uses transition-indexed right eigenvectors and
  precomputed left/right rejection envelopes.
- Modes 2 and 3 are source-level integrated C implementations linked into
  `bin/mc_master`.
