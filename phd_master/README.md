# Master TMcalc Research Engine

A high-performance, modernized Transfer Matrix (TM) engine combined with Monte Carlo sampling tools for calculating the Connective Constant and topological statistics of self-avoiding polygons (SAPs) and multi-polymer systems in arbitrary $L \times M$ lattice tubes.

## 🚀 Quick Start

1.  **Compile all engines**:
    ```bash
    make clean && make
    ```
2.  **Run Transfer Matrix calculation** (e.g., 2x1 Hamiltonian mode):
    ```bash
    ./tm_master -L 2 -M 1 -m 1
    ```
3.  **Run Monte Carlo sampler**:
    ```bash
    ./mc_master -L 2 -M 1 -m 1
    ```
4.  **Run the Audit Suite**:
    ```bash
    python3 audit_engine.py
    ```

## ⚙️ Configuration (CLI Flags)

The engines are fully dynamic and controlled via command-line arguments.

### Transfer Matrix Engine (`tm_master`)

| Flag | Description | Values |
| :--- | :---------- | :----- |
| `-L` | **Lattice Width** | Positive integer (e.g., `-L 3`) |
| `-M` | **Lattice Height** | Positive integer (e.g., `-M 1`) |
| `-m` | **Simulation Mode** | `0`: Std, `1`: Ham, `2`: 2SAP, `3`: 2SAP-Ham |
| `-S` | **Max Sections** | Override max sections memory limit |
| `-K` | **Max Keynum** | Override composite key memory limit |
| `-c` | **Convergence** | Eigenvalue solver threshold (default `1e-8`) |
| `-d` | **Damping** | Enable damped power method (auto for Ham) |
| `-h` | **Help** | Displays usage summary |

### Monte Carlo Sampler (`mc_master`, `mc_2sap`, `mc_2sap_ham`)

| Flag | Description | Values |
| :--- | :---------- | :----- |
| `-L` | **Lattice Width** | Positive integer |
| `-M` | **Lattice Height** | Positive integer |
| `-m` | **Simulation Mode** | `0`: Std, `1`: Ham, `2`: 2SAP, `3`: 2SAP-Ham |
| `-s` | **Span Length** | Maximum walk extent |
| `-r` | **Run Number** | Identifier for this run |
| `-n` | **Num Samples** | Number of configurations to sample |
| `-h` | **Help** | Displays usage summary |

## 🏗️ Build Targets

```bash
make                   # Build all: tm_master, mc_master, mc_2sap, mc_2sap_ham
make tm                # Build only tm_master
make sampler           # Build all MC samplers: mc_master, mc_2sap, mc_2sap_ham
make test              # Run quick test suite (tm_master on 1x1, mode 0)
make verify            # Run full verification suite (bash verify_all.sh)
make clean             # Remove all executables and data/ contents
```

## 🔬 Simulation Modes

*   **Mode 0 (Standard)**: Standard Self-Avoiding Polygons (SAP) on an LxM lattice.
*   **Mode 1 (Hamiltonian)**: Hamiltonian walks that visit every vertex in the lattice section.
*   **Mode 2 (2SAP)**: Two-polymer system (e.g., dual interpenetrating walks).
*   **Mode 3 (2SAP-Hamiltonian)**: Two-polymer system where the union of both walks is Hamiltonian.

## 📊 Data Management & Export

### Transfer Matrix Engine
The `tm_master` engine exports spectral data to the `data/` directory using mode-specific naming conventions:
- **Eigenvectors**: `data/L_Evector_TS_L[L]M[M]_[mode].txt`, `data/R_Evector_TS_L[L]M[M]_[mode].txt`
- **CSR Matrices**: `data/CSR_L[L]M[M]_[mode].bin`

### Monte Carlo Samplers
The MC sampling tools (`mc_master`, `mc_2sap`, `mc_2sap_ham`) export sampled configurations to structured directories:
- **Standard SAPs**: `data/SAPs/MCpolysL[L]M[M]span[S]run[R]num[N].txt`
- **Hamiltonian SAPs**: `data/HamSAPs/MCpolysHamL[L]M[M]span[S]run[R]num[N].txt`
- **Two-polymer systems**: `data/2SAPs/MC2SAPsL[L]M[M]span[S]run[R]num[N].txt`
- **Two-polymer Hamiltonian**: `data/Ham2SAPs/MC2SAPsHamL[L]M[M]span[S]run[R]num[N].txt`

*Mode suffixes: `std`, `ham`, `2sap`, `2sap_ham`.*

## 🛡️ Verification & Auditing

The `audit_engine.py` script provides a rigorous validation suite that performs:
1.  **Spectral Residue Check**: Validates the dominant eigenvector against the exported CSR matrix.
2.  **Analytical Alignment**: Compares code-calculated growth parameters (Alpha/Beta) against authoritative thesis specifications.
3.  **Pole Residue Convergence**: Verifies numerical stability and residue limits.

## 🛠 Project Structure

### Executables Built
- **`tm_master`**: Transfer Matrix calculation engine (spectral solver for growth constants)
- **`mc_master`**: Monte Carlo sampler for standard SAP and Hamiltonian mode configurations
- **`mc_2sap`**: Specialized MC sampler for two-polymer systems (2SAP mode)
- **`mc_2sap_ham`**: Specialized MC sampler for two-polymer Hamiltonian systems (2SAP-Ham mode)

### Source Organization
- `src/MASTER_TMcalc.c`: Core TM engine (state generation via FNV-1a hashing, CSR matrix construction, eigenvalue solver)
- `src/mc_master.c`: Monte Carlo master entry point and runtime configuration
- `src/mc_builder.c`: Dynamic builder for MC sampler configurations
- `src/mc_globals.c`, `src/mc_globals.h`: Shared MC state and parameters
- `src/mc_validation.c`: Configuration validation and sanity checks
- `src/mc_utils.c`: Utility functions (RNG, I/O, memory management)
- `src/mc_memory.c`: Dynamic memory allocation for MC samplers
- `src/mc_deps.c`: Unified dependencies for MC tools (topology, utilities)
- `deps/topology/`: Connectivity validation routines (hinges, spanning walks)
- `deps/utils/`: High-performance utilities (vector allocation, sorting, matrix operations)
- `generated/`: Auto-generated MC sampler code (populated via `scripts/build_2sap_generic.py`)
- `audit_engine.py`: Python-based verification and audit suite
- `data/`: Results directory (automatically created and cleaned via `make clean`)

## 🚀 Performance Features

*   **FNV-1a Hashing**: Dynamic state-space generation (replaces legacy hardcoded templates).
*   **OpenMP Acceleration**: Multi-threaded CSR matrix operations for rapid convergence.
*   **Real-time Telemetry**: Active monitoring of eigenvalue solver delta during runtime.
