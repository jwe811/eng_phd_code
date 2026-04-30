# Master TMcalc Research Engine

A high-performance, modernized Transfer Matrix (TM) engine for calculating the Connective Constant and topological statistics of self-avoiding polygons (SAPs) and multi-polymer systems in arbitrary $L \times M$ lattice tubes.

## 🚀 Quick Start

1.  **Compile the engine**:
    ```bash
    make clean && make
    ```
2.  **Run a simulation** (e.g., 2x1 Hamiltonian mode):
    ```bash
    ./tm_master -L 2 -M 1 -m 1
    ```
3.  **Run the Audit Suite**:
    ```bash
    python3 audit_engine.py
    ```

## ⚙️ Configuration (CLI Flags)

The engine is fully dynamic and controlled via command-line arguments.

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

## 🔬 Simulation Modes

*   **Mode 0 (Standard)**: Standard Self-Avoiding Polygons (SAP) on an LxM lattice.
*   **Mode 1 (Hamiltonian)**: Hamiltonian walks that visit every vertex in the lattice section.
*   **Mode 2 (2SAP)**: Two-polymer system (e.g., dual interpenetrating walks).
*   **Mode 3 (2SAP-Hamiltonian)**: Two-polymer system where the union of both walks is Hamiltonian.

## 📊 Data Management & Export

The engine automatically exports spectral data to the `data/` directory using a mode-specific naming convention to prevent collisions:

- **CSR Matrices**: `data/CSR_L[L]M[M]_[mode].bin`
- **Eigenvectors**: `data/L_Evector_L[L]M[M]_[mode].txt`

*Suffixes: `std`, `ham`, `2sap`, `2sap_ham`.*

## 🛡️ Verification & Auditing

The `audit_engine.py` script provides a rigorous validation suite that performs:
1.  **Spectral Residue Check**: Validates the dominant eigenvector against the exported CSR matrix.
2.  **Analytical Alignment**: Compares code-calculated growth parameters (Alpha/Beta) against authoritative thesis specifications.
3.  **Pole Residue Convergence**: Verifies numerical stability and residue limits.

## 🛠 Project Structure

*   `src/MASTER_TMcalc.c`: Core engine logic (Hashing, Solver, CLI).
*   `audit_engine.py`: Python-based verification and audit suite.
*   `deps/topology/`: Connectivity validation routines.
*   `deps/utils/`: High-performance math and allocation utilities.
*   `data/`: Results directory (automatically cleaned via `make clean`).

## 🚀 Performance Features

*   **FNV-1a Hashing**: Dynamic state-space generation (replaces legacy hardcoded templates).
*   **OpenMP Acceleration**: Multi-threaded CSR matrix operations for rapid convergence.
*   **Real-time Telemetry**: Active monitoring of eigenvalue solver delta during runtime.
