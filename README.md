# SAP and 2SAP Research Toolkit

This repository contains a high-performance C and Python toolkit for studying **Self-Avoiding Polygons (SAPs)** and **Two-Polygon Systems (2SAPs)** in finite $L \times M$ lattice tubes. 

The toolkit supports transfer-matrix construction, spectral calculations, exact exhaustive enumeration, and Monte Carlo sampling. It also includes an interactive web-based workbench for visualization and analysis.

---

## 🚀 Getting Started (From Scratch)

### 1. Prerequisites
Ensure you have the following installed on your system (Ubuntu/WSL2 recommended):
- **C Development**: `gcc`, `make`, and OpenMP support (typically included with `gcc`).
- **Python**: Python 3.10+ with `venv` support.
- **Web Development**: Node.js 18+ and `npm`.

### 2. Installation
Clone the repository and build the core C binaries:
```bash
git clone https://github.com/jwe811/eng_phd_code.git
cd eng_phd_code
make -j$(nproc)
```

### 3. Web Workbench Setup
The interactive workbench requires a Python virtual environment for the backend and `npm` packages for the frontend:

**Backend Setup:**
```bash
cd web/backend
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
cd ../..
```

**Frontend Setup:**
```bash
cd web/frontend
npm install
cd ../..
```

### 4. Launching the Workbench
Once setup is complete, you can launch both the backend and frontend with a single command:
```bash
make web-dev
```
- **Frontend**: [http://127.0.0.1:5173](http://127.0.0.1:5173)
- **Backend API**: [http://127.0.0.1:8000](http://127.0.0.1:8000)

---

## 🛠 Core Workflows (CLI)

### Exhaustive Generation
Generate every possible polygon of an exact span in a specific tube.
```bash
# All standard SAPs of span 2 in a 1x1 tube
bin/creator_all -L 1 -M 1 -m 0 -s 2
```

### Monte Carlo Sampling
Generate random samples using transition-indexed right eigenvectors.
```bash
# Generate 50 random samples at span 5 in a 2x1 tube
bin/mc_master -L 2 -M 1 -m 0 -s 5 -n 50 -r 1 -S 12345
```

### Transfer-Matrix Calculations
Calculate connective constants, amplitudes, and export eigenvectors.
```bash
# Standard SAP connective constant for a 2x1 tube
bin/tm_master -L 2 -M 1 -m 0
```

---

## 🌐 Interactive Workbench
The SAP Workbench provides a graphical interface for the entire research lifecycle:

- **Run Launcher**: Configure and start simulations with live command previews.
- **Job History**: Track background processes, view logs, and manage run results.
- **Data Browser**: Navigate generated files with built-in **Raw Text View**.
- **3D Visualize**: Interactive 3D rendering (Rotate, Pan, Zoom) of any SAP or 2SAP object.
- **Analysis**: Automatic generation of distribution charts for edges, spans, contacts, and topological linking numbers.

---

## 📖 Reference

### Simulation Modes
| Mode | Name | Meaning |
| :--- | :--- | :------ |
| `0` | Standard SAP | One self-avoiding polygon |
| `1` | Hamiltonian SAP | One polygon visiting every vertex in each lattice section |
| `2` | 2SAP | Two polygons sharing the same physical lattice |
| `3` | 2SAP-Hamiltonian | Two polygons whose union is Hamiltonian |

### UofS File Format
Most outputs use the `UofS` direction format. 
- **1-6**: Directions (+x, -x, +y, -y, +z, -z)
- **-111**: End of polygon
- **-999**: End of file

See [docs/uofs_format.md](docs/uofs_format.md) for full coordinate details.

### Output Directory Map
```text
data/
  TransferMatrix/
    TMresults/       CSR matrices and eigenvectors
  MonteCarlo/
    SAPs/            Mode 0 random samples
    2SAPs/           Mode 2 random samples
  CreatorAll/
    All_SAPs/        Mode 0 exhaustive generation
```

### Build Commands
- `make`: Build all core binaries.
- `make check`: Run the full local audit and parity suite.
- `make parity-audit`: Run regression checks against known benchmarks.
- `make clean`: Remove build artifacts.

---

## 🎓 Academic Context
This toolkit was developed as part of a Ph.D. research project on the topological and statistical properties of self-avoiding polygons.

- **Ph.D. Dissertation**: [Link to University of Saskatchewan Library](https://harvest.usask.ca/items/021d9d39-cc85-4584-a7ca-2d594f462496)
- **Coordinate Convention**: Tube direction is always `x`, width is `y`, and height is `z`.

---
*Developed by Dr. Jeremy Eng*