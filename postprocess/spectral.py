"""Spectral post-processing helpers for transfer-matrix exports."""

from __future__ import annotations

import math
import struct
from dataclasses import dataclass
from pathlib import Path


MODE_SUFFIXES = ["std", "ham", "2sap", "2sap_ham"]


@dataclass
class CsrMatrix:
    n_states: int
    row_ptr: list[int]
    out_states: list[int]
    edges: list[int]


def read_csr(path: str | Path) -> CsrMatrix:
    with Path(path).open("rb") as fp:
        header = fp.read(struct.calcsize("QQ"))
        if len(header) != struct.calcsize("QQ"):
            raise ValueError(f"{path} is too short to contain a CSR header")
        n_states, n_transitions = struct.unpack("QQ", header)
        row_ptr = list(struct.unpack("Q" * (n_states + 2), fp.read(8 * (n_states + 2))))
        out_states = list(struct.unpack("Q" * (n_transitions + 1), fp.read(8 * (n_transitions + 1))))
        edges = list(struct.unpack("Q" * (n_transitions + 1), fp.read(8 * (n_transitions + 1))))
    return CsrMatrix(n_states, row_ptr, out_states, edges)


def read_vector(path: str | Path) -> list[float]:
    return [float(line.strip()) for line in Path(path).read_text(encoding="utf-8").splitlines() if line.strip()]


def result_paths(lattice_l: int, lattice_m: int, mode: int, root: str | Path = "data/TransferMatrix/TMresults") -> tuple[Path, Path, Path]:
    suffix = MODE_SUFFIXES[mode]
    root = Path(root)
    return (
        root / f"CSR_L{lattice_l}M{lattice_m}_{suffix}.bin",
        root / f"L_Evector_L{lattice_l}M{lattice_m}_{suffix}.txt",
        root / f"R_Evector_L{lattice_l}M{lattice_m}_{suffix}.txt",
    )


def matvec(csr: CsrMatrix, fugacity: float, vec: list[float]) -> list[float]:
    out = [0.0] * csr.n_states
    for state in range(1, csr.n_states + 1):
        for idx in range(csr.row_ptr[state] + 1, csr.row_ptr[state + 1] + 1):
            out_state = csr.out_states[idx]
            out[state - 1] += (fugacity ** csr.edges[idx]) * vec[out_state - 1]
    return out


def spectral_residual(csr: CsrMatrix, fugacity: float, right: list[float]) -> float:
    image = matvec(csr, fugacity, right)
    return max(abs(a - b) for a, b in zip(image, right))


def beta_from_derivative(csr: CsrMatrix, fugacity: float, left: list[float], right: list[float]) -> float:
    inner = sum(l * r for l, r in zip(left, right))
    deriv = 0.0
    for state in range(1, csr.n_states + 1):
        for idx in range(csr.row_ptr[state] + 1, csr.row_ptr[state + 1] + 1):
            out_state = csr.out_states[idx]
            edge_count = csr.edges[idx]
            deriv += left[state - 1] * edge_count * (fugacity ** edge_count) * right[out_state - 1]
    return deriv / inner


def transition_probability_errors(csr: CsrMatrix, fugacity: float, right: list[float]) -> tuple[float, list[tuple[int, float]]]:
    bad: list[tuple[int, float]] = []
    max_error = 0.0
    for state in range(1, csr.n_states + 1):
        denom = right[state - 1]
        if denom == 0.0:
            continue
        total = 0.0
        for idx in range(csr.row_ptr[state] + 1, csr.row_ptr[state + 1] + 1):
            out_state = csr.out_states[idx]
            total += (fugacity ** csr.edges[idx]) * right[out_state - 1] / denom
        error = abs(total - 1.0)
        max_error = max(max_error, error)
        if error > 1e-8:
            bad.append((state, total))
    return max_error, bad

