#!/usr/bin/env python3
"""Automated parity audit suite for the unified SAP/2SAP engines."""

from __future__ import annotations

import argparse
import hashlib
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ANSI_RE = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")


@dataclass(frozen=True)
class TransferMatrixBenchmark:
    L: int
    M: int
    mode: int
    name: str
    connective_constant: float
    tolerance: float = 1e-5


@dataclass(frozen=True)
class MonteCarloBenchmark:
    L: int
    M: int
    span: int
    samples: int
    mode: int
    name: str
    run: int
    seed: int
    expected_eigenvalue: float
    eigen_tolerance: float
    reject_first: int
    reject_last: int
    output_path: str
    output_sha256: str


TM_BENCHMARKS = (
    TransferMatrixBenchmark(2, 1, 0, "standard", 2.286331),
    TransferMatrixBenchmark(2, 1, 1, "hamiltonian", 1.553873),
    TransferMatrixBenchmark(2, 1, 2, "2sap", 1.529273),
    TransferMatrixBenchmark(2, 1, 3, "2sap_ham", 1.329970),
    TransferMatrixBenchmark(2, 2, 0, "standard", 2.731298),
    TransferMatrixBenchmark(2, 2, 1, "hamiltonian", 1.676259),
    TransferMatrixBenchmark(2, 2, 2, "2sap", 2.077460),
    TransferMatrixBenchmark(2, 2, 3, "2sap_ham", 1.556708),
    TransferMatrixBenchmark(3, 1, 0, "standard", 2.572051),
    TransferMatrixBenchmark(3, 1, 1, "hamiltonian", 1.629231),
    TransferMatrixBenchmark(3, 1, 2, "2sap", 1.907603),
    TransferMatrixBenchmark(3, 1, 3, "2sap_ham", 1.490649),
)


MC_BENCHMARKS = (
    MonteCarloBenchmark(
        L=1,
        M=1,
        span=2,
        samples=1,
        mode=0,
        name="standard",
        run=301,
        seed=227001,
        expected_eigenvalue=7.0,
        eigen_tolerance=1e-14,
        reject_first=0,
        reject_last=1,
        output_path="data/SAPs/MCpolysL1M1span2run301num1.txt",
        output_sha256="610f5528c6f5f2dd520056cf62404a1be49a52de4d9867aba12cadb5b4fcb9e2",
    ),
    MonteCarloBenchmark(
        L=1,
        M=1,
        span=2,
        samples=1,
        mode=1,
        name="hamiltonian",
        run=302,
        seed=227002,
        expected_eigenvalue=3.732050810014727,
        eigen_tolerance=1e-6,
        reject_first=0,
        reject_last=0,
        output_path="data/HamSAPs/MCpolysHamL1M1span2run302num1.txt",
        output_sha256="9c856b6f7816fb88a0e4658d8e2d6bfaf378ebdd580dc962bdec0ada6f58fecc",
    ),
    MonteCarloBenchmark(
        L=2,
        M=1,
        span=2,
        samples=1,
        mode=2,
        name="2sap",
        run=303,
        seed=822100,
        expected_eigenvalue=9.455960990693537,
        eigen_tolerance=5e-14,
        reject_first=2,
        reject_last=3,
        output_path="data/2SAPs/MC2SAPsL2M1span2run303num1.txt",
        output_sha256="06b3ccf6b7b4d8b9a58fa1fd51a292bfdde170d891063d359a0169189ad5370f",
    ),
    MonteCarloBenchmark(
        L=2,
        M=1,
        span=2,
        samples=1,
        mode=3,
        name="2sap_ham",
        run=304,
        seed=931250,
        expected_eigenvalue=5.534148126030995,
        eigen_tolerance=5e-14,
        reject_first=3,
        reject_last=3,
        output_path="data/Ham2SAPs/MC2SAPsHamL2M1span2run304num1.txt",
        output_sha256="11726d6423c1213a180af142e0ef9ec466119279d75b02743a4d84deb53bc617",
    ),
)


def run_command(args: list[str], quiet: bool) -> str:
    proc = subprocess.run(
        args,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    output = proc.stdout
    if proc.returncode != 0:
        if quiet:
            print(output, file=sys.stderr)
        raise RuntimeError(f"{' '.join(args)} failed with exit code {proc.returncode}")
    if not quiet:
        print(output, end="")
    return output


def clean_output(output: str) -> str:
    return ANSI_RE.sub("", output)


def parse_float(pattern: str, output: str, label: str) -> float:
    match = re.search(pattern, output)
    if not match:
        raise RuntimeError(f"could not parse {label}")
    return float(match.group(1))


def parse_int(pattern: str, output: str, label: str) -> int:
    match = re.search(pattern, output)
    if not match:
        raise RuntimeError(f"could not parse {label}")
    return int(match.group(1))


def assert_close(label: str, got: float, expected: float, tolerance: float) -> None:
    diff = abs(got - expected)
    if diff > tolerance:
        raise AssertionError(
            f"{label}: got {got:.15g}, expected {expected:.15g}, "
            f"diff {diff:.3e} exceeds tolerance {tolerance:.3e}"
        )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def audit_transfer_matrix(benchmarks: tuple[TransferMatrixBenchmark, ...], quiet: bool) -> int:
    failures = 0
    print("\nTransfer-Matrix Parity")
    print("mode        lattice   got           expected      status")
    print("----------  --------  ------------  ------------  ------")

    for bench in benchmarks:
        cmd = ["bin/tm_master", "-L", str(bench.L), "-M", str(bench.M), "-m", str(bench.mode)]
        try:
            output = clean_output(run_command(cmd, quiet=True))
            got = parse_float(r"Connective Constant(?:\s*\([^)]*\))?:\s+([0-9.eE+-]+)", output, "connective constant")
            assert_close(
                f"TM {bench.name} L={bench.L} M={bench.M}",
                got,
                bench.connective_constant,
                bench.tolerance,
            )
            status = "PASS"
        except Exception as exc:
            got = float("nan")
            status = f"FAIL: {exc}"
            failures += 1

        print(
            f"{bench.name:<10}  {bench.L}x{bench.M:<5}  "
            f"{got:12.6f}  {bench.connective_constant:12.6f}  {status}"
        )

    return failures


def audit_monte_carlo(benchmarks: tuple[MonteCarloBenchmark, ...], quiet: bool) -> int:
    failures = 0
    print("\nMonte Carlo Sampler Parity")
    print("mode        lattice/span  eigenvalue         rejects   sha      status")
    print("----------  ------------  -----------------  --------  -------  ------")

    for bench in benchmarks:
        cmd = [
            "bin/mc_master",
            "-L",
            str(bench.L),
            "-M",
            str(bench.M),
            "-s",
            str(bench.span),
            "-n",
            str(bench.samples),
            "-m",
            str(bench.mode),
            "-r",
            str(bench.run),
            "-S",
            str(bench.seed),
        ]
        try:
            output = clean_output(run_command(cmd, quiet=True))
            eigen = parse_float(
                r"Calculated (?:2SAP )?dominant eigenvalue:\s+([0-9.eE+-]+)",
                output,
                "dominant eigenvalue",
            )
            reject_first = parse_int(r"Rejected at first step\s+(\d+)\s+times", output, "first rejection count")
            reject_last = parse_int(r"Rejected at last step\s+(\d+)\s+times", output, "last rejection count")
            printed_path_match = re.search(r"printed UofS in file '([^']+)'", output)
            if not printed_path_match:
                raise RuntimeError("could not parse UofS output path")

            output_path = Path(printed_path_match.group(1))
            expected_path = Path(bench.output_path)
            if output_path != expected_path:
                raise AssertionError(f"output path {output_path} != expected {expected_path}")
            full_output_path = ROOT / output_path
            if not full_output_path.is_file():
                raise AssertionError(f"missing generated UofS file {output_path}")

            digest = sha256_file(full_output_path)
            assert_close(f"MC {bench.name} eigenvalue", eigen, bench.expected_eigenvalue, bench.eigen_tolerance)
            if reject_first != bench.reject_first or reject_last != bench.reject_last:
                raise AssertionError(
                    f"rejects got {reject_first}/{reject_last}, "
                    f"expected {bench.reject_first}/{bench.reject_last}"
                )
            if digest != bench.output_sha256:
                raise AssertionError(f"sha256 got {digest}, expected {bench.output_sha256}")
            status = "PASS"
            sha_status = "match"
        except Exception as exc:
            eigen = float("nan")
            reject_first = -1
            reject_last = -1
            sha_status = "mismatch"
            status = f"FAIL: {exc}"
            failures += 1

        print(
            f"{bench.name:<10}  {bench.L}x{bench.M}/s{bench.span:<4} "
            f"{eigen:17.15g}  {reject_first:2d}/{reject_last:<5d}  {sha_status:<7}  {status}"
        )

    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description="Run archival parity audits for TM and MC engines.")
    parser.add_argument("--no-build", action="store_true", help="Skip the initial make all build step.")
    parser.add_argument("--tm-only", action="store_true", help="Run only transfer-matrix benchmarks.")
    parser.add_argument("--mc-only", action="store_true", help="Run only Monte Carlo sampler benchmarks.")
    parser.add_argument("--quiet", action="store_true", help="Suppress command output unless a test fails.")
    args = parser.parse_args()

    if args.tm_only and args.mc_only:
        parser.error("--tm-only and --mc-only cannot be combined")

    if not args.no_build:
        run_command(["make", "all"], quiet=args.quiet)

    failures = 0
    if not args.mc_only:
        failures += audit_transfer_matrix(TM_BENCHMARKS, quiet=args.quiet)
    if not args.tm_only:
        failures += audit_monte_carlo(MC_BENCHMARKS, quiet=args.quiet)

    if failures:
        print(f"\nParity audit failed: {failures} benchmark(s) failed.", file=sys.stderr)
        return 1

    print("\nParity audit passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
