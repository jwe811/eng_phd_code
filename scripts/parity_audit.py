#!/usr/bin/env python3
"""Automated parity audit suite for the unified SAP/2SAP engines."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import itertools
import os
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


@dataclass(frozen=True)
class CreatorAllBenchmark:
    L: int
    M: int
    span: int
    mode: int
    name: str
    expected_count: int
    output_path: str
    output_sha256: str


@dataclass(frozen=True)
class TransferMatrixResult:
    bench: TransferMatrixBenchmark
    got: float
    status: str
    failed: bool


@dataclass(frozen=True)
class MonteCarloResult:
    bench: MonteCarloBenchmark
    eigen: float
    reject_first: int
    reject_last: int
    sha_status: str
    status: str
    failed: bool


@dataclass(frozen=True)
class CreatorAllResult:
    bench: CreatorAllBenchmark
    count: int
    sha_status: str
    validate_status: str
    status: str
    failed: bool


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
        output_path="data/MonteCarlo/SAPs/MCpolysL1M1span2run301num1.txt",
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
        output_path="data/MonteCarlo/HamSAPs/MCpolysHamL1M1span2run302num1.txt",
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
        output_path="data/MonteCarlo/2SAPs/MC2SAPsL2M1span2run303num1.txt",
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
        output_path="data/MonteCarlo/Ham2SAPs/MC2SAPsHamL2M1span2run304num1.txt",
        output_sha256="11726d6423c1213a180af142e0ef9ec466119279d75b02743a4d84deb53bc617",
    ),
)

SLOW_MC_BENCHMARKS = (
    MonteCarloBenchmark(
        L=1,
        M=1,
        span=2,
        samples=2,
        mode=0,
        name="standard-slow",
        run=401,
        seed=227101,
        expected_eigenvalue=7.0,
        eigen_tolerance=1e-14,
        reject_first=1,
        reject_last=5,
        output_path="data/MonteCarlo/SAPs/MCpolysL1M1span2run401num1.txt",
        output_sha256="a4addc34317a876704f178cdde56ecec4faa66c313614a06cf52866d49ac9332",
    ),
    MonteCarloBenchmark(
        L=1,
        M=1,
        span=2,
        samples=2,
        mode=1,
        name="hamiltonian-slow",
        run=402,
        seed=227102,
        expected_eigenvalue=3.732050852321649,
        eigen_tolerance=1e-6,
        reject_first=0,
        reject_last=0,
        output_path="data/MonteCarlo/HamSAPs/MCpolysHamL1M1span2run402num1.txt",
        output_sha256="e5e501396cb1d2c5482e033d175fbb504524afa112a59d8e5eccd664b90827a4",
    ),
    MonteCarloBenchmark(
        L=2,
        M=1,
        span=2,
        samples=2,
        mode=2,
        name="2sap-slow",
        run=403,
        seed=822101,
        expected_eigenvalue=9.455960990693537,
        eigen_tolerance=5e-14,
        reject_first=5,
        reject_last=9,
        output_path="data/MonteCarlo/2SAPs/MC2SAPsL2M1span2run403num1.txt",
        output_sha256="e1495ecffc13b2eb5f5a29980915870f0a5be2f670feb659cad7d9d77df753c6",
    ),
    MonteCarloBenchmark(
        L=2,
        M=1,
        span=2,
        samples=2,
        mode=3,
        name="2sap_ham-slow",
        run=404,
        seed=931251,
        expected_eigenvalue=5.534148126030995,
        eigen_tolerance=5e-14,
        reject_first=4,
        reject_last=2,
        output_path="data/MonteCarlo/Ham2SAPs/MC2SAPsHamL2M1span2run404num1.txt",
        output_sha256="174de6752c3ec70c3fcd4d294e1a9b6341fb5954fd5151b51209203864c07c96",
    ),
)


CREATOR_BENCHMARKS = (
    CreatorAllBenchmark(
        L=1,
        M=1,
        span=2,
        mode=0,
        name="standard",
        expected_count=170,
        output_path="data/CreatorAll/All_SAPs/AllSAPsL1M1span2num1.txt",
        output_sha256="b60e89bf08878691754fe6192deee13ed345beb63344631c30593bb64a098c13",
    ),
    CreatorAllBenchmark(
        L=1,
        M=1,
        span=2,
        mode=1,
        name="hamiltonian",
        expected_count=22,
        output_path="data/CreatorAll/All_HamSAPs/AllHamSAPsL1M1span2num1.txt",
        output_sha256="6cfd9cfbf5f6f9d7dc334fc4cf0be1e4b61e7a230862462d734f029f04f14480",
    ),
    CreatorAllBenchmark(
        L=2,
        M=1,
        span=2,
        mode=2,
        name="2sap",
        expected_count=624,
        output_path="data/CreatorAll/All_2SAPs/All2SAPsL2M1span2num1.txt",
        output_sha256="0fa7a4981d4ab4ceb13b500baa4172bf49adc2c18a0cc84aa8332a650e496732",
    ),
    CreatorAllBenchmark(
        L=2,
        M=1,
        span=2,
        mode=3,
        name="2sap_ham",
        expected_count=92,
        output_path="data/CreatorAll/All_Ham2SAPs/AllHam2SAPsL2M1span2num1.txt",
        output_sha256="760e5ae986da4c11fc4f150e5a940a7ffe83c159a4f91278106bda676d6c39ee",
    ),
)


def run_command(args: list[str], quiet: bool, command_threads: int | None = None) -> str:
    env = None
    if command_threads is not None:
        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(command_threads)
        env.setdefault("OMP_DYNAMIC", "FALSE")

    proc = subprocess.run(
        args,
        cwd=ROOT,
        env=env,
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


def default_jobs(limit: int | None = None) -> int:
    cpus = os.cpu_count() or 1
    jobs = min(4, cpus)
    if limit is not None:
        jobs = min(jobs, limit)
    return max(1, jobs)


def iter_ordered(items, jobs: int, func, *extra_args):
    if jobs <= 1:
        for item in items:
            yield func(item, *extra_args)
        return
    with concurrent.futures.ProcessPoolExecutor(max_workers=jobs) as executor:
        repeated_args = [itertools.repeat(arg) for arg in extra_args]
        yield from executor.map(func, items, *repeated_args)


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


def audit_transfer_matrix_case(bench: TransferMatrixBenchmark, command_threads: int | None) -> TransferMatrixResult:
    cmd = ["bin/tm_master", "-L", str(bench.L), "-M", str(bench.M), "-m", str(bench.mode)]
    try:
        output = clean_output(run_command(cmd, quiet=True, command_threads=command_threads))
        got = parse_float(r"Connective Constant(?:\s*\([^)]*\))?:\s+([0-9.eE+-]+)", output, "connective constant")
        assert_close(
            f"TM {bench.name} L={bench.L} M={bench.M}",
            got,
            bench.connective_constant,
            bench.tolerance,
        )
        return TransferMatrixResult(bench, got, "PASS", False)
    except Exception as exc:
        return TransferMatrixResult(bench, float("nan"), f"FAIL: {exc}", True)


def audit_transfer_matrix(benchmarks: tuple[TransferMatrixBenchmark, ...], jobs: int, command_threads: int | None) -> int:
    failures = 0
    print("\nTransfer-Matrix Parity")
    print("mode        lattice   got           expected      status")
    print("----------  --------  ------------  ------------  ------")

    for result in iter_ordered(benchmarks, jobs, audit_transfer_matrix_case, command_threads):
        bench = result.bench
        if result.failed:
            failures += 1

        print(
            f"{bench.name:<10}  {bench.L}x{bench.M:<5}  "
            f"{result.got:12.6f}  {bench.connective_constant:12.6f}  {result.status}"
        )

    return failures


def audit_monte_carlo_case(bench: MonteCarloBenchmark, command_threads: int | None) -> MonteCarloResult:
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
        output = clean_output(run_command(cmd, quiet=True, command_threads=command_threads))
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
        return MonteCarloResult(bench, eigen, reject_first, reject_last, "match", "PASS", False)
    except Exception as exc:
        return MonteCarloResult(bench, float("nan"), -1, -1, "mismatch", f"FAIL: {exc}", True)


def audit_monte_carlo(benchmarks: tuple[MonteCarloBenchmark, ...], jobs: int, command_threads: int | None) -> int:
    failures = 0
    print("\nMonte Carlo Sampler Parity")
    print("mode        lattice/span  eigenvalue         rejects   sha      status")
    print("----------  ------------  -----------------  --------  -------  ------")

    for result in iter_ordered(benchmarks, jobs, audit_monte_carlo_case, command_threads):
        bench = result.bench
        if result.failed:
            failures += 1

        print(
            f"{bench.name:<10}  {bench.L}x{bench.M}/s{bench.span:<4} "
            f"{result.eigen:17.15g}  {result.reject_first:2d}/{result.reject_last:<5d}  {result.sha_status:<7}  {result.status}"
        )

    return failures


def audit_creator_all_case(bench: CreatorAllBenchmark, command_threads: int | None) -> CreatorAllResult:
    cmd = [
        "bin/creator_all",
        "-L",
        str(bench.L),
        "-M",
        str(bench.M),
        "-s",
        str(bench.span),
        "-m",
        str(bench.mode),
    ]
    try:
        output = clean_output(run_command(cmd, quiet=True, command_threads=command_threads))
        count = parse_int(r"Generated\s+(\d+)\s+", output, "CreatorAll count")
        if count != bench.expected_count:
            raise AssertionError(f"count got {count}, expected {bench.expected_count}")

        output_path = ROOT / bench.output_path
        if not output_path.is_file():
            raise AssertionError(f"missing generated UofS file {bench.output_path}")

        digest = sha256_file(output_path)
        if digest != bench.output_sha256:
            raise AssertionError(f"sha256 got {digest}, expected {bench.output_sha256}")

        validation_output = run_command(
            [
                "python3",
                "scripts/uofs_tool.py",
                "validate",
                bench.output_path,
                "-L",
                str(bench.L),
                "-M",
                str(bench.M),
                "-s",
                str(bench.span),
            ],
            quiet=True,
        )
        failed = parse_int(r"failed_objects:\s+(\d+)", validation_output, "CreatorAll validation failures")
        if failed != 0:
            raise AssertionError(f"UofS validation reported {failed} failed objects")

        return CreatorAllResult(bench, count, "match", "ok", "PASS", False)
    except Exception as exc:
        return CreatorAllResult(bench, -1, "mismatch", "fail", f"FAIL: {exc}", True)


def audit_creator_all(benchmarks: tuple[CreatorAllBenchmark, ...], jobs: int, command_threads: int | None) -> int:
    failures = 0
    print("\nCreatorAll Parity")
    print("mode        lattice/span  count     sha      validate  status")
    print("----------  ------------  --------  -------  --------  ------")

    for result in iter_ordered(benchmarks, jobs, audit_creator_all_case, command_threads):
        bench = result.bench
        if result.failed:
            failures += 1

        print(
            f"{bench.name:<10}  {bench.L}x{bench.M}/s{bench.span:<4} "
            f"{result.count:8d}  {result.sha_status:<7}  {result.validate_status:<8}  {result.status}"
        )

    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description="Run archival parity audits for TM, MC, and CreatorAll engines.")
    parser.add_argument("--no-build", action="store_true", help="Skip the initial make all build step.")
    parser.add_argument("--tm-only", action="store_true", help="Run only transfer-matrix benchmarks.")
    parser.add_argument("--mc-only", action="store_true", help="Run only Monte Carlo sampler benchmarks.")
    parser.add_argument("--creator-only", action="store_true", help="Run only CreatorAll benchmarks.")
    parser.add_argument("--slow", action="store_true", help="Include slower deterministic multi-sample MC benchmarks.")
    parser.add_argument("--quiet", action="store_true", help="Suppress command output unless a test fails.")
    parser.add_argument("-j", "--jobs", type=int, default=default_jobs(), help="Benchmark subprocesses to run concurrently. Default: up to 4.")
    parser.add_argument(
        "--command-threads",
        type=int,
        help="OMP_NUM_THREADS for each audited binary. Default: 1 when --jobs > 1, otherwise inherit the environment.",
    )
    args = parser.parse_args()

    only_flags = [args.tm_only, args.mc_only, args.creator_only]
    if sum(1 for flag in only_flags if flag) > 1:
        parser.error("--tm-only, --mc-only, and --creator-only cannot be combined")

    if not args.no_build:
        run_command(["make", "all"], quiet=args.quiet)

    if args.jobs < 1:
        parser.error("--jobs must be at least 1")
    if args.command_threads is not None and args.command_threads < 1:
        parser.error("--command-threads must be at least 1")

    command_threads = args.command_threads
    if command_threads is None and args.jobs > 1:
        command_threads = 1

    failures = 0
    mc_benchmarks = MC_BENCHMARKS + (SLOW_MC_BENCHMARKS if args.slow else ())

    if args.creator_only:
        failures += audit_creator_all(CREATOR_BENCHMARKS, jobs=args.jobs, command_threads=command_threads)
    elif args.tm_only:
        failures += audit_transfer_matrix(TM_BENCHMARKS, jobs=args.jobs, command_threads=command_threads)
    elif args.mc_only:
        failures += audit_monte_carlo(mc_benchmarks, jobs=args.jobs, command_threads=command_threads)
    else:
        failures += audit_transfer_matrix(TM_BENCHMARKS, jobs=args.jobs, command_threads=command_threads)
        failures += audit_monte_carlo(mc_benchmarks, jobs=args.jobs, command_threads=command_threads)
        failures += audit_creator_all(CREATOR_BENCHMARKS, jobs=args.jobs, command_threads=command_threads)

    if failures:
        print(f"\nParity audit failed: {failures} benchmark(s) failed.", file=sys.stderr)
        return 1

    print("\nParity audit passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
