#!/usr/bin/env python3
"""Small repeatable runtime benchmarks for the main engines."""

from __future__ import annotations

import argparse
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class BenchmarkCase:
    name: str
    args: tuple[str, ...]


CASES = (
    BenchmarkCase("tm-standard-2x1", ("bin/tm_master", "-L", "2", "-M", "1", "-m", "0")),
    BenchmarkCase("tm-2sap-2x1", ("bin/tm_master", "-L", "2", "-M", "1", "-m", "2")),
    BenchmarkCase("mc-2sap-2x1-s2", ("bin/mc_master", "-L", "2", "-M", "1", "-m", "2", "-s", "2", "-n", "1", "-r", "903", "-S", "822100")),
    BenchmarkCase("creator-2sap-2x1-s2", ("bin/creator_all", "-L", "2", "-M", "1", "-m", "2", "-s", "2")),
)


def run_case(case: BenchmarkCase, repeats: int) -> tuple[float, float]:
    times: list[float] = []
    for _ in range(repeats):
        start = time.perf_counter()
        proc = subprocess.run(
            case.args,
            cwd=ROOT,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT,
            check=False,
        )
        elapsed = time.perf_counter() - start
        if proc.returncode != 0:
            raise RuntimeError(f"{' '.join(case.args)} failed with exit code {proc.returncode}")
        times.append(elapsed)
    return min(times), sum(times) / len(times)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run small engine runtime benchmarks.")
    parser.add_argument("-r", "--repeats", type=int, default=3, help="Number of repeats per case. Default: 3.")
    args = parser.parse_args()

    if args.repeats < 1:
        parser.error("--repeats must be at least 1")

    print("Benchmarks")
    print("case                    best_s    avg_s")
    print("----------------------  --------  --------")
    for case in CASES:
        best, average = run_case(case, args.repeats)
        print(f"{case.name:<22}  {best:8.4f}  {average:8.4f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
