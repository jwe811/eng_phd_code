#!/usr/bin/env python3
"""Small repeatable runtime benchmarks for the main engines."""

from __future__ import annotations

import argparse
import concurrent.futures
import itertools
import os
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


def _command_env(command_threads: int | None):
    if command_threads is None:
        return None
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(command_threads)
    env.setdefault("OMP_DYNAMIC", "FALSE")
    return env


def run_case(case: BenchmarkCase, repeats: int, command_threads: int | None = None) -> tuple[float, float]:
    times: list[float] = []
    env = _command_env(command_threads)
    for _ in range(repeats):
        start = time.perf_counter()
        proc = subprocess.run(
            case.args,
            cwd=ROOT,
            env=env,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT,
            check=False,
        )
        elapsed = time.perf_counter() - start
        if proc.returncode != 0:
            raise RuntimeError(f"{' '.join(case.args)} failed with exit code {proc.returncode}")
        times.append(elapsed)
    return min(times), sum(times) / len(times)


def default_jobs(limit: int) -> int:
    return max(1, min(limit, 4, os.cpu_count() or 1))


def main() -> int:
    parser = argparse.ArgumentParser(description="Run small engine runtime benchmarks.")
    parser.add_argument("-r", "--repeats", type=int, default=3, help="Number of repeats per case. Default: 3.")
    parser.add_argument("-j", "--jobs", type=int, default=default_jobs(len(CASES)), help="Benchmark cases to run concurrently. Default: up to 4.")
    parser.add_argument(
        "--command-threads",
        type=int,
        help="OMP_NUM_THREADS for each benchmark binary. Default: 1 when --jobs > 1, otherwise inherit the environment.",
    )
    args = parser.parse_args()

    if args.repeats < 1:
        parser.error("--repeats must be at least 1")
    if args.jobs < 1:
        parser.error("--jobs must be at least 1")
    if args.command_threads is not None and args.command_threads < 1:
        parser.error("--command-threads must be at least 1")

    command_threads = args.command_threads
    if command_threads is None and args.jobs > 1:
        command_threads = 1

    print("Benchmarks")
    print("case                    best_s    avg_s")
    print("----------------------  --------  --------")
    if args.jobs == 1:
        results = [run_case(case, args.repeats, command_threads) for case in CASES]
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) as executor:
            results = list(
                executor.map(
                    run_case,
                    CASES,
                    itertools.repeat(args.repeats),
                    itertools.repeat(command_threads),
                )
            )

    for case, (best, average) in zip(CASES, results):
        print(f"{case.name:<22}  {best:8.4f}  {average:8.4f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
