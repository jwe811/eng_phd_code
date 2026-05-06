#!/usr/bin/env python3
"""Smoke tests for clear CLI failures on invalid inputs."""

from __future__ import annotations

import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class FailureCase:
    name: str
    args: tuple[str, ...]
    expected_text: str


FAILURE_CASES = (
    FailureCase(
        name="tm-negative-lattice",
        args=("bin/tm_master", "-L", "-1", "-M", "1", "-m", "0"),
        expected_text="must be >= 0",
    ),
    FailureCase(
        name="mc-invalid-mode",
        args=("bin/mc_master", "-L", "1", "-M", "1", "-m", "9", "-s", "2", "-n", "1"),
        expected_text="Invalid Simulation Mode",
    ),
    FailureCase(
        name="creator-span-too-small",
        args=("bin/creator_all", "-L", "1", "-M", "1", "-m", "0", "-s", "1"),
        expected_text="Span (-s) must be at least 2",
    ),
    FailureCase(
        name="mc-unsupported-2sap-lattice",
        args=("bin/mc_master", "-L", "3", "-M", "3", "-m", "2", "-s", "2", "-n", "1"),
        expected_text="Unsupported 2SAP L and M values",
    ),
)


def main() -> int:
    failures = 0

    print("CLI Smoke Tests")
    print("case                         status")
    print("---------------------------  ------")
    for case in FAILURE_CASES:
        proc = subprocess.run(
            case.args,
            cwd=ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        output = proc.stdout
        if proc.returncode == 0:
            print(f"{case.name:<27}  FAIL: exited 0")
            failures += 1
        elif case.expected_text not in output:
            print(f"{case.name:<27}  FAIL: missing '{case.expected_text}'")
            failures += 1
        else:
            print(f"{case.name:<27}  PASS")

    if failures:
        print(f"\nCLI smoke tests failed: {failures} case(s).", file=sys.stderr)
        return 1

    print("\nCLI smoke tests passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
