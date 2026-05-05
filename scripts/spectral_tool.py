#!/usr/bin/env python3
"""Post-process transfer-matrix spectral exports."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from postprocess.spectral import (
    beta_from_derivative,
    read_csr,
    read_vector,
    result_paths,
    spectral_residual,
    transition_probability_errors,
)


def load_problem(args):
    csr_path, left_path, right_path = result_paths(args.L, args.M, args.mode, args.root)
    return read_csr(csr_path), read_vector(left_path), read_vector(right_path)


def cmd_audit(args):
    csr, left, right = load_problem(args)
    residual = spectral_residual(csr, args.fugacity, right)
    beta = beta_from_derivative(csr, args.fugacity, left, right)
    max_error, bad = transition_probability_errors(csr, args.fugacity, right)
    print(f"states: {csr.n_states}")
    print(f"spectral_residual: {residual:.12e}")
    print(f"beta_derivative_ratio: {beta:.12e}")
    print(f"max_transition_probability_error: {max_error:.12e}")
    print(f"bad_transition_rows: {len(bad)}")
    if args.fail_on_bad and (residual > args.tolerance or max_error > args.tolerance):
        raise SystemExit(1)


def cmd_transitions(args):
    csr, _left, right = load_problem(args)
    max_error, bad = transition_probability_errors(csr, args.fugacity, right)
    print(f"max_transition_probability_error: {max_error:.12e}")
    for state, total in bad[: args.limit]:
        print(f"{state} {total:.17g}")
    if len(bad) > args.limit:
        print(f"... {len(bad) - args.limit} more rows")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    for name, help_text, func in [
        ("audit", "Compute residual, beta derivative ratio, and transition normalization", cmd_audit),
        ("transition-check", "Print transition-probability normalization failures", cmd_transitions),
    ]:
        p = sub.add_parser(name, help=help_text)
        p.add_argument("-L", type=int, required=True)
        p.add_argument("-M", type=int, required=True)
        p.add_argument("-m", "--mode", type=int, choices=range(4), required=True)
        p.add_argument("-x", "--fugacity", type=float, required=True)
        p.add_argument("--root", default="data/TransferMatrix/TMresults")
        p.add_argument("--tolerance", type=float, default=1e-8)
        p.add_argument("--fail-on-bad", action="store_true")
        p.add_argument("--limit", type=int, default=20)
        p.set_defaults(func=func)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
