#!/usr/bin/env python3
"""Knot/link oriented post-processing tools."""

from __future__ import annotations

import argparse
import concurrent.futures
import itertools
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from postprocess.bfacf import shrink_label, shrink_length
from postprocess.topology import linking_number
from postprocess.uofs import read_objects


def _linking_number_line(args):
    idx, obj, keep_going = args
    try:
        value = linking_number(obj.polygons[0], obj.polygons[1])
        return f"{idx} {value:g}"
    except ValueError as exc:
        if keep_going:
            return f"{idx} unresolved {exc}"
        raise


def cmd_shrink_id(args):
    objects = read_objects(args.input, args.polys_per_object)
    with Path(args.output).open("w", encoding="utf-8") as fp:
        for obj in objects:
            labels = [
                shrink_label(poly, iterations=args.iterations, seed=args.seed + idx, unknot_threshold=args.threshold)
                for idx, poly in enumerate(obj.polygons)
            ]
            fp.write(" ".join(labels) + "\n")


def cmd_shrink_lengths(args):
    objects = read_objects(args.input, args.polys_per_object)
    for obj_index, obj in enumerate(objects, start=1):
        lengths = [
            shrink_length(poly, iterations=args.iterations, seed=args.seed + idx)
            for idx, poly in enumerate(obj.polygons)
        ]
        print(obj_index, *lengths)


def cmd_tally(args):
    counts: Counter[str] = Counter()
    for path in args.inputs:
        with Path(path).open("r", encoding="utf-8") as fp:
            for line in fp:
                label = line.strip()
                if label:
                    counts[label] += 1
    for label, count in sorted(counts.items()):
        print(f"{label} {count}")


def cmd_linking_number(args):
    objects = read_objects(args.input, 2)
    if args.jobs < 1:
        raise SystemExit("--jobs must be at least 1")
    jobs = min(args.jobs, max(1, len(objects)))
    work = zip(range(1, len(objects) + 1), objects, itertools.repeat(args.keep_going))
    if jobs == 1 or len(objects) <= 1:
        lines = [_linking_number_line(item) for item in work]
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=jobs) as executor:
            lines = list(executor.map(_linking_number_line, work))
    for line in lines:
        print(line)


def cmd_split_by_label(args):
    objects = read_objects(args.input, args.polys_per_object)
    labels = [line.strip() for line in Path(args.labels).read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(labels) != len(objects):
        raise SystemExit(f"label count {len(labels)} does not match object count {len(objects)}")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    buckets = {}
    for obj, label in zip(objects, labels):
        buckets.setdefault(label, []).append(obj)
    from postprocess.uofs import write_objects

    for label, bucket in buckets.items():
        safe = label.replace("/", "_").replace(" ", "_")
        write_objects(output_dir / f"{safe}.txt", bucket)
        print(f"{label} {len(bucket)}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("shrink-id", help="Run BFACF-style shrink knot/unknot labeling")
    p.add_argument("input")
    p.add_argument("output")
    p.add_argument("-p", "--polys-per-object", type=int, default=1)
    p.add_argument("--iterations", type=int)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--threshold", type=int, default=24)
    p.set_defaults(func=cmd_shrink_id)

    p = sub.add_parser("shrink-lengths", help="Print shrunk lengths for each object")
    p.add_argument("input")
    p.add_argument("-p", "--polys-per-object", type=int, default=1)
    p.add_argument("--iterations", type=int)
    p.add_argument("--seed", type=int, default=1)
    p.set_defaults(func=cmd_shrink_lengths)

    p = sub.add_parser("tally-labels", help="Tally knot/link label files")
    p.add_argument("inputs", nargs="+")
    p.set_defaults(func=cmd_tally)

    p = sub.add_parser("linking-number", help="Compute projection-based linking numbers for 2SAP UofS objects")
    p.add_argument("input")
    p.add_argument("--keep-going", action="store_true")
    p.add_argument("-j", "--jobs", type=int, default=1, help="Worker processes for per-object linking numbers")
    p.set_defaults(func=cmd_linking_number)

    p = sub.add_parser("split-by-label", help="Split a UofS file into one file per label")
    p.add_argument("input")
    p.add_argument("labels")
    p.add_argument("output_dir")
    p.add_argument("-p", "--polys-per-object", type=int, default=1)
    p.set_defaults(func=cmd_split_by_label)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
