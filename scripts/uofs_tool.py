#!/usr/bin/env python3
"""Modern UofS post-processing command line tools."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from postprocess.analysis import aggregate_contacts, edge_count_distribution, span_distribution
from postprocess.uofs import (
    canonical_object,
    read_objects,
    validate_object,
    write_lp,
    write_objects,
    write_points,
    write_printable_dir,
)


def load(args):
    return read_objects(args.input, args.polys_per_object)


def cmd_summary(args):
    objects = load(args)
    poly_count = sum(len(obj.polygons) for obj in objects)
    print(f"objects: {len(objects)}")
    print(f"polygons: {poly_count}")
    print(f"total_edges: {sum(obj.total_length for obj in objects)}")
    if objects:
        print(f"average_edges_per_object: {sum(obj.total_length for obj in objects) / len(objects):.6f}")
        print(f"average_{args.span_axis}_span: {sum(obj.axis_span(args.span_axis) for obj in objects) / len(objects):.6f}")


def cmd_to_points(args):
    write_points(args.output, load(args))


def cmd_to_printable(args):
    write_printable_dir(args.output_dir, load(args))


def cmd_to_lp(args):
    write_lp(args.output, load(args))


def cmd_validate(args):
    objects = load(args)
    failures = 0
    for idx, obj in enumerate(objects, start=1):
        errors = validate_object(obj, args.L, args.M, args.span, args.hamiltonian, args.span_axis)
        if errors:
            failures += 1
            print(f"object {idx}: FAIL")
            for error in errors:
                print(f"  {error}")
    print(f"validated_objects: {len(objects)}")
    print(f"failed_objects: {failures}")
    if failures:
        raise SystemExit(1)


def cmd_count_edges(args):
    dist = edge_count_distribution(load(args))
    for length, count in sorted(dist.items()):
        print(f"{length} {count}")


def cmd_count_spans(args):
    dist = span_distribution(load(args), args.span_axis)
    for span, count in sorted(dist.items()):
        print(f"{span} {count}")


def cmd_canonicalize(args):
    objects = [canonical_object(obj, unordered=args.unordered) for obj in load(args)]
    write_objects(args.output, objects)


def cmd_unique(args):
    seen = {}
    unique = []
    for obj in load(args):
        canon = canonical_object(obj, unordered=args.unordered)
        key = tuple((poly.start, poly.directions) for poly in canon.polygons)
        seen[key] = seen.get(key, 0) + 1
        if seen[key] == 1:
            unique.append(canon)
    if args.output:
        write_objects(args.output, unique)
    print(f"unique_objects: {len(unique)}")
    print(f"duplicate_objects: {sum(seen.values()) - len(unique)}")


def cmd_contacts(args):
    contacts = aggregate_contacts(load(args), args.radius, jobs=args.jobs)
    if args.output:
        with Path(args.output).open("w", newline="", encoding="utf-8") as fp:
            writer = csv.writer(fp)
            writer.writerow(["cyclic_distance", "contacts"])
            for distance, count in sorted(contacts.items()):
                writer.writerow([distance, count])
    else:
        for distance, count in sorted(contacts.items()):
            print(f"{distance} {count}")


def add_common(parser):
    parser.add_argument("input", help="Input UofS file")
    parser.add_argument("-p", "--polys-per-object", type=int, default=1, help="1 for SAP, 2 for 2SAP")
    parser.add_argument("--span-axis", choices=["x", "y", "z"], default="x", help="Coordinate axis used as tube span")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("summary", help="Print object, polygon, edge, and span summaries")
    add_common(p)
    p.set_defaults(func=cmd_summary)

    p = sub.add_parser("to-points", help="Convert UofS directions to point-coordinate format")
    add_common(p)
    p.add_argument("output")
    p.set_defaults(func=cmd_to_points)

    p = sub.add_parser("to-printable", help="Write one coordinate file per polygon")
    add_common(p)
    p.add_argument("output_dir")
    p.set_defaults(func=cmd_to_printable)

    p = sub.add_parser("to-lp", help="Convert UofS to length-points format")
    add_common(p)
    p.add_argument("output")
    p.set_defaults(func=cmd_to_lp)

    p = sub.add_parser("validate", help="Validate closure, self-avoidance, bounds, and Hamiltonian coverage")
    add_common(p)
    p.add_argument("-L", type=int)
    p.add_argument("-M", type=int)
    p.add_argument("-s", "--span", type=int)
    p.add_argument("--hamiltonian", action="store_true")
    p.set_defaults(func=cmd_validate)

    p = sub.add_parser("count-edges", help="Count objects by total edge length")
    add_common(p)
    p.set_defaults(func=cmd_count_edges)

    p = sub.add_parser("count-spans", help="Count objects by z-span")
    add_common(p)
    p.set_defaults(func=cmd_count_spans)

    p = sub.add_parser("canonicalize", help="Write canonicalized UofS objects")
    add_common(p)
    p.add_argument("output")
    p.add_argument("--unordered", action="store_true", help="Sort polygons inside each object, useful for 2SAP")
    p.set_defaults(func=cmd_canonicalize)

    p = sub.add_parser("unique", help="Count and optionally write unique canonical objects")
    add_common(p)
    p.add_argument("-o", "--output")
    p.add_argument("--unordered", action="store_true", help="Sort polygons inside each object, useful for 2SAP")
    p.set_defaults(func=cmd_unique)

    p = sub.add_parser("contacts", help="Compute contact counts by cyclic distance")
    add_common(p)
    p.add_argument("-r", "--radius", type=float, default=1.0)
    p.add_argument("-o", "--output")
    p.add_argument("-j", "--jobs", type=int, default=1, help="Worker processes for per-object contact counts")
    p.set_defaults(func=cmd_contacts)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
