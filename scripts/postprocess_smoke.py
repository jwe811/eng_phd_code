#!/usr/bin/env python3
"""Smoke tests for revived post-processing modules."""

from __future__ import annotations

import tempfile
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from postprocess.analysis import aggregate_contacts, edge_count_distribution
from postprocess.bfacf import shrink_label
from postprocess.topology import linking_number
from postprocess.uofs import (
    Polygon,
    canonical_object,
    read_objects,
    validate_object,
    write_lp,
    write_objects,
    write_points,
)


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        source = root / "square.txt"
        source.write_text("UofS\n0 0 0\n1\n3\n2\n4\n-111\n-999\n", encoding="utf-8")

        objects = read_objects(source)
        assert len(objects) == 1
        assert objects[0].total_length == 4
        assert validate_object(objects[0], 1, 1, 0, hamiltonian=True, span_axis="z") == []
        assert edge_count_distribution(objects)[4] == 1
        assert shrink_label(objects[0].polygons[0]) == "unknot"

        points = root / "points.txt"
        write_points(points, objects)
        assert points.read_text(encoding="utf-8").startswith("0 0 0\n1 0 0\n")

        lp = root / "square_lp.txt"
        write_lp(lp, objects)
        assert lp.read_text(encoding="utf-8").strip() == "4 0 0 0 1 0 0 1 1 0 0 1 0"

        canon = canonical_object(objects[0])
        out = root / "canon.txt"
        write_objects(out, [canon])
        assert len(read_objects(out)) == 1
        assert sum(aggregate_contacts(objects).values()) == 0
        assert linking_number(objects[0].polygons[0], objects[0].polygons[0]) == 0

        # Regression test for a 2SAP pair whose coordinate projections are degenerate
        poly_a = Polygon((2, 0, 0), (2, 2, 3, 1, 1, 5, 4, 6))
        poly_b = Polygon((2, 2, 0), (2, 2, 5, 4, 4, 1, 3, 3, 1, 6))
        assert linking_number(poly_a, poly_b) == 0

    print("postprocess smoke tests passed")


if __name__ == "__main__":
    main()
