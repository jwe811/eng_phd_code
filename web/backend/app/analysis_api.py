from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

from postprocess.analysis import aggregate_contacts, edge_count_distribution, span_distribution
from postprocess.bfacf import shrink_label
from postprocess.topology import linking_number
from postprocess.uofs import read_objects, validate_object

from .uofs_api import infer_polys_per_object, summarize_uofs


FILENAME_GEOMETRY = re.compile(r"L(?P<L>\d+)M(?P<M>\d+)span(?P<span>\d+)")


def infer_geometry(path: Path) -> tuple[int | None, int | None, int | None]:
    match = FILENAME_GEOMETRY.search(path.name)
    if not match:
        return None, None, None
    return int(match.group("L")), int(match.group("M")), int(match.group("span"))


def _counter_rows(counter: Counter[int]) -> list[dict[str, int]]:
    return [{"value": value, "count": count} for value, count in sorted(counter.items())]


def run_analysis(path: Path, action: str) -> dict:
    polys_per_object = infer_polys_per_object(str(path))
    objects = read_objects(path, polys_per_object)

    if action == "summary":
        return {"action": action, "path": str(path), "metrics": summarize_uofs(path, polys_per_object), "rows": []}

    if action == "validate":
        lattice_l, lattice_m, span = infer_geometry(path)
        hamiltonian = "ham" in str(path).lower()
        failures = []
        for idx, obj in enumerate(objects, start=1):
            errors = validate_object(obj, lattice_l, lattice_m, span, hamiltonian=hamiltonian, span_axis="x")
            if errors:
                failures.append({"object": idx, "errors": errors[:5]})
        return {
            "action": action,
            "path": str(path),
            "metrics": {
                "validated_objects": len(objects),
                "failed_objects": len(failures),
                "inferred_L": lattice_l,
                "inferred_M": lattice_m,
                "inferred_span": span,
                "hamiltonian": hamiltonian,
            },
            "rows": failures[:50],
        }

    if action == "count_edges":
        return {"action": action, "path": str(path), "metrics": {"objects": len(objects)}, "rows": _counter_rows(edge_count_distribution(objects))}

    if action == "count_spans":
        return {"action": action, "path": str(path), "metrics": {"objects": len(objects), "axis": "x"}, "rows": _counter_rows(span_distribution(objects, "x"))}

    if action == "contacts":
        contacts = aggregate_contacts(objects, radius=1.0, jobs=1)
        return {"action": action, "path": str(path), "metrics": {"objects": len(objects), "radius": 1.0}, "rows": _counter_rows(contacts)}

    if action == "linking_number":
        if polys_per_object != 2:
            raise ValueError("linking number requires a 2SAP file")
        rows = []
        for idx, obj in enumerate(objects[:200], start=1):
            try:
                rows.append({"object": idx, "linking_number": linking_number(obj.polygons[0], obj.polygons[1])})
            except ValueError as exc:
                rows.append({"object": idx, "error": str(exc)})
        return {"action": action, "path": str(path), "metrics": {"objects": len(objects), "reported_objects": len(rows)}, "rows": rows}

    if action == "shrink_labels":
        rows = []
        for obj_idx, obj in enumerate(objects[:50], start=1):
            labels = [
                shrink_label(poly, seed=1 + poly_idx)
                for poly_idx, poly in enumerate(obj.polygons)
            ]
            rows.append({"object": obj_idx, "labels": labels})
        return {"action": action, "path": str(path), "metrics": {"objects": len(objects), "reported_objects": len(rows)}, "rows": rows}

    raise ValueError(f"unknown analysis action: {action}")

