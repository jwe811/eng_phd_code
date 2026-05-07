from __future__ import annotations

from postprocess.uofs import read_objects, validate_object


def infer_polys_per_object(path: str) -> int:
    lowered = path.lower()
    return 2 if "2sap" in lowered else 1


def summarize_uofs(path, polys_per_object: int | None = None) -> dict:
    count = polys_per_object or infer_polys_per_object(str(path))
    objects = read_objects(path, count)
    total_edges = sum(obj.total_length for obj in objects)
    return {
        "path": str(path),
        "objects": len(objects),
        "polygons": sum(len(obj.polygons) for obj in objects),
        "total_edges": total_edges,
        "average_edges_per_object": total_edges / len(objects) if objects else 0,
        "average_x_span": sum(obj.axis_span("x") for obj in objects) / len(objects) if objects else 0,
        "polys_per_object": count,
    }


def object_points(path, index: int, polys_per_object: int | None = None) -> dict:
    count = polys_per_object or infer_polys_per_object(str(path))
    objects = read_objects(path, count)
    if index < 0 or index >= len(objects):
        raise IndexError(index)
    obj = objects[index]
    errors = validate_object(obj, span_axis="x")
    return {
        "path": str(path),
        "index": index,
        "object_count": len(objects),
        "validation_errors": errors,
        "polygons": [
            {
                "start": poly.start,
                "length": poly.length,
                "points": poly.points(),
            }
            for poly in obj.polygons
        ],
    }

