import re
from postprocess.uofs import read_objects, validate_object
from postprocess.analysis import contact_distribution
from postprocess.bfacf import shrink_label
from postprocess.topology import linking_number


def infer_polys_per_object(path: str) -> int:
    lowered = path.lower()
    return 2 if "2sap" in lowered else 1


def infer_lattice_params(path: str) -> dict:
    match = re.search(r"L(\d+)M(\d+)span(\d+)", path)
    if match:
        return {
            "L": int(match.group(1)),
            "M": int(match.group(2)),
            "span": int(match.group(3)),
        }
    return {"L": None, "M": None, "span": None}


def summarize_uofs(path, polys_per_object: int | None = None) -> dict:
    count = polys_per_object or infer_polys_per_object(str(path))
    objects = read_objects(path, count)
    total_edges = sum(obj.total_length for obj in objects)
    res = {
        "path": str(path),
        "objects": len(objects),
        "polygons": sum(len(obj.polygons) for obj in objects),
        "total_edges": total_edges,
        "average_edges_per_object": total_edges / len(objects) if objects else 0,
        "average_x_span": sum(obj.axis_span("x") for obj in objects) / len(objects) if objects else 0,
        "polys_per_object": count,
    }
    res.update(infer_lattice_params(str(path)))
    return res


def object_points(path, index: int, polys_per_object: int | None = None) -> dict:
    count = polys_per_object or infer_polys_per_object(str(path))
    objects = read_objects(path, count)
    if index < 0 or index >= len(objects):
        raise IndexError(index)
    obj = objects[index]
    params = infer_lattice_params(str(path))
    errors = validate_object(
        obj,
        lattice_l=params["L"],
        lattice_m=params["M"],
        span=params["span"],
        span_axis="x",
    )
    
    # Calculate extra metrics
    contacts = sum(contact_distribution(obj, radius=1.0).values())
    labels = [shrink_label(p) for p in obj.polygons]
    ln = None
    if len(obj.polygons) == 2:
        try:
            ln = linking_number(obj.polygons[0], obj.polygons[1])
        except Exception:
            ln = "error"

    res = {
        "path": str(path),
        "index": index,
        "object_count": len(objects),
        "validation_errors": errors,
        "span_x": obj.axis_span("x"),
        "contacts": contacts,
        "shrink_labels": labels,
        "linking_number": ln,
        "polygons": [
            {
                "start": poly.start,
                "length": poly.length,
                "points": poly.points(),
            }
            for poly in obj.polygons
        ],
    }
    res.update(params)
    return res

