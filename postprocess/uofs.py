"""UofS polygon file parsing, writing, and geometry helpers.

The Monte Carlo and CreatorAll programs write one polygon as a starting
coordinate followed by direction codes, terminated by -111. A file terminates
with -999. For 2SAP files, two consecutive polygons form one object.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence


DIR_DELTA = {
    1: (1, 0, 0),
    2: (-1, 0, 0),
    3: (0, 1, 0),
    4: (0, -1, 0),
    5: (0, 0, 1),
    6: (0, 0, -1),
}
OPPOSITE = {1: 2, 2: 1, 3: 4, 4: 3, 5: 6, 6: 5}


@dataclass(frozen=True)
class Polygon:
    start: tuple[int, int, int]
    directions: tuple[int, ...]

    def points(self) -> list[tuple[int, int, int]]:
        x, y, z = self.start
        pts = [(x, y, z)]
        for direction in self.directions:
            if direction not in DIR_DELTA:
                raise ValueError(f"Invalid direction code {direction}")
            dx, dy, dz = DIR_DELTA[direction]
            x += dx
            y += dy
            z += dz
            pts.append((x, y, z))
        return pts

    def translated(self, delta: tuple[int, int, int]) -> "Polygon":
        return Polygon(tuple(self.start[i] + delta[i] for i in range(3)), self.directions)

    @property
    def length(self) -> int:
        return len(self.directions)

    @property
    def z_span(self) -> int:
        pts = self.points()
        zs = [p[2] for p in pts]
        return max(zs) - min(zs)

    def axis_span(self, axis: str = "x") -> int:
        idx = {"x": 0, "y": 1, "z": 2}[axis]
        coords = [p[idx] for p in self.points()]
        return max(coords) - min(coords)


@dataclass(frozen=True)
class UofSObject:
    polygons: tuple[Polygon, ...]

    @property
    def total_length(self) -> int:
        return sum(poly.length for poly in self.polygons)

    @property
    def z_span(self) -> int:
        points = [pt for poly in self.polygons for pt in poly.points()]
        zs = [p[2] for p in points]
        return max(zs) - min(zs)

    def axis_span(self, axis: str = "x") -> int:
        idx = {"x": 0, "y": 1, "z": 2}[axis]
        points = [pt for poly in self.polygons for pt in poly.points()]
        coords = [p[idx] for p in points]
        return max(coords) - min(coords)


def read_polygons(path: str | Path) -> list[Polygon]:
    polygons: list[Polygon] = []
    with Path(path).open("r", encoding="utf-8") as fp:
        first = fp.readline()
        if first and first.strip() != "UofS":
            fp.seek(0)
        while True:
            line = fp.readline()
            if not line:
                break
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if parts[0] == "-999":
                break
            if len(parts) != 3:
                raise ValueError(f"Expected polygon start coordinate in {path}, got: {line}")
            start = tuple(int(part) for part in parts)
            directions: list[int] = []
            while True:
                direction_line = fp.readline()
                if not direction_line:
                    raise ValueError(f"Unterminated polygon in {path}")
                direction = int(direction_line.strip())
                if direction == -111:
                    break
                if direction not in DIR_DELTA:
                    raise ValueError(f"Invalid direction code {direction} in {path}")
                directions.append(direction)
            polygons.append(Polygon(start, tuple(directions)))
    return polygons


def group_objects(polygons: Sequence[Polygon], polygons_per_object: int) -> list[UofSObject]:
    if polygons_per_object <= 0:
        raise ValueError("polygons_per_object must be positive")
    if len(polygons) % polygons_per_object != 0:
        raise ValueError(
            f"File contains {len(polygons)} polygons, not divisible by {polygons_per_object}"
        )
    return [
        UofSObject(tuple(polygons[i : i + polygons_per_object]))
        for i in range(0, len(polygons), polygons_per_object)
    ]


def read_objects(path: str | Path, polygons_per_object: int = 1) -> list[UofSObject]:
    return group_objects(read_polygons(path), polygons_per_object)


def write_objects(path: str | Path, objects: Iterable[UofSObject]) -> None:
    with Path(path).open("w", encoding="utf-8") as fp:
        fp.write("UofS\n")
        for obj in objects:
            for poly in obj.polygons:
                fp.write(f"{poly.start[0]} {poly.start[1]} {poly.start[2]}\n")
                for direction in poly.directions:
                    fp.write(f"{direction}\n")
                fp.write("-111\n")
        fp.write("-999\n")


def write_points(path: str | Path, objects: Iterable[UofSObject]) -> None:
    with Path(path).open("w", encoding="utf-8") as fp:
        for obj in objects:
            for poly in obj.polygons:
                for x, y, z in poly.points():
                    fp.write(f"{x} {y} {z}\n")
                fp.write("-111\n")
        fp.write("-999\n")


def write_printable_dir(output_dir: str | Path, objects: Iterable[UofSObject]) -> None:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    idx = 1
    for obj in objects:
        for poly in obj.polygons:
            write_points(out / f"poly_{idx}.txt", [UofSObject((poly,))])
            idx += 1


def write_lp(path: str | Path, objects: Iterable[UofSObject]) -> None:
    with Path(path).open("w", encoding="utf-8") as fp:
        for obj in objects:
            for poly in obj.polygons:
                pts = poly.points()
                if pts[-1] == pts[0]:
                    pts = pts[:-1]
                fp.write(str(poly.length))
                for x, y, z in pts:
                    fp.write(f" {x} {y} {z}")
                fp.write("\n")


def iter_lp_polygons(path: str | Path) -> Iterator[list[tuple[int, int, int]]]:
    with Path(path).open("r", encoding="utf-8") as fp:
        for line in fp:
            parts = [int(part) for part in line.split()]
            if not parts:
                continue
            length = parts[0]
            coords = parts[1:]
            if len(coords) % 3 != 0:
                raise ValueError(f"LP line has non-triplet coordinates: {line.strip()}")
            pts = [(coords[i], coords[i + 1], coords[i + 2]) for i in range(0, len(coords), 3)]
            if len(pts) not in (length, length + 1):
                raise ValueError(f"LP length {length} is inconsistent with {len(pts)} points")
            yield pts


def _direction_between(a: tuple[int, int, int], b: tuple[int, int, int]) -> int:
    delta = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
    for direction, known in DIR_DELTA.items():
        if delta == known:
            return direction
    raise ValueError(f"Points {a} and {b} are not unit lattice neighbors")


def canonical_polygon(poly: Polygon) -> Polygon:
    pts = poly.points()
    if pts[-1] != pts[0]:
        raise ValueError("Cannot canonicalize an open polygon")
    cycle = pts[:-1]
    variants: list[Polygon] = []
    n = len(cycle)
    for reverse in (False, True):
        ordered = list(reversed(cycle)) if reverse else list(cycle)
        for offset in range(n):
            rotated = ordered[offset:] + ordered[:offset]
            min_x = min(p[0] for p in rotated)
            min_y = min(p[1] for p in rotated)
            min_z = min(p[2] for p in rotated)
            shifted = [(p[0] - min_x, p[1] - min_y, p[2] - min_z) for p in rotated]
            closed = shifted + [shifted[0]]
            directions = tuple(_direction_between(closed[i], closed[i + 1]) for i in range(n))
            variants.append(Polygon(shifted[0], directions))
    return min(variants, key=lambda p: (p.start, p.directions))


def _canonical_polygon_fixed_space(poly: Polygon, delta: tuple[int, int, int]) -> Polygon:
    pts = poly.points()
    if pts[-1] != pts[0]:
        raise ValueError("Cannot canonicalize an open polygon")
    cycle = pts[:-1]
    variants: list[Polygon] = []
    n = len(cycle)
    for reverse in (False, True):
        ordered = list(reversed(cycle)) if reverse else list(cycle)
        for offset in range(n):
            rotated = ordered[offset:] + ordered[:offset]
            shifted = [(p[0] + delta[0], p[1] + delta[1], p[2] + delta[2]) for p in rotated]
            closed = shifted + [shifted[0]]
            directions = tuple(_direction_between(closed[i], closed[i + 1]) for i in range(n))
            variants.append(Polygon(shifted[0], directions))
    return min(variants, key=lambda p: (p.start, p.directions))


def canonical_object(obj: UofSObject, unordered: bool = False) -> UofSObject:
    if len(obj.polygons) == 1:
        polys = (canonical_polygon(obj.polygons[0]),)
    else:
        all_points = [pt for poly in obj.polygons for pt in poly.points()[:-1]]
        min_x = min(p[0] for p in all_points)
        min_y = min(p[1] for p in all_points)
        min_z = min(p[2] for p in all_points)
        delta = (-min_x, -min_y, -min_z)
        polys = tuple(_canonical_polygon_fixed_space(poly, delta) for poly in obj.polygons)
    if unordered:
        polys = tuple(sorted(polys, key=lambda p: (p.start, p.directions)))
    return UofSObject(polys)


def validate_object(
    obj: UofSObject,
    lattice_l: int | None = None,
    lattice_m: int | None = None,
    span: int | None = None,
    hamiltonian: bool = False,
    span_axis: str = "x",
) -> list[str]:
    errors: list[str] = []
    seen_vertices: set[tuple[int, int, int]] = set()
    for poly_idx, poly in enumerate(obj.polygons, start=1):
        try:
            pts = poly.points()
        except ValueError as exc:
            return [str(exc)]
        if pts[-1] != pts[0]:
            errors.append(f"polygon {poly_idx} is not closed")
        interior = pts[:-1] if pts[-1] == pts[0] else pts
        if len(set(interior)) != len(interior):
            errors.append(f"polygon {poly_idx} repeats a vertex before closure")
        if span_axis not in {"x", "y", "z"}:
            raise ValueError("span_axis must be x, y, or z")
        axis_limits = {"x": None, "y": None, "z": None}
        axis_limits[span_axis] = span
        transverse = [axis for axis in ("x", "y", "z") if axis != span_axis]
        axis_limits[transverse[0]] = lattice_l
        axis_limits[transverse[1]] = lattice_m
        for x, y, z in interior:
            values = {"x": x, "y": y, "z": z}
            for axis, limit in axis_limits.items():
                if limit is not None and not (0 <= values[axis] <= limit):
                    errors.append(f"polygon {poly_idx} has {axis}={values[axis]} outside [0,{limit}]")
            if x == -999:
                errors.append("unexpected sentinel-like coordinate in polygon")
        overlap = seen_vertices.intersection(interior)
        if overlap:
            errors.append(f"polygons overlap at {min(overlap)}")
        seen_vertices.update(interior)
    if hamiltonian:
        if lattice_l is None or lattice_m is None or span is None:
            errors.append("Hamiltonian validation requires L, M, and span")
        else:
            expected = (lattice_l + 1) * (lattice_m + 1) * (span + 1)
            if len(seen_vertices) != expected:
                errors.append(f"Hamiltonian object visits {len(seen_vertices)} vertices, expected {expected}")
    return errors
