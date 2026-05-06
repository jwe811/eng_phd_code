"""File-based SAP/2SAP post-processing analyses."""

from __future__ import annotations

import concurrent.futures
import itertools
import os
from collections import Counter

from .uofs import UofSObject


def edge_count_distribution(objects: list[UofSObject]) -> Counter[int]:
    return Counter(obj.total_length for obj in objects)


def span_distribution(objects: list[UofSObject], axis: str = "x") -> Counter[int]:
    return Counter(obj.axis_span(axis) for obj in objects)


def contact_distribution(obj: UofSObject, radius: float = 1.0) -> Counter[int]:
    counts: Counter[int] = Counter()
    r2 = radius * radius
    for poly in obj.polygons:
        pts = poly.points()
        if pts[-1] == pts[0]:
            pts = pts[:-1]
        n = len(pts)
        for i in range(n):
            for j in range(i + 1, n):
                cyclic_distance = min(j - i, n - (j - i))
                if cyclic_distance <= 1:
                    continue
                dx = pts[i][0] - pts[j][0]
                dy = pts[i][1] - pts[j][1]
                dz = pts[i][2] - pts[j][2]
                if dx * dx + dy * dy + dz * dz <= r2:
                    counts[cyclic_distance] += 1
    return counts


def _contact_distribution_star(args: tuple[UofSObject, float]) -> Counter[int]:
    obj, radius = args
    return contact_distribution(obj, radius)


def resolve_jobs(jobs: int | None, item_count: int) -> int:
    if jobs is None:
        jobs = os.cpu_count() or 1
    if jobs < 1:
        raise ValueError("jobs must be at least 1")
    return min(jobs, max(1, item_count))


def aggregate_contacts(objects: list[UofSObject], radius: float = 1.0, jobs: int = 1) -> Counter[int]:
    total: Counter[int] = Counter()
    worker_count = resolve_jobs(jobs, len(objects))
    if worker_count == 1 or len(objects) <= 1:
        for obj in objects:
            total.update(contact_distribution(obj, radius))
        return total

    with concurrent.futures.ProcessPoolExecutor(max_workers=worker_count) as executor:
        for counts in executor.map(_contact_distribution_star, zip(objects, itertools.repeat(radius))):
            total.update(counts)
    return total
