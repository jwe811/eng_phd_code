"""Topology-oriented post-processing helpers."""

from __future__ import annotations

from .uofs import Polygon


def _segments(poly: Polygon):
    pts = poly.points()
    for a, b in zip(pts, pts[1:]):
        yield a, b


def _project(pt, drop_axis: str):
    axes = {"x": (1, 2, 0), "y": (0, 2, 1), "z": (0, 1, 2)}[drop_axis]
    return (pt[axes[0]], pt[axes[1]], pt[axes[2]])


def _cross2(ax, ay, bx, by):
    return ax * by - ay * bx


def _linking_projection(poly_a: Polygon, poly_b: Polygon, drop_axis: str, eps: float = 1e-9) -> float | None:
    total = 0
    for a0, a1 in _segments(poly_a):
        ax0, ay0, ah0 = _project(a0, drop_axis)
        ax1, ay1, ah1 = _project(a1, drop_axis)
        adx, ady = ax1 - ax0, ay1 - ay0
        for b0, b1 in _segments(poly_b):
            bx0, by0, bh0 = _project(b0, drop_axis)
            bx1, by1, bh1 = _project(b1, drop_axis)
            bdx, bdy = bx1 - bx0, by1 - by0
            denom = _cross2(adx, ady, bdx, bdy)
            if abs(denom) < eps:
                continue
            qpx, qpy = bx0 - ax0, by0 - ay0
            t = _cross2(qpx, qpy, bdx, bdy) / denom
            u = _cross2(qpx, qpy, adx, ady) / denom
            if eps < t < 1.0 - eps and eps < u < 1.0 - eps:
                ah = ah0 + t * (ah1 - ah0)
                bh = bh0 + u * (bh1 - bh0)
                if abs(ah - bh) < eps:
                    return None
                total += (1 if denom > 0 else -1) * (1 if bh > ah else -1)
            elif -eps <= t <= 1.0 + eps and -eps <= u <= 1.0 + eps:
                return None
    return total / 2.0


def linking_number(poly_a: Polygon, poly_b: Polygon) -> float:
    for axis in ("z", "y", "x"):
        value = _linking_projection(poly_a, poly_b, axis)
        if value is not None:
            return value
    raise ValueError("Could not find a nondegenerate coordinate projection for linking number")

