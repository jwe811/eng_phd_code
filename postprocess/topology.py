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


def _cross3(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _dot(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _normalize(vec: tuple[float, float, float]) -> tuple[float, float, float]:
    norm = (vec[0] ** 2 + vec[1] ** 2 + vec[2] ** 2) ** 0.5
    if norm == 0.0:
        raise ValueError("Zero projection normal")
    return (vec[0] / norm, vec[1] / norm, vec[2] / norm)


def _build_basis(normal: tuple[float, float, float]) -> tuple[tuple[float, float, float], ...]:
    w = _normalize(normal)
    if abs(w[0]) < 0.9:
        tangent = (1.0, 0.0, 0.0)
    else:
        tangent = (0.0, 1.0, 0.0)
    u = _normalize(_cross3(w, tangent))
    v = _cross3(w, u)
    return u, v, w


def _project(pt: tuple[int, int, int], basis: tuple[tuple[float, float, float], ...]) -> tuple[float, float, float]:
    u, v, w = basis
    return (_dot(pt, u), _dot(pt, v), _dot(pt, w))


def _linking_projection(
    poly_a: Polygon,
    poly_b: Polygon,
    basis: tuple[tuple[float, float, float], ...],
    eps: float = 1e-9,
) -> float | None:
    total = 0
    for a0, a1 in _segments(poly_a):
        ax0, ay0, ah0 = _project(a0, basis)
        ax1, ay1, ah1 = _project(a1, basis)
        adx, ady = ax1 - ax0, ay1 - ay0
        for b0, b1 in _segments(poly_b):
            bx0, by0, bh0 = _project(b0, basis)
            bx1, by1, bh1 = _project(b1, basis)
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
    candidate_normals = [
        (0.0, 0.0, 1.0),
        (0.0, 1.0, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, 1.0, 0.0),
        (1.0, 0.0, 1.0),
        (0.0, 1.0, 1.0),
        (1.0, 1.0, 1.0),
        (1.0, 2.0, 0.0),
        (2.0, 1.0, 0.0),
        (1.0, 0.0, 2.0),
        (0.0, 1.0, 2.0),
        (2.0, 0.0, 1.0),
        (0.0, 2.0, 1.0),
        (1.0, 2.0, 3.0),
        (1.0, 3.0, 2.0),
        (2.0, 3.0, 5.0),
        (3.0, 4.0, 5.0),
        (3.0, 5.0, 7.0),
        (5.0, 7.0, 11.0),
        (2.0, 5.0, 7.0),
    ]
    for normal in candidate_normals:
        basis = _build_basis(normal)
        value = _linking_projection(poly_a, poly_b, basis)
        if value is not None:
            return value
    raise ValueError("Could not find a nondegenerate projection for linking number")

