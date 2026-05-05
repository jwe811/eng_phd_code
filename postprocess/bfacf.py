"""BFACF-style shrink heuristic revived from the archive.

This is intended as a conservative post-processing classifier. It performs the
archive's core -2 shrink move and corner-swap move with a deterministic random
seed option. It reports "unknot" when a polygon shrinks below the traditional
24-edge threshold; otherwise it reports "knot_or_unresolved".
"""

from __future__ import annotations

import random

from .uofs import OPPOSITE, Polygon


def _can_shrink(a: int, b: int, c: int) -> bool:
    return OPPOSITE.get(a) == c and b != a and b != c


def _is_right_angle(a: int, b: int) -> bool:
    return a != b and OPPOSITE.get(a) != b


def shrink_length(poly: Polygon, iterations: int | None = None, seed: int = 1) -> int:
    rng = random.Random(seed)
    directions = list(poly.directions)
    if iterations is None:
        iterations = 300 * max(1, len(directions))
    for _ in range(iterations):
        n = len(directions)
        if n < 24:
            break
        chosen = rng.randrange(n)
        one = (chosen - 1) % n
        two = chosen
        three = (chosen + 1) % n
        if _can_shrink(directions[one], directions[two], directions[three]):
            for idx in sorted((one, three), reverse=True):
                directions.pop(idx)
            continue
        nxt = (chosen + 1) % len(directions)
        if _is_right_angle(directions[chosen], directions[nxt]):
            directions[chosen], directions[nxt] = directions[nxt], directions[chosen]
    return len(directions)


def shrink_label(poly: Polygon, iterations: int | None = None, seed: int = 1, unknot_threshold: int = 24) -> str:
    return "unknot" if shrink_length(poly, iterations, seed) < unknot_threshold else "knot_or_unresolved"

