"""Remark 2.2 of the note: attach the label of a step to the two swapped VALUES
instead of the two swapped positions, and see where the failure of SCNP goes."""
from functools import lru_cache
from itertools import combinations


def covers(v):
    n = len(v)
    for a in range(n):
        for b in range(a + 1, n):
            between = any(v[a] < v[c] < v[b] for c in range(a + 1, b))
            if v[a] < v[b] and not between:
                t = list(v)
                t[a], t[b] = t[b], t[a]
                yield tuple(t), (a + 1, b + 1), (v[a], v[b])


@lru_cache(None)
def reaches(v, w):
    return v == w or any(reaches(y, w) for y, _, _ in covers(v))


def chains(u, w, by_values):
    if u == w:
        yield []
        return
    for v, positions, values in covers(u):
        if reaches(v, w):
            for rest in chains(v, w, by_values):
                yield [values if by_values else positions] + rest


def support(steps, n):
    exps = {(0,) * (n - 1)}
    for a, b in steps:
        exps = {e[:k - 1] + (e[k - 1] + 1,) + e[k:]
                for e in exps for k in range(a, b)}
    return exps


u = (2, 1, 4, 3, 6, 5)
for w in [(5, 4, 1, 6, 3, 2), (3, 6, 5, 2, 1, 4)]:
    name = "".join(map(str, w))
    print(name, "contains 4231:",
          any(d < b < c < a for a, b, c, d in combinations(w, 4)))
    for by_values in (False, True):
        supports = [support(c, 6) for c in chains(u, w, by_values)]
        total = set().union(*supports)
        print("  214365 <=", name,
              "labels on values:   " if by_values else "labels on positions:",
              "monomials", len(total),
              "| best chain", max(map(len, supports)),
              "| SCNP", any(s == total for s in supports))
