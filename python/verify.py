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
                yield tuple(t), (a + 1, b + 1)


@lru_cache(None)
def reaches(v, w):
    return v == w or any(reaches(y, w) for y, _ in covers(v))


def chains(u, w):
    if u == w:
        yield []
        return
    for v, step in covers(u):
        if reaches(v, w):
            for rest in chains(v, w):
                yield [step] + rest


def support(steps, n):
    exps = {(0,) * (n - 1)}
    for a, b in steps:
        exps = {e[:k - 1] + (e[k - 1] + 1,) + e[k:]
                for e in exps for k in range(a, b)}
    return exps


u, w = (2, 1, 4, 3, 6, 5), (5, 4, 1, 6, 3, 2)
alpha, beta = (0, 2, 0, 2, 3), (3, 2, 0, 2, 0)
supports = [support(c, len(u)) for c in chains(u, w)]
total = set().union(*supports)
print("w contains 4231:",
      any(d < b < c < a for a, b, c, d in combinations(w, 4)))
print("chains:", len(supports))
print("monomials of D:", len(total))
print("most monomials from one chain:", max(map(len, supports)))
print("alpha, beta in D:", alpha in total, beta in total)
print("chains with both alpha and beta:",
      sum(alpha in s and beta in s for s in supports))
