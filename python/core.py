"""Core definitions, written from scratch for the paper.

Permutations are tuples of values, positions are 1-based in all printed output.
A step of a chain swaps positions a < b (a cover in Bruhat order) and carries the
label x_a + ... + x_{b-1}.  A window is a pair (p, q), 1 <= p < q <= n, and
g_C(p, q) = number of steps (a, b) of C with p <= a and b <= q.
"""
from itertools import permutations, combinations
from functools import lru_cache

def P(s):
    return tuple(int(c, 36) for c in s) if isinstance(s, str) else tuple(s)

def S(v):
    return "".join("0123456789ABCDEFGHIJ"[x] for x in v)

def length(v):
    n = len(v)
    return sum(1 for i in range(n) for j in range(i + 1, n) if v[i] > v[j])

def inv(v):
    n = len(v)
    return frozenset((i + 1, j + 1) for i in range(n) for j in range(i + 1, n) if v[i] > v[j])

def covers_up(v):
    """All (a, b, x) with x = v t_ab covering v; a < b are 1-based positions."""
    n = len(v)
    out = []
    for a in range(n):
        for b in range(a + 1, n):
            if v[a] < v[b] and not any(v[a] < v[c] < v[b] for c in range(a + 1, b)):
                x = list(v); x[a], x[b] = x[b], x[a]
                out.append((a + 1, b + 1, tuple(x)))
    return out

def bruhat_leq(u, w):
    """Tableau criterion: for every k, the sorted first k entries of u are entrywise <= those of w."""
    n = len(u)
    for k in range(1, n):
        a = sorted(u[:k]); b = sorted(w[:k])
        if any(x > y for x, y in zip(a, b)):
            return False
    return True

def windows(n):
    return [(p, q) for p in range(1, n + 1) for q in range(p + 1, n + 1)]

def step_vector(n, a, b):
    return tuple(1 if (p <= a and b <= q) else 0 for (p, q) in windows(n))

def add(x, y):
    return tuple(i + j for i, j in zip(x, y))

def leq(x, y):
    return all(i <= j for i, j in zip(x, y))

def pareto_insert(S_, x):
    for y in S_:
        if leq(y, x):
            return S_
    return [y for y in S_ if not leq(x, y)] + [x]

def interval(u, w):
    """Elements of [u, w] grouped by length, found by going up from u."""
    levels = [{u}]
    seen = {u}
    while True:
        nxt = set()
        for v in levels[-1]:
            for a, b, x in covers_up(v):
                if x not in seen and bruhat_leq(x, w):
                    seen.add(x); nxt.add(x)
        if not nxt:
            break
        levels.append(nxt)
    return levels

def minimal_g(u, w):
    """Pareto-minimal g-vectors of the chains from u to w (list of tuples indexed by windows(n))."""
    n = len(u)
    W = windows(n)
    levels = interval(u, w)
    assert levels[-1] == {w}
    Sv = {w: [tuple([0] * len(W))]}
    for lev in reversed(levels[:-1]):
        for v in lev:
            cur = []
            for a, b, x in covers_up(v):
                if x in Sv:
                    e = step_vector(n, a, b)
                    for g in Sv[x]:
                        cur = pareto_insert(cur, add(e, g))
            Sv[v] = cur
    return Sv[u]

def scnp(u, w):
    return len(minimal_g(u, w)) == 1

# ---------------------------------------------------------------- brute force
def all_chains(u, w):
    """Yield every saturated chain from u to w as a list of steps (a, b)."""
    if u == w:
        yield []
        return
    for a, b, x in covers_up(u):
        if bruhat_leq(x, w):
            for rest in all_chains(x, w):
                yield [(a, b)] + rest

def chain_support(n, steps):
    """Exponent vectors of the product of (x_a + ... + x_{b-1}) over the steps."""
    supp = {tuple([0] * (n - 1))}
    for a, b in steps:
        new = set()
        for e in supp:
            for k in range(a, b):
                f = list(e); f[k - 1] += 1
                new.add(tuple(f))
        supp = new
    return supp

def brute(u, w):
    """(number of chains, |supp D|, largest single-chain support, SCNP?) by expanding every chain."""
    n = len(u)
    seen = {}
    total = set(); nch = 0
    for ch in all_chains(u, w):
        nch += 1
        key = tuple(sorted(ch))
        if key not in seen:
            seen[key] = chain_support(n, ch)
            total |= seen[key]
    best = max(len(s) for s in seen.values())
    return nch, len(total), best, best == len(total), len(seen)

def contains(w, pat):
    k = len(pat)
    for idx in combinations(range(len(w)), k):
        vals = [w[i] for i in idx]
        if all((vals[i] < vals[j]) == (pat[i] < pat[j]) for i in range(k) for j in range(i + 1, k)):
            return idx
    return None
