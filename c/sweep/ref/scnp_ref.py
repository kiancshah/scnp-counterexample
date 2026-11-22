#!/usr/bin/env python3
"""Independent reference for the SCNP check (slow; small n only).

Enumerates every saturated chain u = v0 < v1 < ... < vl = w in Bruhat order
(covers v < v*t_ab, i.e. swapping POSITIONS a < b), expands each chain product
m_C = prod (y_a + ... + y_{b-1}) explicitly, and compares the support of
D_{u,w} = (1/l!) sum_C m_C with the largest single-chain support.

  scnp_ref.py pair 214365 541632     -> chains, |supp D|, best single chain
  scnp_ref.py list 5                 -> every failing pair of S_5 ("u w len")
"""
import itertools
import sys


def length(p):
    return sum(p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p)))


def leq(u, w):
    """Tableau criterion: sorted prefixes of u are entrywise <= those of w."""
    return all(all(x <= y for x, y in zip(sorted(u[:k]), sorted(w[:k])))
               for k in range(1, len(u)))


def covers(v):
    n = len(v)
    for a in range(n):
        for b in range(a + 1, n):
            if v[a] < v[b] and not any(v[a] < v[c] < v[b] for c in range(a + 1, b)):
                x = list(v)
                x[a], x[b] = x[b], x[a]
                yield tuple(x), a, b


def chains(u, w):
    out = []

    def rec(v, steps):
        if v == w:
            out.append(steps)
            return
        for x, a, b in covers(v):
            if leq(x, w):
                rec(x, steps + [(a, b)])

    rec(u, [])
    return out


def support(steps, nvar):
    mons = {(0,) * nvar}
    for a, b in steps:
        mons = {m[:k] + (m[k] + 1,) + m[k + 1:] for m in mons for k in range(a, b)}
    return mons


def analyse(u, w):
    nvar = len(u) - 1
    cs = chains(u, w)
    sups = [support(s, nvar) for s in cs]
    D = set().union(*sups)
    return len(cs), len(D), max(len(s) for s in sups)


def parse(s):
    return tuple(10 if c in "0Aa" else int(c) for c in s)


def fmt(p):
    return "".join("A" if x == 10 else str(x) for x in p)


def main():
    if len(sys.argv) == 4 and sys.argv[1] == "pair":
        u, w = parse(sys.argv[2]), parse(sys.argv[3])
        assert leq(u, w), "u is not <= w"
        c, s, b = analyse(u, w)
        print(f"chains={c} supp={s} best_chain={b} scnp={'yes' if b == s else 'NO'}")
    elif len(sys.argv) == 3 and sys.argv[1] == "list":
        n = int(sys.argv[2])
        P = list(itertools.permutations(range(1, n + 1)))
        for u in P:
            for w in P:
                if leq(u, w):
                    c, s, b = analyse(u, w)
                    if b < s:
                        print(fmt(u), fmt(w), length(w) - length(u))
    else:
        print(__doc__)
        sys.exit(1)


if __name__ == "__main__":
    main()
