#!/usr/bin/env python3
"""Independent reference for `scnp greedy` (slow; n <= 5).

Lists every chain of every interval [u, w], decides greediness straight from
An-Tung-Zhang, Def. 3.5, and dominance by expanding supports.

  greedy_ref.py 2143 4312     -> the chains of one interval, marked GREEDY / supports
  greedy_ref.py 5             -> totals for S_5
"""
import itertools
import sys

from scnp_ref import chains, covers, leq, parse, support


def cover_positions(lo, hi):
    for x, a, b in covers(lo):
        if x == hi:
            return a, b
    return None


def walk(u, steps):
    v, seq = list(u), [tuple(u)]
    for a, b in steps:
        v[a], v[b] = v[b], v[a]
        seq.append(tuple(v))
    return seq


def is_greedy(u, seq):
    """No step w_{i-1} < w_i = w_{i-1} t_ab can be replaced by a step into w_i
    from some w' >= u with positions (a, b'), b' > b, or (a', b), a' < a."""
    n = len(u)
    for prev, cur in zip(seq, seq[1:]):
        a, b = cover_positions(prev, cur)
        alternatives = [(a, b2) for b2 in range(b + 1, n)] + [(a2, b) for a2 in range(a)]
        for p, q in alternatives:
            x = list(cur)
            x[p], x[q] = x[q], x[p]
            x = tuple(x)
            if cover_positions(x, cur) and leq(u, x):
                return False
    return True


def analyse(u, w):
    cs = chains(u, w)
    sups = [support(s, len(u) - 1) for s in cs]
    D = set().union(*sups)
    greedy = [is_greedy(u, walk(u, s)) for s in cs]
    dom = [len(S) == len(D) for S in sups]
    return cs, sups, D, greedy, dom


def main():
    if len(sys.argv) == 3:
        u, w = parse(sys.argv[1]), parse(sys.argv[2])
        cs, sups, D, greedy, dom = analyse(u, w)
        for s, S, g in zip(cs, sups, greedy):
            seq = " < ".join("".join(map(str, p)) for p in walk(u, s))
            labels = " ".join(f"x{a + 1}..x{b}" for a, b in s)
            print(f"{seq}   labels {labels}   supp {len(S)}/{len(D)}{'   GREEDY' if g else ''}")
        return
    n = int(sys.argv[1])
    P = list(itertools.permutations(range(1, n + 1)))
    tot = {"pairs": 0, "scnp": 0, "some greedy dominant": 0, "every greedy dominant": 0,
           "no greedy chain": 0}
    for u in P:
        for w in P:
            if leq(u, w):
                _, _, _, greedy, dom = analyse(u, w)
                tot["pairs"] += 1
                tot["scnp"] += any(dom)
                tot["some greedy dominant"] += any(g and d for g, d in zip(greedy, dom))
                tot["every greedy dominant"] += all(d for g, d in zip(greedy, dom) if g)
                tot["no greedy chain"] += not any(greedy)
    print(f"n={n}: " + ", ".join(f"{k} {v}" for k, v in tot.items()))


if __name__ == "__main__":
    main()
