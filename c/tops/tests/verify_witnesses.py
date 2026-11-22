#!/usr/bin/env python3
"""Re-decide witness pairs "u w len" (e.g. runs/wit/all11.sample) by brute
force, i.e. by expanding every chain product (reference/scnp_reference.py),
without the H-vector criterion used by the fast checker.

At n = 11 the reference's interval() (up-set of u intersected with the
down-set of w, both over all of S_n) is too slow, so the interval is built
here by walking up from u and keeping x <= w by the tableau criterion
[LIT: Bjorner-Brenti, Combinatorics of Coxeter Groups, Thm. 2.1.5].

  verify_witnesses.py FILE [MAXLEN]
"""
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "reference"))
from scnp_reference import label, parse_perm, product_support, steps_up  # noqa


def tableau_leq(u, w):
    n = len(u)
    for k in range(1, n):
        a, b = sorted(u[:k]), sorted(w[:k])
        if any(x > y for x, y in zip(a, b)):
            return False
    return True


def interval(u, w):
    seen, todo = {u}, [u]
    while todo:
        v = todo.pop()
        for x, _, _ in steps_up(v):
            if x not in seen and tableau_leq(x, w):
                seen.add(x)
                todo.append(x)
    return seen


def scnp_brute(u, w):
    inside = interval(u, w)
    memo = {}

    def from_node(v):
        if v in memo:
            return memo[v]
        if v == w:
            res = Counter({(): 1})
        else:
            res = Counter()
            for x, a, b in steps_up(v):
                if x in inside:
                    for rest, c in from_node(x).items():
                        res[tuple(sorted(rest + (label(a, b),)))] += c
        memo[v] = res
        return res

    sups = [product_support(ms, len(u) - 1) for ms in from_node(u)]
    union = set().union(*sups)
    return any(s == union for s in sups)


def main():
    path = sys.argv[1]
    maxlen = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    checked = agree = 0
    for line in open(path):
        ut, wt, d = line.split()[:3]
        if int(d) > maxlen:
            continue
        checked += 1
        good = scnp_brute(parse_perm(ut), parse_perm(wt))
        if not good:
            agree += 1
        else:
            print("DISAGREE: fast checker says bad, brute force says SCNP:", ut, wt)
    print("%d witness pairs of length <= %d re-decided by brute force; "
          "%d confirmed bad" % (checked, maxlen, agree))


if __name__ == "__main__":
    main()
