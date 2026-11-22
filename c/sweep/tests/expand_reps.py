#!/usr/bin/env python3
"""Expand orbit representatives ("u w len" lines, as written by
`scnp sweep --list reps`) to the full symmetry orbits:
(u,w), (w0 w, w0 u), (w0 u w0, w0 w w0), (w w0, u w0).  Prints unique lines."""
import sys

def parse(s):
    return [10 if c in "0Aa" else int(c) for c in s]

def fmt(p):
    return "".join("A" if x == 10 else str(x) for x in p)

for line in sys.stdin:
    a, b, d = line.split()[:3]
    u, w = parse(a), parse(b)
    n = len(u)
    lw0 = lambda p: [n + 1 - x for x in p]          # w0 * p
    rw0 = lambda p: p[::-1]                         # p * w0
    seen = set()
    for x, y in ((u, w), (lw0(w), lw0(u)), (lw0(rw0(u)), lw0(rw0(w))), (rw0(w), rw0(u))):
        key = (fmt(x), fmt(y))
        if key not in seen:
            seen.add(key)
            print(key[0], key[1], d)
