#!/usr/bin/env python3
"""
scnp_criterion.py: the fast SCNP test, written out in plain Python.

This is the method the C checker uses, without any optimisation, so that it
can be read next to scnp_reference.py (which expands polynomials and is
obviously correct but slow). Both are compared in ../tests/.

THE TEST (writeup/scnp_explained.tex, Section 8):

  A "block" is a set of consecutive variable indices R = {i, ..., j} with
  1 <= i <= j <= n-1. There are n(n-1)/2 blocks. For a chain C let

      z_C(R) = number of labels of C whose block lies inside R.

  (1) Comparing two chains [writeup Cor. 8.4, from Hall's theorem and
      Postnikov, "Permutohedra, associahedra, and beyond", Prop. 6.3]:
          supp m_C' is contained in supp m_C  <=>  z_C <= z_C' entrywise.
  (2) Hence D_{u,w} has SCNP iff some chain C attains the entrywise
      minimum H = min over all chains of z_C, in every entry at once.
  (3) The minimum and the "attained" flag can be computed one step at a
      time. Here we go DOWN from a fixed top w: for every v <= w,
          H(v)  = entrywise min of z_C over chains C from v up to w,
          ok(v) = True iff one chain from v to w attains H(v).
      H(w) = 0 and ok(w) = True. For v < w, over the steps v -> x with x <= w
      (label block B, and delta_B(R) = 1 if B lies inside R, else 0):
          H(v)  = min over steps of  H(x) + delta_B,
          ok(v) = some step has ok(x) and H(x) + delta_B == H(v).
      Then (v, w) is bad exactly when ok(v) is False. [writeup Thm. 8.7
      states this recursion going up from u; going down from w is the same
      argument with the chains read backwards.]

  So ONE pass down from w decides SCNP for every pair (u, w), u <= w.

Usage:
  python3 scnp_criterion.py top 541632      # every bad u below w
  python3 scnp_criterion.py table 6         # all pairs of S_6 (a few seconds)
"""

import itertools
import sys
from collections import Counter

from scnp_reference import (atz_prediction, format_perm, length, parse_perm,
                            steps_down)


def all_blocks(n):
    """The blocks {i..j}, 1 <= i <= j <= n-1, as pairs (i, j), in a fixed order."""
    return [(i, j) for i in range(1, n) for j in range(i, n)]


def delta(label_block, blocks):
    """The vector delta_B: entry 1 at every block R containing B."""
    lo, hi = label_block
    return tuple(1 if (i <= lo and hi <= j) else 0 for (i, j) in blocks)


def add(h, d):
    return tuple(x + y for x, y in zip(h, d))


def entrywise_min(h, g):
    return tuple(min(x, y) for x, y in zip(h, g))


def sweep_down(w):
    """Run the recursion (3) from the top w.

    Returns two dicts over all v <= w: H[v] (a tuple, one entry per block)
    and ok[v] (bool). Permutations are processed one length-level at a time,
    from l(w) down to 0, so every x is finished before its lower neighbours
    use it.
    """
    n = len(w)
    blocks = all_blocks(n)
    H = {w: (0,) * len(blocks)}
    ok = {w: True}
    level = [w]
    while level:
        # candidates[v] = list of (H(x) + delta_B, ok(x)) over steps v -> x
        candidates = {}
        for x in level:
            for v, a, b in steps_down(x):
                cand = add(H[x], delta((a, b - 1), blocks))
                candidates.setdefault(v, []).append((cand, ok[x]))
        next_level = []
        for v, cands in candidates.items():
            h = cands[0][0]
            for c, _ in cands[1:]:
                h = entrywise_min(h, c)
            H[v] = h
            ok[v] = any(okx and c == h for c, okx in cands)
            next_level.append(v)
        level = next_level
    return H, ok


def bad_bottoms(w):
    """All u <= w for which (u, w) is bad."""
    _, ok = sweep_down(w)
    return sorted(v for v in ok if not ok[v])


def cmd_top(w_text):
    w = parse_perm(w_text)
    bad = bad_bottoms(w)
    print("w = %s: %d bad bottoms u" % (w_text, len(bad)))
    for u in bad:
        print("  u = %s  length %d" % (format_perm(u), length(w) - length(u)))


def cmd_table(n):
    pairs_by_len = Counter()
    bad_by_len = Counter()
    wrong = []
    for w in itertools.permutations(range(1, n + 1)):
        _, ok = sweep_down(w)
        for v, good in ok.items():
            d = length(w) - length(v)
            pairs_by_len[d] += 1
            if not good:
                bad_by_len[d] += 1
        is_bad_top = not all(ok.values())
        if is_bad_top != atz_prediction(w):
            wrong.append(format_perm(w))
    print("S_%d: %d pairs, %d bad pairs" % (n, sum(pairs_by_len.values()),
                                            sum(bad_by_len.values())))
    for d in sorted(pairs_by_len):
        print("  length %2d: %8d pairs %6d bad" % (d, pairs_by_len[d],
                                                    bad_by_len[d]))
    print("tops where Conj. 5.2(ii) is wrong: %d %s" % (len(wrong), wrong))


def main(argv):
    if len(argv) == 3 and argv[1] == "top":
        cmd_top(argv[2])
    elif len(argv) == 3 and argv[1] == "table":
        cmd_table(int(argv[2]))
    else:
        print(__doc__)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
