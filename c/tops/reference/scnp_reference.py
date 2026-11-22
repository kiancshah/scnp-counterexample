#!/usr/bin/env python3
"""
scnp_reference.py: a slow, simple reference implementation of the SCNP
question for Postnikov-Stanley polynomials.

Every function here follows a definition directly and makes no use of any
clever criterion. It is meant to be read line by line and to be trusted,
not to be fast. The fast checker (../fast/scnp_fast.c) is tested against it.

Conventions (identical to the rest of this repository, see the main README
and writeup/scnp_explained.tex):

  * A permutation w of {1,...,n} is a tuple in one-line notation,
    w = (w(1), ..., w(n)). For n = 10 the letter "A" stands for 10.
  * Positions are 1-based in everything the user sees.
  * A STEP (Bruhat cover) v -> x swaps the entries at POSITIONS a < b of v,
    and it is allowed exactly when
        (i)  v(a) < v(b), and
        (ii) no position c with a < c < b has v(a) < v(c) < v(b).
  * The LABEL of that step is the linear form x_a + x_{a+1} + ... + x_{b-1}.
    We store a label as the pair (a, b-1): its "block" of variable indices.
  * D_{u,w} = (1 / l!) * sum over all saturated chains C from u to w of
    m_C = product of the labels of C. The constant 1/l! never matters here.
  * SCNP holds for (u, w) when some single chain C has
    supp(m_C) = supp(D_{u,w}). A pair where SCNP fails is called BAD.
  * w is a BAD TOP when some u <= w makes (u, w) bad.

The conjecture under study (An-Tung-Zhang, Conjecture 5.2(ii) in arXiv v1 of
"Newton polytopes of dual Schubert polynomials", as recorded in this
repository's writeup):  w is a bad top  <=>  w contains the pattern 4231.

Usage:
  python3 scnp_reference.py pair 214365 541632   # one pair, full details
  python3 scnp_reference.py table 5              # every pair of S_5
  python3 scnp_reference.py selftest             # small known facts
"""

import itertools
import sys
from collections import Counter


# ---------------------------------------------------------------------------
# 1. Permutations
# ---------------------------------------------------------------------------

def parse_perm(text):
    """'214365' -> (2, 1, 4, 3, 6, 5).  'A' (or '0') means 10, 'B' 11, 'C' 12."""
    values = []
    for ch in text:
        if ch in "Aa0":
            values.append(10)
        elif ch in "Bb":
            values.append(11)
        elif ch in "Cc":
            values.append(12)
        else:
            values.append(int(ch))
    if sorted(values) != list(range(1, len(values) + 1)):
        raise ValueError("not a permutation: " + text)
    return tuple(values)


def format_perm(w):
    """(2, 1, 4, 3, 6, 5) -> '214365'.  10 is written 'A', 11 'B', 12 'C'."""
    return "".join("123456789ABC"[x - 1] for x in w)


def length(w):
    """Number of inversions: pairs of positions i < j with w(i) > w(j)."""
    n = len(w)
    count = 0
    for i in range(n):
        for j in range(i + 1, n):
            if w[i] > w[j]:
                count += 1
    return count


def swap_positions(w, a, b):
    """w with the entries at 1-based positions a and b exchanged."""
    x = list(w)
    x[a - 1], x[b - 1] = x[b - 1], x[a - 1]
    return tuple(x)


# ---------------------------------------------------------------------------
# 2. Steps (Bruhat covers) and the Bruhat order
# ---------------------------------------------------------------------------

def is_step(v, a, b):
    """True if swapping positions a < b of v is a step going up."""
    if not v[a - 1] < v[b - 1]:                       # condition (i)
        return False
    for c in range(a + 1, b):                          # condition (ii)
        if v[a - 1] < v[c - 1] < v[b - 1]:
            return False
    return True


def steps_up(v):
    """All steps out of v going up, as triples (x, a, b) with x = v t_ab."""
    n = len(v)
    result = []
    for a in range(1, n + 1):
        for b in range(a + 1, n + 1):
            if is_step(v, a, b):
                result.append((swap_positions(v, a, b), a, b))
    return result


def steps_down(x):
    """All steps INTO x from below, as triples (v, a, b) with v -> x a step.

    We simply try every swap and keep those that are steps going up
    from the result back to x; this reuses the single definition is_step.
    """
    n = len(x)
    result = []
    for a in range(1, n + 1):
        for b in range(a + 1, n + 1):
            v = swap_positions(x, a, b)
            if is_step(v, a, b):
                result.append((v, a, b))
    return result


def up_set(u):
    """Every permutation reachable from u by a sequence of steps (u included).

    This IS the definition of the Bruhat order used here: v >= u iff v is in
    up_set(u). No tableau criterion is used.
    """
    seen = {u}
    todo = [u]
    while todo:
        v = todo.pop()
        for x, _, _ in steps_up(v):
            if x not in seen:
                seen.add(x)
                todo.append(x)
    return seen


def down_set(w):
    """Every permutation from which w is reachable by steps (w included)."""
    seen = {w}
    todo = [w]
    while todo:
        x = todo.pop()
        for v, _, _ in steps_down(x):
            if v not in seen:
                seen.add(v)
                todo.append(v)
    return seen


def interval(u, w):
    """The Bruhat interval [u, w] as a set (empty if u is not <= w)."""
    return up_set(u) & down_set(w)


# ---------------------------------------------------------------------------
# 3. Chains, labels and supports
# ---------------------------------------------------------------------------

def label(a, b):
    """Label of the step swapping positions a < b: block (a, b-1),
    standing for the linear form x_a + ... + x_{b-1}."""
    return (a, b - 1)


def all_chains(u, w):
    """Every saturated chain from u to w, as a list of labels in order.

    Only practical for short intervals; it is here so the slow path can be
    compared with the multiset path below.
    """
    inside = interval(u, w)
    chains = []

    def extend(v, labels_so_far):
        if v == w:
            chains.append(list(labels_so_far))
            return
        for x, a, b in steps_up(v):
            if x in inside:
                labels_so_far.append(label(a, b))
                extend(x, labels_so_far)
                labels_so_far.pop()

    if inside:
        extend(u, [])
    return chains


def label_multisets(u, w):
    """The set of label MULTISETS of the chains from u to w, with the number
    of chains having each multiset.

    Two chains with the same multiset of labels have the same product m_C
    (multiplication is commutative), so the support of D_{u,w} only depends
    on this set. A multiset is stored as a sorted tuple of labels.

    Computed from the top down: for each v in [u, w], the chains from v to w
    are a first step v -> x followed by a chain from x to w.
    """
    inside = interval(u, w)
    memo = {}

    def from_node(v):
        if v in memo:
            return memo[v]
        if v == w:
            result = Counter({(): 1})
        else:
            result = Counter()
            for x, a, b in steps_up(v):
                if x in inside:
                    lab = label(a, b)
                    for rest, count in from_node(x).items():
                        result[tuple(sorted(rest + (lab,)))] += count
        memo[v] = result
        return result

    if not inside:
        return Counter()
    return from_node(u)


def product_support(labels, nvars):
    """Support of the product of the labels, by explicit expansion.

    A monomial is an exponent vector (e_1, ..., e_nvars). Expanding a product
    of sums means picking one variable from each factor; all coefficients
    are positive, so nothing cancels and the support is exactly the set of
    exponent vectors obtained this way.
    """
    monomials = {(0,) * nvars}
    for (lo, hi) in labels:
        new = set()
        for m in monomials:
            for k in range(lo, hi + 1):            # pick x_k from this factor
                e = list(m)
                e[k - 1] += 1
                new.add(tuple(e))
        monomials = new
    return monomials


def analyse_pair(u, w):
    """Decide SCNP for u <= w by brute force.

    Returns a dict with the number of chains, |supp D|, the size of the best
    single-chain support, the number of distinct label multisets, and
    whether SCNP holds.
    """
    nvars = len(u) - 1
    multisets = label_multisets(u, w)
    if not multisets:
        raise ValueError("u is not <= w")
    supports = {ms: product_support(ms, nvars) for ms in multisets}
    supp_D = set()
    for s in supports.values():
        supp_D |= s                       # supp D = union of chain supports
    best = max(len(s) for s in supports.values())
    dominant = [ms for ms, s in supports.items() if s == supp_D]
    return {
        "chains": sum(multisets.values()),
        "multisets": len(multisets),
        "supp_D": len(supp_D),
        "best_chain": best,
        "scnp": len(dominant) > 0,
        "dominant_multisets": dominant,
        "supp_D_set": supp_D,
        "supports": supports,
    }


def is_scnp(u, w):
    return analyse_pair(u, w)["scnp"]


# ---------------------------------------------------------------------------
# 4. Patterns and the conjecture
# ---------------------------------------------------------------------------

def contains_pattern(w, p):
    """True if some len(p) entries of w (left to right) are in the same
    relative order as the entries of p."""
    k = len(p)
    for positions in itertools.combinations(range(len(w)), k):
        entries = [w[i] for i in positions]
        if all((entries[i] < entries[j]) == (p[i] < p[j])
               for i in range(k) for j in range(i + 1, k)):
            return True
    return False


PATTERN_4231 = (4, 2, 3, 1)


def atz_prediction(w):
    """What An-Tung-Zhang Conj. 5.2(ii) predicts: w is a bad top iff w
    contains 4231."""
    return contains_pattern(w, PATTERN_4231)


# ---------------------------------------------------------------------------
# 5. Command line
# ---------------------------------------------------------------------------

def show_monomial(e):
    parts = []
    for k, power in enumerate(e, start=1):
        if power == 1:
            parts.append("x%d" % k)
        elif power > 1:
            parts.append("x%d^%d" % (k, power))
    return "*".join(parts) if parts else "1"


def cmd_pair(u_text, w_text):
    u, w = parse_perm(u_text), parse_perm(w_text)
    r = analyse_pair(u, w)
    print("u = %s, w = %s, l(w) - l(u) = %d, |[u,w]| = %d"
          % (u_text, w_text, length(w) - length(u), len(interval(u, w))))
    print("chains = %d, distinct label multisets = %d"
          % (r["chains"], r["multisets"]))
    print("|supp D| = %d, best single chain = %d, SCNP = %s"
          % (r["supp_D"], r["best_chain"], "yes" if r["scnp"] else "NO"))
    if not r["scnp"]:
        # Show, for one best chain, which monomials it misses.
        best_ms = max(r["supports"], key=lambda ms: len(r["supports"][ms]))
        missing = sorted(r["supp_D_set"] - r["supports"][best_ms])
        print("a best chain has labels",
              " ".join("[%d..%d]" % lab for lab in best_ms))
        print("and misses", ", ".join(show_monomial(e) for e in missing))


def cmd_table(n):
    """Brute force over every pair u <= w of S_n."""
    perms = list(itertools.permutations(range(1, n + 1)))
    pairs_by_len = Counter()
    bad_by_len = Counter()
    bad_tops = set()
    for u in perms:
        for w in up_set(u):
            d = length(w) - length(u)
            pairs_by_len[d] += 1
            if not is_scnp(u, w):
                bad_by_len[d] += 1
                bad_tops.add(w)
    print("S_%d: %d pairs, %d bad pairs" % (n, sum(pairs_by_len.values()),
                                            sum(bad_by_len.values())))
    for d in sorted(pairs_by_len):
        print("  length %2d: %8d pairs %6d bad" % (d, pairs_by_len[d],
                                                    bad_by_len[d]))
    wrong = [w for w in perms if (w in bad_tops) != atz_prediction(w)]
    print("bad tops: %d; tops where Conj. 5.2(ii) is wrong: %d %s"
          % (len(bad_tops), len(wrong), [format_perm(w) for w in wrong]))


def cmd_selftest():
    # Example 2.13 of the writeup: 2! D_{213}^{321} = 2 x1 x2 + x2^2.
    r = analyse_pair(parse_perm("213"), parse_perm("321"))
    assert r["supp_D_set"] == {(1, 1), (0, 2)} and r["scnp"]
    # The smallest bad pair (writeup Section 4.1).
    r = analyse_pair(parse_perm("1324"), parse_perm("4231"))
    assert (r["chains"], r["supp_D"], r["best_chain"], r["scnp"]) == \
        (24, 6, 5, False), r
    # Kian Shah's counterexample: 2688 chains, |supp D| = 81, best 73.
    r = analyse_pair(parse_perm("214365"), parse_perm("541632"))
    assert (r["chains"], r["supp_D"], r["best_chain"], r["scnp"]) == \
        (2688, 81, 73, False), r
    assert not contains_pattern(parse_perm("541632"), PATTERN_4231)
    # The chain list and the multiset count agree on a small interval.
    u, w = parse_perm("1324"), parse_perm("4231")
    assert len(all_chains(u, w)) == sum(label_multisets(u, w).values())
    print("selftest ok")


def main(argv):
    if len(argv) == 4 and argv[1] == "pair":
        cmd_pair(argv[2], argv[3])
    elif len(argv) == 3 and argv[1] == "table":
        cmd_table(int(argv[2]))
    elif len(argv) == 2 and argv[1] == "selftest":
        cmd_selftest()
    else:
        print(__doc__)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
