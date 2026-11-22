"""Re-run every check that fits in a few minutes of plain Python.

Run:  python3 verify_small.py          (n <= 6, about 3 minutes)
      python3 verify_small.py --more   (also n = 7 for the construction; about 1 more minute)

What is checked, each time against the definitions in core.py:
  1. the window criterion for SCNP agrees with expanding every chain (all pairs of S_4 and S_5);
  2. the table of bad pairs and bad permutations for n = 4, 5, 6, and the pattern statement;
  3. the explicit bottom u for every w containing 4231 (Theorem 5.4);
  4. the weak order theorem: Inv(u) inside Inv(w) forces SCNP with least vector Delta (Theorem 1.5);
  5. the escape formula g = Delta + esc on every chain (Proposition 3.3).
"""
import sys, time
from core import *
from tokens import check_escape_formula, inside
from construction import check as check_construction

PATTERNS = [P("4231"), P("541632"), P("5276143"), P("5472163")]


def all_bottoms(w):
    """For a fixed top w: the Pareto-minimal vectors g of the chains v -> w, for every v <= w."""
    n = len(w)
    ident = tuple(range(1, n + 1))
    levels = interval(ident, w)
    Sv = {w: [tuple([0] * len(windows(n)))]}
    for lev in reversed(levels[:-1]):
        for v in lev:
            cur = []
            for a, b, x in covers_up(v):
                if x in Sv:
                    e = step_vector(n, a, b)
                    for g in Sv[x]:
                        cur = pareto_insert(cur, add(e, g))
            Sv[v] = cur
    return Sv


def delta_vector(u, w):
    Iu, Iw = inv(u), inv(w)
    return tuple(sum(1 for t in Iw if inside(t, W)) - sum(1 for t in Iu if inside(t, W)) for W in windows(len(u)))


if __name__ == "__main__":
    more = "--more" in sys.argv
    t0 = time.time()

    print("1. criterion against brute force")
    cnt = bad = 0
    for n in (4, 5):
        perms = list(permutations(range(1, n + 1)))
        for u in perms:
            for w in perms:
                if bruhat_leq(u, w):
                    cnt += 1
                    a = scnp(u, w)
                    assert a == brute(u, w)[3]
                    bad += (not a)
    print(f"   agree on all {cnt} pairs of S_4 and S_5 ({bad} bad pairs)")

    print("2. table, and 4. weak order theorem")
    print("   n | pairs | bad pairs | bad w | w containing 4231 | bad w = w containing one of the four patterns")
    for n in (4, 5, 6):
        pairs = badpairs = weak = 0
        badw = set()
        for w in permutations(range(1, n + 1)):
            Sv = all_bottoms(w)
            for v, vecs in Sv.items():
                pairs += 1
                if len(vecs) > 1:
                    badpairs += 1
                    badw.add(w)
                if inv(v) <= inv(w):
                    weak += 1
                    assert vecs == [delta_vector(v, w)], (v, w)
        c4231 = sum(1 for w in permutations(range(1, n + 1)) if contains(w, PATTERNS[0]) is not None)
        cany = {w for w in permutations(range(1, n + 1)) if any(len(p) <= n and contains(w, p) is not None for p in PATTERNS)}
        print(f"   {n} | {pairs} | {badpairs} | {len(badw)} | {c4231} | {badw == cany}    (weak order pairs, all good with g = Delta: {weak})")

    print("3. the explicit bottom for tops containing 4231")
    for n in range(4, 8 if more else 7):
        c = 0
        for w in permutations(range(1, n + 1)):
            if contains(w, PATTERNS[0]) is None:
                continue
            u = check_construction(w)
            assert not scnp(u, w)
            c += 1
        print(f"   n = {n}: {c} permutations, the pair (u, w) is bad every time")

    print("5. escape formula")
    tot = 0
    for n in (3, 4, 5):
        perms = list(permutations(range(1, n + 1)))
        for u in perms:
            for w in perms:
                if bruhat_leq(u, w) and length(w) - length(u) <= 6:
                    tot += check_escape_formula(u, w)
    print(f"   holds on {tot} chains (n <= 5, length <= 6)")
    print(f"all checks passed in {time.time() - t0:.0f} s")
