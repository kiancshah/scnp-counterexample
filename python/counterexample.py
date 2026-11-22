"""Everything about the counterexample u = 214365 <= w = 541632 (Section 6 of the paper), recomputed from the definitions.

Run:  python3 counterexample.py      (about 10 seconds)
"""
from core import *
from tokens import run_chain, esc, g_of, inside

u, w = P("214365"), P("541632")
n = 6
SIX = [(1, 2), (5, 6), (1, 4), (3, 6), (1, 5), (2, 6)]
BOUND = {(1, 2): 0, (5, 6): 0, (1, 4): 2, (3, 6): 2, (1, 5): 4, (2, 6): 4}
X = [(1, 3), (1, 6), (2, 3), (4, 6), (2, 6), (4, 5), (2, 5)]
Y = [(7 - b, 7 - a) for a, b in X]


def walk(steps):
    v = list(u)
    seq = [S(v)]
    for a, b in steps:
        assert (a, b) in [(x, y) for x, y, _ in covers_up(tuple(v))], "not a cover"
        v[a - 1], v[b - 1] = v[b - 1], v[a - 1]
        seq.append(S(v))
    assert tuple(v) == w
    return seq


if __name__ == "__main__":
    print("u =", S(u), " w =", S(w), " u <= w:", bruhat_leq(u, w))
    print("Inv(u) =", sorted(inv(u)))
    print("Inv(w) =", sorted(inv(w)))
    print("w contains 4231:", contains(w, (4, 2, 3, 1)) is not None, "| w contains 3412:", contains(w, (3, 4, 1, 2)) is not None)
    Iu, Iw = inv(u), inv(w)
    delta = {W: sum(1 for t in Iw if inside(t, W)) - sum(1 for t in Iu if inside(t, W)) for W in windows(n)}
    print("Delta on the six windows:", {W: delta[W] for W in SIX})
    for name, ch in (("X", X), ("Y", Y)):
        g = g_of(n, ch)
        print(f"chain {name}:", " < ".join(walk(ch)))
        print(f"   steps {ch}")
        print(f"   g on the six windows: {[g[W] for W in SIX]}   |supp m_{name}| = {len(chain_support(n, ch))}")
        print("   contains x2^2 x4^2 x5^3:", (0, 2, 0, 2, 3) in chain_support(n, ch),
              "| contains x1^3 x2^2 x4^2:", (3, 2, 0, 2, 0) in chain_support(n, ch))
    chains = list(all_chains(u, w))
    ok = [c for c in chains if all(g_of(n, c)[W] <= BOUND[W] for W in SIX)]
    print("number of chains:", len(chains))
    print("chains meeting all six bounds (Lemma 6.2 says none):", len(ok))
    nch, total, best, good, multisets = brute(u, w)
    print(f"|supp D| = {total}, largest single-chain support = {best}, distinct step multisets = {multisets}, SCNP: {good}")
    mins = minimal_g(u, w)
    print("minimal vectors g:", len(mins))
    Wn = windows(n)
    for m in mins:
        cs = [c for c in chains if tuple(g_of(n, c)[W] for W in Wn) == m]
        ends, moved = set(), set()
        for c in cs:
            _, toks = run_chain(u, c)
            for t in toks:
                if t['kind'] == 'initial' and t['start'] == (3, 4):
                    ends.add(t['hist'][-1])
                elif t['hist'][-1] != t['start']:
                    moved.add((t['kind'], t['start'], t['hist'][-1]))
        print(f"   {len(cs)} chains, first step {sorted({c[0] for c in cs})}, lost token ends at {sorted(ends)}, other moved tokens {sorted(moved)}")
    uu, ww = P("236145"), P("563412")
    print("mirror pair 236145 <= 563412:", brute(uu, ww)[:4], "| 236145 contains 1324:", contains(uu, (1, 3, 2, 4)) is not None)
