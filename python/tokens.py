"""Token bookkeeping along a chain, straight from the one-step lemma."""
from core import *

def run_chain(u, steps):
    """Follow tokens along the chain given by its steps.
    Returns (w, tokens) where each token is a dict(kind='initial'|'born', start=(i,j), path=[(time,place),...])."""
    v = tuple(u); n = len(v)
    tokens = [dict(kind='initial', start=t, born=0, hist=[t]) for t in sorted(inv(v))]
    at = {t['start']: t for t in tokens}
    for step, (a, b) in enumerate(steps, 1):
        al, be = v[a - 1], v[b - 1]
        assert al < be and not any(al < v[c - 1] < be for c in range(a + 1, b)), "not a cover"
        moves = {}
        for c in range(1, a):
            if al < v[c - 1] < be: moves[(c, a)] = (c, b)
        for c in range(b + 1, n + 1):
            if al < v[c - 1] < be: moves[(b, c)] = (a, c)
        new_at = {}
        for place, tok in at.items():
            np_ = moves.get(place, place)
            new_at[np_] = tok
        x = list(v); x[a - 1], x[b - 1] = x[b - 1], x[a - 1]; v = tuple(x)
        born = dict(kind='born', start=(a, b), born=step, hist=[None] * step + [(a, b)])
        new_at[(a, b)] = born; tokens.append(born)
        for place, tok in new_at.items():
            if tok is not born: tok['hist'].append(place)
        at = new_at
        assert set(at) == set(inv(v)), "one-step lemma violated"
    return v, tokens

def inside(place, W):
    return place is not None and W[0] <= place[0] and place[1] <= W[1]

def esc(tokens, W):
    return sum(1 for t in tokens if any(inside(p, W) for p in t['hist']) and not inside(t['hist'][-1], W))

def g_of(n, steps):
    return {W: sum(1 for a, b in steps if W[0] <= a and b <= W[1]) for W in windows(n)}

def check_escape_formula(u, w):
    n = len(u); cnt = 0
    Iu, Iw = inv(u), inv(w)
    for ch in all_chains(u, w):
        v, toks = run_chain(u, ch)
        assert v == w
        g = g_of(n, ch)
        for W in windows(n):
            d = sum(1 for t in Iw if inside(t, W)) - sum(1 for t in Iu if inside(t, W))
            assert g[W] == d + esc(toks, W), (S(u), S(w), ch, W)
        cnt += 1
    return cnt

if __name__ == "__main__":
    import itertools, random
    tot = 0; pairs = 0
    for n in (3, 4, 5):
        perms = list(permutations(range(1, n + 1)))
        for u in perms:
            for w in perms:
                if bruhat_leq(u, w) and length(w) - length(u) <= 7:
                    tot += check_escape_formula(u, w); pairs += 1
    print("escape formula verified on", tot, "chains of", pairs, "pairs (n<=5, length<=7)")
    print("escape formula on 214365<=541632:", check_escape_formula(P("214365"), P("541632")), "chains")
    # example chains for 1324 <= 4231
    u, w = P("1324"), P("4231")
    for m in minimal_g(u, w):
        for ch in all_chains(u, w):
            g = g_of(4, ch)
            if tuple(g[W] for W in windows(4)) == m:
                v, toks = run_chain(u, ch)
                moved = [(t['kind'], t['start'], t['hist'][-1]) for t in toks if t['hist'][-1] != t['start']]
                seq = [S(u)]; x = list(u)
                for a, b in ch:
                    x[a-1], x[b-1] = x[b-1], x[a-1]; seq.append(S(x))
                print("  ", " < ".join(seq), "steps", ch, "moved:", moved, "escapes:", [W for W in windows(4) if esc(toks, W)])
    print("2143 <= 4312 good:", scnp(P("2143"), P("4312")), " inv(2143)", sorted(inv(P("2143"))), "inv(4312)", sorted(inv(P("4312"))))
