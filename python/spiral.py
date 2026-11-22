"""Boolean bad intervals: atoms, the atom graph, and competing chains as rotations of one cyclic order.

Run:  python3 spiral.py            (the two intervals of size 7; a few seconds)
      python3 spiral.py --all      (also the three Boolean intervals of size 10; several minutes)

Words used here.
  atom            an element covering u inside [u, w]; it is u with two positions a < b swapped, written (a, b)
  atom graph      vertices = positions 1..n, one edge {a, b} for every atom (a, b)
  Boolean         [u, w] is a Boolean lattice; then a chain is the same as an order in which the atoms are "filled"
  competing       a chain whose vector g is minimal among all chains from u to w
"""
from core import *
import sys


def structure(u, w):
    n = len(u)
    levels = interval(u, w)
    nodes = [v for lev in levels for v in lev]
    atoms = sorted(levels[1])
    atom_swap = {}
    for a, b, x in covers_up(u):
        if x in levels[1]:
            atom_swap[x] = (a, b)
    below = {u: frozenset()}
    for r in range(1, len(levels)):
        for x in levels[r]:
            if r == 1:
                below[x] = frozenset([x])
                continue
            s = set()
            for v in levels[r - 1]:
                if any(y == x for _, _, y in covers_up(v)):
                    s |= below[v]
            below[x] = frozenset(s)
    L = len(levels) - 1
    boolean = (len(atoms) == L and len(nodes) == 2 ** L and len(set(below.values())) == len(nodes))
    return dict(n=n, levels=levels, nodes=nodes, atoms=atoms, atom_swap=atom_swap, below=below, L=L, boolean=boolean)


def chains_with_g(u, w, target):
    """All chains from u to w whose vector g equals target. Partial chains that already exceed target are dropped."""
    n = len(u)
    out = []

    def rec(v, g, steps):
        if v == w:
            if g == target:
                out.append(list(steps))
            return
        for a, b, x in covers_up(v):
            if not bruhat_leq(x, w):
                continue
            g2 = add(g, step_vector(n, a, b))
            if leq(g2, target):
                steps.append((a, b))
                rec(x, g2, steps)
                steps.pop()

    rec(u, tuple([0] * len(windows(n))), [])
    return out


def atom_orders(u, w, st, chain):
    """For a Boolean interval: the order in which the chain fills the atoms."""
    v = u
    cur = frozenset()
    order = []
    for a, b in chain:
        x = list(v)
        x[a - 1], x[b - 1] = x[b - 1], x[a - 1]
        v = tuple(x)
        new = st['below'][v] - cur
        assert len(new) == 1
        order.append(st['atom_swap'][next(iter(new))])
        cur = st['below'][v]
    return order


def chain_from_atom_order(u, st, order):
    """For a Boolean interval: the steps of the chain that fills the atoms in the given order."""
    inv_swap = {v: k for k, v in st['atom_swap'].items()}
    byset = {st['below'][x]: x for x in st['nodes']}
    cur = frozenset()
    v = u
    steps = []
    for s in order:
        cur = cur | {inv_swap[s]}
        x = byset[cur]
        d = [i + 1 for i in range(len(v)) if v[i] != x[i]]
        assert len(d) == 2
        steps.append(tuple(d))
        v = x
    return steps


def g_vector(n, steps):
    g = tuple([0] * len(windows(n)))
    for a, b in steps:
        g = add(g, step_vector(n, a, b))
    return g


def report(us, ws):
    u, w = P(us), P(ws)
    n = len(u)
    st = structure(u, w)
    mins = minimal_g(u, w)
    print(f"{us} <= {ws}: length {st['L']}, {len(st['nodes'])} elements, Boolean: {st['boolean']}")
    print("  atoms (edges of the atom graph):", sorted(st['atom_swap'].values()))
    per = [chains_with_g(u, w, m) for m in mins]
    print("  competing vectors:", len(mins), "| chains attaining each:", [len(c) for c in per])
    if not st['boolean']:
        for cs in per:
            print("    first steps:", sorted({c[0] for c in cs}))
        return
    orders = [[atom_orders(u, w, st, c) for c in cs] for cs in per]
    base = orders[0][0]
    k = len(base)
    rots = {tuple(base[i:] + base[:i]): i for i in range(k)}
    if not all(tuple(o) in rots for os_ in orders for o in os_):
        print("  the competing chains are NOT rotations of one cyclic order")
        return
    used = sorted(rots[tuple(o)] for os_ in orders for o in os_)
    print("  cyclic order of the atoms:", base)
    print("  rotations that compete start at:", [base[i] for i in used])
    G = {r: g_vector(n, chain_from_atom_order(u, st, base[r:] + base[:r])) for r in range(k)}
    for r in range(k):
        if r in used:
            continue
        nxt = (r + 1) % k
        a = sorted(chain_from_atom_order(u, st, base[r:] + base[:r]))
        b = sorted(chain_from_atom_order(u, st, base[nxt:] + base[:nxt]))
        only_a = [x for x in a if x not in b]
        only_b = [x for x in b if x not in a]
        print(f"  missing rotation, start {base[r]}: its g is >= the g of the next rotation (start {base[nxt]}):",
              leq(G[nxt], G[r]), "| steps", only_a, "replaced by", only_b)


if __name__ == "__main__":
    pairs = [("1324", "4231"), ("214365", "541632"), ("1432765", "5276143"), ("3216547", "5472163")]
    if "--all" in sys.argv:
        pairs += [("3217564A98", "54821A9376"), ("3217654A98", "5482A19376"), ("3217654A98", "64821A9375")]
    for us, ws in pairs:
        report(us, ws)
