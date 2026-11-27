"""Section 4.3 of the note: the rotations of the six cubes of Table 3, the swaps
they perform (Table 4), and what changes when the first swap is moved to the end."""

CUBES = [  # bottom u, top w, the atoms in cyclic order around the loop
    ("1324", "4231", "12 24 34 13"),
    ("1432765", "5276143", "37 36 35 45 14 12 27"),
    ("3216547", "5472163", "15 25 35 34 47 67 16"),
    ("3217564A98", "54821A9376", "4A 6A 69 68 78 17 15 25 35 34"),
    ("3217654A98", "5482A19376", "17 16 26 36 34 4A 5A 59 58 78"),
    ("3217654A98", "64821A9375", "4A 6A 69 68 78 17 15 25 35 34"),
]


def perm(word):
    return tuple(int(c, 16) for c in word)


def show(pair):
    return "t" + "".join("%X" % x for x in pair)


def swap(v, a, b):
    t = list(v)
    t[a - 1], t[b - 1] = t[b - 1], t[a - 1]
    return tuple(t)


def leq(u, w):
    """Bruhat order: compare, for every i and j, how many of the first i entries are >= j."""
    n = len(u)
    return all(sum(x >= j for x in u[:i]) <= sum(x >= j for x in w[:i])
               for i in range(1, n) for j in range(1, n + 1))


def covers(v):
    n = len(v)
    for a in range(n):
        for b in range(a + 1, n):
            between = any(v[a] < v[c] < v[b] for c in range(a + 1, b))
            if v[a] < v[b] and not between:
                yield swap(v, a + 1, b + 1), (a + 1, b + 1)


def rotation(u, w, atoms):
    """Add the atoms in the given order; return the swaps and the values they exchange."""
    v, swaps, values = u, [], []
    for a, b in atoms:
        atom = swap(u, a, b)
        (x, step), = [(x, s) for x, s in covers(v) if leq(x, w) and leq(atom, x)]
        swaps.append(step)
        values.append(tuple(sorted((v[step[0] - 1], v[step[1] - 1]))))
        v = x
    assert v == w
    return swaps, values


for bottom, top, order in CUBES:
    u, w = perm(bottom), perm(top)
    atoms = [perm(a) for a in order.split()]
    k = len(atoms)
    rows = [rotation(u, w, atoms[r:] + atoms[:r]) for r in range(k)]
    print("\n%s <= %s" % (bottom, top))
    for r in range(k):
        swaps, values = rows[r]
        nxt, nxt_values = rows[(r + 1) % k]
        moved = swaps[1:] + swaps[:1]
        changed = [i for i in range(k) if moved[i] != nxt[i]]
        if changed == [k - 1]:
            kind = "slide %s -> %s, same values: %s" % (
                show(moved[-1]), show(nxt[-1]), values[0] == nxt_values[-1])
        elif k - 1 not in changed and sorted(moved) == sorted(nxt):
            kind = "same swaps in another order"
        elif k - 1 not in changed and len(changed) == 2:
            kind = "uncrossing %s -> %s" % (
                ",".join(show(moved[i]) for i in changed),
                ",".join(show(nxt[i]) for i in changed))
        else:
            kind = "mixed"
        first = swaps[0]
        print("  start (%s)  %s  %s| %s" % (
            order.split()[r], " ".join(show(s) for s in swaps),
            "neighbour swap first " if first[1] - first[0] == 1 else "", kind))
