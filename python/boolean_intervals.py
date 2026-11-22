"""Observation 9.2 of the paper, checked exhaustively for S_n.

  (a) every bad pair (u, w) contains a bad pair (u', w'), u <= u' <= w' <= w, whose interval is Boolean;
  (b) if [u, w] is Boolean and (u, w) is bad, then the atom graph of (u, w) has a cycle.

Run:  cd python
      python3 boolean_intervals.py 4 5 6 7     (a few seconds)
      python3 boolean_intervals.py 8           (15 minutes on two cores, about 1.5 GB of memory)

Input: the list of all bad pairs of S_n, one pair "u w length" per line, as written by the C sweep.
  n = 4, 5, 6, 7:  data/sweep/scnp_fails_S{n}.txt.gz, stored in this repository.
  n = 8:           not stored (29 131 995 lines).  Make it first, in a few seconds, with
                     cd c/sweep && make && mkdir -p runs && ./scnp sweep -n 8 -o runs/s8 --list all --gzip

Boolean is tested strictly: [u, w] has length l, exactly l atoms and 2^l elements, and the map
v -> {atoms below v} is injective.  Counting elements alone is not enough.

The Bruhat order of S_n is kept as bitsets: Python integers in which bit j of down[i] is set
when perm j <= perm i.
"""
import glob
import gzip
import os
import sys
from collections import Counter, defaultdict
from itertools import permutations
from multiprocessing import Pool

from core import P, S, length


def bits(m):
    """indices of the set bits of m"""
    while m:
        low = m & -m
        yield low.bit_length() - 1
        m ^= low


class Bruhat:
    def __init__(self, n):
        self.n = n
        self.perms = list(permutations(range(1, n + 1)))
        self.index = {p: i for i, p in enumerate(self.perms)}
        self.length = [length(p) for p in self.perms]
        self.cov = [self.upper_covers(p) for p in self.perms]
        order = sorted(range(len(self.perms)), key=lambda i: self.length[i])
        self.down = [1 << i for i in range(len(self.perms))]
        for i in order:                  # the elements below i are complete before i is used
            for j, _ in self.cov[i]:
                self.down[j] |= self.down[i]
        self.up = [1 << i for i in range(len(self.perms))]
        for i in reversed(order):        # the elements above the covers of i are complete first
            for j, _ in self.cov[i]:
                self.up[i] |= self.up[j]

    def upper_covers(self, p):
        """(index of p t_ab, (a, b)) for every cover; positions a < b are 1-based"""
        out = []
        for a in range(self.n):
            for b in range(a + 1, self.n):
                if p[a] < p[b] and not any(p[a] < p[c] < p[b] for c in range(a + 1, b)):
                    x = list(p)
                    x[a], x[b] = x[b], x[a]
                    out.append((self.index[tuple(x)], (a + 1, b + 1)))
        return out

    def leq(self, v, w):
        return self.down[w] >> v & 1

    def atoms(self, u, w):
        """the elements covering u inside [u, w], as (index, swapped positions)"""
        return [(x, s) for x, s in self.cov[u] if self.leq(x, w)]

    def boolean(self, u, w):
        """the atoms of [u, w] if the interval is Boolean, else None"""
        l = self.length[w] - self.length[u]
        atoms = self.atoms(u, w)
        if l < 1 or len(atoms) != l:
            return None
        inside = self.up[u] & self.down[w]
        if inside.bit_count() != 2 ** l:
            return None
        # split [u, w] by "is atom x below?" for every atom: injective iff 2^l nonempty blocks
        blocks = [inside]
        for x, _ in atoms:
            above = self.up[x] & inside
            blocks = [c for b in blocks for c in (b & above, b & ~above) if c]
        return atoms if len(blocks) == 2 ** l else None


def has_cycle(edges):
    """does the graph with these edges contain a cycle?  (union-find)"""
    parent = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            x = parent[x]
        return x

    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra == rb:
            return True
        parent[ra] = rb
    return False


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def fails_files(n):
    stored = os.path.join(ROOT, "data", "sweep", f"scnp_fails_S{n}.txt.gz")
    if os.path.exists(stored):
        return [stored]
    return sorted(glob.glob(os.path.join(ROOT, "c", "sweep", "runs", f"s{n}", "chunks", "*.txt.gz")))


def batches(n, size=20000):
    batch = []
    for f in fails_files(n):
        for line in gzip.open(f, "rt"):
            batch.append(line)
            if len(batch) == size:
                yield batch
                batch = []
    if batch:
        yield batch


B = None


def study(batch):
    """the Boolean pairs among a batch of bad pairs, each with 'does the atom graph have a cycle'"""
    found = []
    for line in batch:
        us, ws, _ = line.split()
        u, w = B.index[P(us)], B.index[P(ws)]
        atoms = B.boolean(u, w)
        if atoms is not None:
            found.append((u, w, len(atoms), has_cycle([s for _, s in atoms])))
    return found, len(batch)


def run(n, show=10):
    global B
    if not fails_files(n):
        print(f"S_{n}: no list of bad pairs found (see the top of this file)")
        return
    B = Bruhat(n)
    N = len(B.perms)

    # (b): the Boolean pairs among the bad pairs, and their atom graphs
    bad_boolean, n_bad = [], 0
    with Pool() as pool:
        for found, m in pool.imap(study, batches(n)):
            bad_boolean += found
            n_bad += m
    forest = [(u, w) for u, w, _, cyc in bad_boolean if not cyc]
    by_rank = Counter(l for _, _, l, _ in bad_boolean)
    print(f"S_{n}: {n_bad} bad pairs, {len(bad_boolean)} of them Boolean")
    print(f"   (b) bad Boolean pairs whose atom graph has no cycle: {len(forest)}")
    for u, w in forest[:show]:
        print(f"       {S(B.perms[u])} <= {S(B.perms[w])}")
    print("   bad Boolean pairs by rank: " + ", ".join(f"{l}: {by_rank[l]}" for l in sorted(by_rank)))

    # (a): contains[u] = set of w such that [u, w] contains a bad Boolean pair
    #      = up-sets of the tops of the bad Boolean pairs with bottom u, and contains[x] for x covering u
    own = [0] * N
    for u, w, _, _ in bad_boolean:
        own[u] |= B.up[w]
    contains = [0] * N
    for u in sorted(range(N), key=lambda i: -B.length[i]):
        c = own[u]
        for x, _ in B.cov[u]:
            c |= contains[x]
        contains[u] = c
    rows = defaultdict(list)
    for f in fails_files(n):
        for line in gzip.open(f, "rt"):
            us, ws, _ = line.split()
            rows[B.index[P(us)]].append(B.index[P(ws)])
    pairs = n_contain = bad_without = good_with = 0
    examples = []
    for u in range(N):
        bad_u = 0
        for w in rows.get(u, ()):
            bad_u |= 1 << w
        pairs += B.up[u].bit_count()
        n_contain += contains[u].bit_count()
        missing = bad_u & ~contains[u]
        bad_without += missing.bit_count()
        good_with += (contains[u] & ~bad_u).bit_count()
        for w in bits(missing):
            if len(examples) < show:
                examples.append((u, w))
    print(f"   {pairs} pairs u <= w, {n_contain} of them contain a bad Boolean pair")
    print(f"   (a) bad pairs containing no bad Boolean pair: {bad_without}")
    for u, w in examples:
        print(f"       {S(B.perms[u])} <= {S(B.perms[w])}")
    print(f"   good pairs containing a bad Boolean pair (the converse of (a) fails): {good_with}")


if __name__ == "__main__":
    for n in map(int, sys.argv[1:]):
        run(n)
