from sage.combinat.permutation import Permutation
from sage.rings.polynomial.polynomial_ring_constructor import (
    PolynomialRing)
from sage.rings.rational_field import QQ

R = PolynomialRing(QQ, 'x', 5)
x = R.gens()
u = Permutation([2, 1, 4, 3, 6, 5])
w = Permutation([5, 4, 1, 6, 3, 2])


def label(p, q):
    d = [i for i in range(len(p)) if p[i] != q[i]]
    return sum(x[i] for i in range(d[0], d[1]))


def weights(p):
    if p == w:
        yield R.one()
        return
    for q in p.bruhat_succ():
        if q.bruhat_lequal(w):
            for m in weights(q):
                yield label(p, q) * m


supports = [frozenset(m.exponents()) for m in weights(u)]
total = frozenset().union(*supports)
print("w avoids 4231:", w.avoids([4, 2, 3, 1]))
print("chains:", len(supports))
print("monomials of D:", len(total))
print("most monomials from one chain:", max(len(s) for s in supports))
print("some chain has all monomials:", any(s == total for s in supports))
