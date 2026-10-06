"""
Sanity tests for fflib.py.  Run:  python -m pytest code/test_fflib.py   (or python code/test_fflib.py)

The numba kernels are checked against an independent, slow, pure-python
factorisation (trial division by all monic irreducibles of degree <= n/2).
"""
import itertools
import os
import sys
from collections import Counter
from math import factorial

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fflib  # noqa: E402


def _polymod(a, m, p):
    a = a[:]
    while len(a) >= len(m):
        c = a[-1] % p
        if c:
            s = len(a) - len(m)
            for i, y in enumerate(m):
                a[s + i] = (a[s + i] - c * y) % p
        a.pop()
    while a and a[-1] % p == 0:
        a.pop()
    return a


def _polydiv(a, m, p):
    a = a[:]
    q = [0] * (len(a) - len(m) + 1)
    while len(a) >= len(m):
        c = a[-1] % p
        s = len(a) - len(m)
        q[s] = c
        for i, y in enumerate(m):
            a[s + i] = (a[s + i] - c * y) % p
        a.pop()
    return q


def _irreducibles(k, p):
    out = []
    for tail in itertools.product(range(p), repeat=k):
        P = list(tail) + [1]
        ok = True
        for j in range(1, k // 2 + 1):
            for g in _irreducibles_cache(j, p):
                if not _polymod(P, g, p):
                    ok = False
                    break
            if not ok:
                break
        if ok:
            out.append(P)
    return out


_cache = {}


def _irreducibles_cache(k, p):
    if (k, p) not in _cache:
        _cache[(k, p)] = _irreducibles(k, p)
    return _cache[(k, p)]


def slow_pattern(P, p):
    """Factor degrees by trial division; None if not squarefree."""
    n = len(P) - 1
    degs = []
    cur = P[:]
    for k in range(1, n + 1):
        for g in _irreducibles_cache(k, p):
            while len(cur) - 1 >= k and not _polymod(cur, g, p):
                cur = _polydiv(cur, g, p)
                degs.append(k)
        if len(cur) == 1:
            break
    if len(set(map(tuple, []))) != 0:
        pass
    return tuple(sorted(degs, reverse=True))


def slow_is_squarefree(P, p):
    n = len(P) - 1
    for k in range(1, n // 2 + 1):
        for g in _irreducibles_cache(k, p):
            g2 = fflib.poly_mul_py(g, g, p)
            if len(g2) - 1 <= n and not _polymod(P, g2, p):
                return False
    return True


def test_patterns_against_trial_division():
    for p, d, f in [(5, 2, [1, 0, 1]), (7, 2, [3, 1, 1]), (3, 2, [0, 0, 1]), (7, 3, [1, 2, 0, 1])]:
        codes = fflib.interval_codes(f, p)
        sq = fflib.f_squared(f, p)
        for idx in range(0, len(codes), max(1, len(codes) // 400)):
            P = sq[:]
            r = idx
            for i in range(d + 1):
                P[i] = (P[i] + r % p) % p
                r //= p
            fast = fflib.decode_code(int(codes[idx]), 2 * d)
            if slow_is_squarefree(P, p):
                assert fast == slow_pattern(P, p), (p, d, f, P, fast, slow_pattern(P, p))
            else:
                assert fast is None, (p, d, f, P, fast)


def test_small_counts():
    # d = 1: all monic quadratics, (q^2 - q)/2 irreducible
    for p in (3, 5, 7, 11, 13):
        assert fflib.count_irreducible([0, 1], p) == (p * p - p) // 2
    # the v5 example d=2, f=t^2+1, q=11
    assert fflib.count_irreducible([1, 0, 1], 11) == 330


def test_chebotarev_totals():
    # number of squarefree polynomials in I_f is q^{d+1} - (#discriminant points)
    p, f = 11, [2, 3, 0, 1]
    codes = fflib.interval_codes(f, p)
    cnt = Counter(fflib.decode_code(int(c), 6) for c in codes)
    assert sum(cnt.values()) == p ** 4


def test_twisted_identity_small():
    for p, f in [(5, [1, 0, 1]), (7, [2, 1, 1]), (7, [1, 2, 0, 1]), (5, [0, 1, 3, 1])]:
        d = len(f) - 1
        n_irr = fflib.count_irreducible(f, p)
        tw = fflib.twisted_count(f, p)
        assert tw["gen"] == 2 * d * n_irr, (p, f, tw, n_irr)


def test_closed_forms():
    for p in (5, 7, 11, 13):
        for f in ([0, 0, 0, 1], [1, 2, 3, 1], [0, 1, 0, 1], [4, 0, 1, 1]):
            assert fflib.count_irreducible(f, p) == fflib.closed_form(f, p), (p, f)
    for p in (3, 5, 7, 11):
        assert fflib.count_irreducible([1, 1, 1], p) == fflib.closed_form([1, 1, 1], p)


def test_partition_probabilities_sum_to_one():
    from fractions import Fraction

    def z(lam):
        c = Counter(lam)
        r = 1
        for k, m in c.items():
            r *= k ** m * factorial(m)
        return r

    def partitions(n, mx=None):
        if mx is None:
            mx = n
        if n == 0:
            yield ()
            return
        for k in range(min(n, mx), 0, -1):
            for rest in partitions(n - k, k):
                yield (k,) + rest

    assert sum(Fraction(1, z(l)) for l in partitions(6)) == 1


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
