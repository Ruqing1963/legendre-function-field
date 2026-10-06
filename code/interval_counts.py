"""
Experiment 1 -- exact counts N_irr(f, q) and full factorisation statistics on I_f.

Outputs
  data/interval_counts.csv          one row per (d, p, f)
  data/factorization_types.csv      empirical frequency of each factorisation type vs P(lambda) = 1/z_lambda
"""
import csv
import os
import random
import sys
import time
from collections import Counter
from math import factorial

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fflib  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")

# (d, list of primes) -- all with p > 2d, plus a few small-characteristic cases (p <= 2d) flagged as such
GRID = {
    1: [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97],
    2: [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47],
    3: [3, 5, 7, 11, 13, 17, 19, 23, 29, 31],
    4: [3, 5, 7, 11, 13, 17, 19],
    5: [3, 7, 11, 13],
    6: [5, 7, 13],
}
FACTOR_STATS = [(2, 47), (3, 31), (4, 19), (5, 13), (6, 13)]
SEED = 20260220


def z_lambda(lam):
    c = Counter(lam)
    r = 1
    for k, m in c.items():
        r *= k ** m * factorial(m)
    return r


def sample_fs(d, p, rng):
    """Three choices of f: t^d, a fixed 'structured' f, one random f."""
    fs = [[0] * d + [1]]
    structured = [1] + [0] * (d - 1) + [1] if d >= 2 else [1, 1]  # t^d + 1
    if structured not in fs:
        fs.append(structured)
    fs.append([rng.randrange(p) for _ in range(d)] + [1])
    return fs


def main():
    os.makedirs(DATA, exist_ok=True)
    rng = random.Random(SEED)
    rows = []
    stats_rows = []
    t0 = time.time()
    for d, primes in GRID.items():
        n = 2 * d
        irr = fflib.irreducible_code(n)
        for p in primes:
            for f in sample_fs(d, p, rng):
                codes = fflib.interval_codes(f, p)
                q = p
                N = int((codes == irr).sum())
                nonsq = int((codes < 0).sum())
                main_term = q ** (d + 1) / n
                cf = fflib.closed_form(f, p)
                rows.append({
                    "d": d, "p": p, "f_coeffs_low_first": " ".join(map(str, f)),
                    "p_gt_2d": int(p > 2 * d),
                    "size_I_f": q ** (d + 1), "N_irr": N, "main_term": round(main_term, 4),
                    "ratio_2dN_over_q^(d+1)": round(n * N / q ** (d + 1), 8),
                    "err_over_q^d": round((n * N - q ** (d + 1)) / q ** d, 6),
                    "err_over_q^(d+1/2)": round((n * N - q ** (d + 1)) / q ** (d + 0.5), 6),
                    "nonsquarefree": nonsq,
                    "closed_form": "" if cf is None else cf,
                    "closed_form_match": "" if cf is None else int(cf == N),
                })
                if (d, p) in FACTOR_STATS and f == [0] * d + [1]:
                    cnt = Counter(int(c) for c in codes if c >= 0)
                    tot = sum(cnt.values())
                    seen = {}
                    for code, k in cnt.items():
                        seen[fflib.decode_code(code, n)] = k
                    # list every partition of n (include those never observed)
                    def parts(m, mx=None):
                        mx = m if mx is None else mx
                        if m == 0:
                            yield ()
                            return
                        for k in range(min(m, mx), 0, -1):
                            for rest in parts(m - k, k):
                                yield (k,) + rest
                    for lam in parts(n):
                        k = seen.get(lam, 0)
                        stats_rows.append({
                            "d": d, "p": p, "f_coeffs_low_first": " ".join(map(str, f)),
                            "partition": "+".join(map(str, lam)), "count": k,
                            "frequency": round(k / tot, 8),
                            "S_2d_probability": round(1 / z_lambda(lam), 8),
                            "squarefree_total": tot,
                        })
            print(f"d={d} p={p} done  ({time.time() - t0:.1f}s)", flush=True)

    with open(os.path.join(DATA, "interval_counts.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(DATA, "factorization_types.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(stats_rows[0].keys()))
        w.writeheader()
        w.writerows(stats_rows)
    print(f"wrote {len(rows)} count rows, {len(stats_rows)} factorisation rows")


if __name__ == "__main__":
    main()
