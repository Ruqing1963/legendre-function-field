"""
Experiment 2 -- exhaustive search over ALL f of degree d for small primes p.

N_irr(f) depends on f only through the coefficients of t^{2d-1}, ..., t^{d+1} of f^2,
i.e. through (c_{d-1}, ..., c_1).  When p does not divide d, the substitution
t -> t - c_{d-1}/d preserves the count, so we may take c_{d-1} = 0.
We record min / max of N_irr over all classes, including small characteristic
p <= 2d, which lies outside the hypotheses of the theorems.

Output: data/exhaustive_min_counts.csv
"""
import csv
import itertools
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fflib  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")

CASES = [
    (2, [3, 5, 7, 11, 13]),
    (3, [3, 5, 7, 11, 13]),
    (4, [3, 5, 7, 11, 13]),
    (5, [3, 5, 7]),
    (6, [3, 5]),
]


def classes(d, p):
    if d % p != 0 and d >= 2:  # translation t -> t + b needs d invertible mod p
        # c_{d-1} = 0 ; free c_{d-2}..c_1
        for tail in itertools.product(range(p), repeat=d - 2):
            yield [0] + list(tail) + [0, 1]  # c_0 = 0, c_1..c_{d-2} = tail, c_{d-1} = 0
    else:
        for tail in itertools.product(range(p), repeat=d - 1):
            yield [0] + list(tail) + [1]


def main():
    os.makedirs(DATA, exist_ok=True)
    rows = []
    t0 = time.time()
    for d, primes in CASES:
        irr = fflib.irreducible_code(2 * d)
        for p in primes:
            vals = []
            argmin = None
            for f in classes(d, p):
                f = f[: d + 1]
                N = int((fflib.interval_codes(f, p) == irr).sum())
                vals.append(N)
                if argmin is None or N < argmin[0]:
                    argmin = (N, f)
            q = p
            rows.append({
                "d": d, "p": p, "p_gt_2d": int(p > 2 * d), "classes_checked": len(vals),
                "min_N_irr": min(vals), "max_N_irr": max(vals),
                "main_term": round(q ** (d + 1) / (2 * d), 3),
                "min_ratio": round(min(vals) * 2 * d / q ** (d + 1), 6),
                "argmin_f_low_first": " ".join(map(str, argmin[1])),
                "always_positive": int(min(vals) > 0),
            })
            print(f"d={d} p={p}: classes={len(vals)} min={min(vals)} max={max(vals)} ({time.time()-t0:.1f}s)",
                  flush=True)
    with open(os.path.join(DATA, "exhaustive_min_counts.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
