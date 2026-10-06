"""
Experiment 3 -- verify the twisted-variety identity (Theorem 4.3 of the revised paper)

    #W(F_q) = 2d * N_irr + E,   0 <= E <= sum_{e | 2d, e < 2d} q^e,

by enumerating every beta in F_{q^{2d}} and testing whether its characteristic
polynomial lies in I_f, and compare |#W(F_q) - q^{d+1}| with the explicit
Cafure--Matera bound (delta-1)(delta-2) q^{d+1/2} + 5 delta^{13/3} q^d, delta = (d-1)!.

Output: data/twisted_identity.csv
"""
import csv
import os
import random
import sys
import time
from math import factorial

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fflib  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")

CASES = [(1, p) for p in (3, 5, 7, 11)] + \
        [(2, p) for p in (5, 7, 11, 13)] + \
        [(3, p) for p in (5, 7, 11)] + \
        [(4, 11)]
SEED = 7


def main():
    os.makedirs(DATA, exist_ok=True)
    rng = random.Random(SEED)
    rows = []
    t0 = time.time()
    for d, p in CASES:
        fs = [[0] * d + [1], [rng.randrange(p) for _ in range(d)] + [1]]
        m = fflib.find_irreducible(2 * d, p)
        for f in fs:
            q = p
            N = fflib.count_irreducible(f, p)
            tw = fflib.twisted_count(f, p, m=m)
            delta = factorial(d - 1)
            cm_bound = (delta - 1) * (delta - 2) * q ** (d + 0.5) + 5 * delta ** (13 / 3) * q ** d
            e_bound = sum(q ** e for e in range(1, 2 * d) if (2 * d) % e == 0)
            rows.append({
                "d": d, "p": p, "f_coeffs_low_first": " ".join(map(str, f)),
                "W_points": tw["W"], "generators": tw["gen"], "E_subfield": tw["E"],
                "N_irr": N, "identity_holds": int(tw["gen"] == 2 * d * N),
                "E_bound": e_bound, "E_within_bound": int(0 <= tw["E"] <= e_bound),
                "W_minus_q^(d+1)": tw["W"] - q ** (d + 1),
                "CM_bound": round(cm_bound, 2),
                "CM_hypothesis_q>2(d+2)delta^2": int(q > 2 * (d + 2) * delta ** 2),
                "within_CM_bound": int(abs(tw["W"] - q ** (d + 1)) <= cm_bound),
            })
            print(rows[-1], f"({time.time()-t0:.1f}s)", flush=True)
    with open(os.path.join(DATA, "twisted_identity.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
