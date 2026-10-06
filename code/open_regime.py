"""
Experiment 5 -- the open range q <= d.

For each (q, d) in the grid, every Legendre interval I_f over F_q is examined, up to the
affine substitutions t -> lambda t + b (orbits from fqlib.legendre_orbits). For each orbit
representative we record N_irr, the prime-power correction E, and
    Lambda_sum = 2d N_irr + E  ( = sum of the von Mangoldt function over I_f ).
Cases marked "sampled" examine a random subset of orbits only.

Outputs
  data/open_regime_orbits.csv    one row per (q, d, orbit)
  data/open_regime_summary.csv   one row per (q, d)

Usage: python code/open_regime.py [--only q,d] [--max-tests N]
"""
import argparse
import csv
import math
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fqlib  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")

# (q, d, mode); mode "all" = every orbit, ("sample", k) = k random orbits
GRID = (
    [(2, d, "all") for d in range(2, 21)]
    + [(3, d, "all") for d in range(3, 11)]
    + [(4, d, "all") for d in range(4, 11)]
    + [(5, d, "all") for d in range(5, 8)]
    + [(3, 11, ("sample", 30)), (5, 8, ("sample", 20)), (7, 7, ("sample", 40)), (8, 8, ("sample", 6))]
)
SEED = 20261006


def summarize(q, d, rows, mode, n_intervals_total):
    n = 2 * d
    w = [r["orbit_size"] for r in rows]
    W = sum(w)
    lam = [r["Lambda_sum"] for r in rows]
    mean = sum(wi * x for wi, x in zip(w, lam)) / W
    var = sum(wi * (x - q ** (d + 1)) ** 2 for wi, x in zip(w, lam)) / W
    nmin = min(r["N_irr"] for r in rows)
    nmax = max(r["N_irr"] for r in rows)
    return {
        "q": q, "d": d, "p": fqlib.char_of(q), "mode": "all" if mode == "all" else f"sample{mode[1]}",
        "intervals_total": n_intervals_total, "orbits_examined": len(rows), "intervals_examined": W,
        "min_N_irr": nmin, "max_N_irr": nmax,
        "main_term": q ** (d + 1) / n,
        "min_ratio": nmin * n / q ** (d + 1),
        "max_ratio": nmax * n / q ** (d + 1),
        "mean_Lambda_over_q^(d+1)": mean / q ** (d + 1),
        "var_Lambda_over_q^(d+1)": var / q ** (d + 1),
        "KR_prediction_d-2": d - 2,
        "min_Z": min((x - q ** (d + 1)) / math.sqrt(max(d - 2, 1) * q ** (d + 1)) for x in lam),
        "classical_bound_applies(q>d)": int(q > d),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None)
    ap.add_argument("--max-tests", type=float, default=1.5e9)
    ap.add_argument("--full", action="store_true", help="with --only: examine every orbit")
    args = ap.parse_args()
    os.makedirs(DATA, exist_ok=True)
    rng = random.Random(SEED)
    orbit_rows, summary_rows = [], []
    grid = GRID
    if args.only:
        q0, d0 = map(int, args.only.split(","))
        grid = [g for g in GRID if g[0] == q0 and g[1] == d0] or [(q0, d0, "all")]
        if args.full:
            grid = [(q0, d0, "all")]
    t0 = time.time()
    for q, d, mode in grid:
        if q ** (d - 1) > 2e6:
            print(f"skip q={q} d={d}: too many f to enumerate", flush=True)
            continue
        tabs = fqlib.field_tables(q)
        orbits = fqlib.legendre_orbits(q, d)
        n_int = sum(s for _, s in orbits)
        chosen = orbits if mode == "all" else rng.sample(orbits, min(mode[1], len(orbits)))
        cost = len(chosen) * q ** (d + 1)
        if cost > args.max_tests:
            print(f"skip q={q} d={d}: {cost:.2e} tests exceeds budget", flush=True)
            continue
        rows = []
        for f, size in chosen:
            N = fqlib.interval_count(f, q, tabs)
            E = fqlib.prime_power_E(f, q, tabs)
            rows.append({"q": q, "d": d, "f_coeffs_low_first": " ".join(map(str, f)), "orbit_size": size,
                         "N_irr": N, "E": E, "Lambda_sum": 2 * d * N + E})
        orbit_rows += rows
        s = summarize(q, d, rows, mode, n_int)
        summary_rows.append(s)
        print(f"q={q} d={d} {s['mode']}: orbits={len(rows)} min_N={s['min_N_irr']} min_ratio={s['min_ratio']:.4f} "
              f"var/q^(d+1)={s['var_Lambda_over_q^(d+1)']:.3f} (KR {d-2})  [{time.time()-t0:.0f}s]", flush=True)

    suffix = "" if not args.only else "_" + args.only.replace(",", "_")
    with open(os.path.join(DATA, f"open_regime_orbits{suffix}.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(orbit_rows[0].keys()))
        w.writeheader()
        w.writerows(orbit_rows)
    with open(os.path.join(DATA, f"open_regime_summary{suffix}.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()))
        w.writeheader()
        for r in summary_rows:
            w.writerow({k: (round(v, 6) if isinstance(v, float) else v) for k, v in r.items()})


if __name__ == "__main__":
    main()
