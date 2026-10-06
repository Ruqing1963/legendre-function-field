"""
Experiment 4 -- explicit thresholds q_0(d) such that N_irr >= 1 for all q > q_0(d) (p > 2d).

  * v5_claimed      : (8d)^(2d+6), stated in the Zenodo v5 preprint.  NOT supported by a valid
                      argument (it relied on a Betti-number bound that is not in the cited source).
  * revised_rigorous: Corollary 5.2 of the revised paper, from Cafure--Matera, Thm 7.1, applied to the
                      twisted root variety W of dimension d+1 and degree delta = (d-1)!:
                      q_0 = max{ 2(d+2) delta^2 , (A + sqrt(B))^2 },  A = (delta-1)(delta-2), B = 5 delta^(13/3) + 2.
                      For d <= 3 the exact formulas of Section 6 give N_irr >= 1 for every admissible q.
  * katz_betti      : Proposition 5.4, Katz's explicit Betti bound + Deligne: q > 9 (d+1)^(6d-2).
  * sawin          : Proposition 5.6, from Sawin (Duke Math. J. 170 (2021), Cor. 4.7) with n = 2d, m = d-1:
                      |2d N_irr - q^(d+1)| <= C q^((d+2)/2) + d q + d^2,  C = 6 (2d+2)^(3d-1)  (p > 2d),
                      so N_irr >= 1 once q > (C+1)^(2/d).  Defined for d >= 2.

Output: data/thresholds.csv  (log10 values, since numbers are astronomically large)
"""
import csv
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")


def log10_v5(d):
    return (2 * d + 6) * math.log10(8 * d)


def log10_revised(d):
    lg_delta = math.lgamma(d) / math.log(10)  # log10((d-1)!)
    delta = 10 ** lg_delta if lg_delta < 300 else float("inf")
    if delta != float("inf"):
        A = (delta - 1) * (delta - 2)
        B = 5 * delta ** (13 / 3) + 2 if lg_delta * 13 / 3 < 300 else float("inf")
        if B != float("inf"):
            t1 = 2 * (d + 2) * delta ** 2
            t2 = (A + math.sqrt(B)) ** 2
            return math.log10(max(t1, t2))
    # huge-delta regime: A ~ delta^2, sqrt(B) ~ sqrt(5) delta^(13/6), so
    # (A + sqrt B)^2 ~ delta^4 (1 + sqrt(5) delta^(1/6))^2
    return 4 * lg_delta + 2 * math.log10(1 + math.sqrt(5) * 10 ** (lg_delta / 6))


def log10_katz(d):
    """Proposition 5.3: sigma_c(W) <= 3 (d+1)^(3d-1)  (Katz 2001, remark after Thm 12),
    Deligne => |#W - q^(d+1)| <= (sigma-1) q^(d+1/2); need s^2 - (sigma-1) s - 2 > 0, s = sqrt(q).
    Sufficient: q > sigma^2.  Returns log10(sigma^2) = 2 log10(3 (d+1)^(3d-1))."""
    return 2 * (math.log10(3) + (3 * d - 1) * math.log10(d + 1))


def log10_sawin(d):
    """q > (C + 1)^(2/d) with C = 6 (2d+2)^(3d-1); None for d = 1 (Sawin needs h < n)."""
    if d < 2:
        return None
    lgC = math.log10(6) + (3 * d - 1) * math.log10(2 * d + 2)
    return (2 / d) * (lgC + math.log10(1 + 10 ** (-lgC)))


def main():
    os.makedirs(DATA, exist_ok=True)
    rows = []
    for d in range(1, 61):
        rev = log10_revised(d)
        note = ""
        if d <= 3:
            note = "exact formula: N_irr >= 1 for all admissible q"
        rows.append({
            "d": d,
            "log10_v5_claimed": round(log10_v5(d), 4),
            "log10_katz_betti": round(log10_katz(d), 4),
            "log10_revised_rigorous": round(rev, 4),
            "log10_sawin": "" if log10_sawin(d) is None else round(log10_sawin(d), 4),
            "revised_q0_if_small": int(math.ceil(10 ** rev)) if rev < 15 else "",
            "note": note,
        })
    with open(os.path.join(DATA, "thresholds.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    for r in rows[:12]:
        print(r)


if __name__ == "__main__":
    main()
