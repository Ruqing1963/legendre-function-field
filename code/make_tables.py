"""
Build the LaTeX tables in paper/tables/ and the summary in results/ from data/*.csv.
"""
import json
import os

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
TAB = os.path.join(ROOT, "paper", "tables")
RES = os.path.join(ROOT, "results")


def fpoly(coeffs):
    c = list(map(int, str(coeffs).split()))
    terms = []
    for i in range(len(c) - 1, -1, -1):
        a = c[i]
        if a == 0:
            continue
        mon = "" if i == 0 else ("t" if i == 1 else f"t^{{{i}}}")
        coef = "" if (a == 1 and i > 0) else str(a)
        terms.append(coef + mon)
    return "$" + "+".join(terms) + "$"


def counts_table(c):
    rows = []
    pick = c[(c.p_gt_2d == 1)]
    for d in sorted(pick.d.unique()):
        ps = sorted(pick[pick.d == d].p.unique())
        chosen = [ps[0], ps[len(ps) // 2], ps[-1]] if len(ps) >= 3 else ps
        for p in sorted(set(chosen)):
            sub = pick[(pick.d == d) & (pick.p == p)]
            for _, r in sub.iloc[[0, -1]].drop_duplicates().iterrows():
                rows.append(r)
    lines = [r"\begin{tabular}{rrlrrrl}", r"\toprule",
             r"$d$ & $q$ & $f$ & $N_{\rm irr}$ & $q^{d+1}/2d$ & $\frac{2dN-q^{d+1}}{q^d}$ & formula \\", r"\midrule"]
    last_d = None
    for r in rows:
        if last_d is not None and r.d != last_d:
            lines.append(r"\addlinespace")
        last_d = r.d
        cf = "" if pd.isna(r.closed_form) or r.closed_form == "" else (r"\checkmark" if int(r.closed_form_match) == 1 else "MISMATCH")
        lines.append(f"{int(r.d)} & {int(r.p)} & {fpoly(r.f_coeffs_low_first)} & {int(r.N_irr)} & "
                     f"{r.main_term:.1f} & {r['err_over_q^d']:+.3f} & {cf} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines)


def exhaustive_table(e):
    lines = [r"\begin{tabular}{rrcrrrr}", r"\toprule",
             r"$d$ & $q$ & $p>2d$ & classes of $f$ & $\min N_{\rm irr}$ & $\max N_{\rm irr}$ & $q^{d+1}/2d$ \\",
             r"\midrule"]
    last = None
    for _, r in e.iterrows():
        if last is not None and r.d != last:
            lines.append(r"\addlinespace")
        last = r.d
        lines.append(f"{int(r.d)} & {int(r.p)} & {'yes' if r.p_gt_2d else 'no'} & {int(r.classes_checked)} & "
                     f"{int(r.min_N_irr)} & {int(r.max_N_irr)} & {r.main_term:.1f} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines)


def thresholds_table(t):
    lines = [r"\begin{tabular}{rrrrr}", r"\toprule",
             r"$d$ & v5 claim (unproved) & Prop.~5.4 (Katz) & Cor.~5.3 (Cafure--Matera) & Prop.~5.5 (Sawin) \\",
             r"\midrule"]
    for _, r in t[t.d.isin([1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40, 50])].iterrows():
        vals = {"katz": r.log10_katz_betti, "cm": r.log10_revised_rigorous}
        if str(r.log10_sawin) not in ("", "nan"):
            vals["sawin"] = float(r.log10_sawin)
        best = min(vals, key=vals.get)
        cell = {k: (rf"\textbf{{{v:.1f}}}" if (k == best and r.d > 3) else f"{v:.1f}") for k, v in vals.items()}
        star = r"$^{\ast}$" if r.d <= 3 else ""
        lines.append(f"{int(r.d)} & {r.log10_v5_claimed:.1f} & {cell['katz']}{star} & {cell['cm']}{star} & "
                     f"{cell.get('sawin', '--')}{star if 'sawin' in cell else ''} \\\\")
    lines += [r"\bottomrule",
              r"\multicolumn{5}{l}{\footnotesize $^{\ast}$\,not needed: Theorem~D gives $N_{\rm irr}\ge1$ for all admissible $q$.}",
              r"\end{tabular}"]
    return "\n".join(lines)


def main():
    os.makedirs(TAB, exist_ok=True)
    os.makedirs(RES, exist_ok=True)
    c = pd.read_csv(os.path.join(DATA, "interval_counts.csv"))
    e = pd.read_csv(os.path.join(DATA, "exhaustive_min_counts.csv"))
    t = pd.read_csv(os.path.join(DATA, "thresholds.csv"))
    tw = pd.read_csv(os.path.join(DATA, "twisted_identity.csv"))
    ft = pd.read_csv(os.path.join(DATA, "factorization_types.csv"), dtype={"partition": str})
    with open(os.path.join(TAB, "counts.tex"), "w") as fh:
        fh.write(counts_table(c))
    with open(os.path.join(TAB, "exhaustive.tex"), "w") as fh:
        fh.write(exhaustive_table(e))
    with open(os.path.join(TAB, "thresholds.tex"), "w") as fh:
        fh.write(thresholds_table(t))

    cf = c[c.closed_form.notna() & (c.closed_form.astype(str) != "")]
    ft = ft.assign(absdiff=(ft.frequency - ft.S_2d_probability).abs())
    tv = 0.5 * ft.groupby(["d", "p"])["absdiff"].sum()
    summary = {
        "interval_counts": {
            "rows": int(len(c)),
            "d_range": [int(c.d.min()), int(c.d.max())],
            "min_N_irr": int(c.N_irr.min()),
            "closed_form_checked": int(len(cf)),
            "closed_form_mismatches": int((cf.closed_form_match.astype(int) == 0).sum()),
            "max_abs_err_over_q^d_by_d (p>2d)": {
                int(d): float(g["err_over_q^d"].abs().max()) for d, g in c[c.p_gt_2d == 1].groupby("d")},
        },
        "twisted_identity": {
            "rows": int(len(tw)),
            "identity_holds_all": bool(tw.identity_holds.all()),
            "E_within_bound_all": bool(tw.E_within_bound.all()),
            "within_CM_bound_all": bool(tw.within_CM_bound.all()),
        },
        "exhaustive": {
            "cases": int(len(e)),
            "all_positive": bool(e.always_positive.all()),
            "smallest_min": int(e.min_N_irr.min()),
        },
        "factorization_types_total_variation_distance": {f"d={d},p={p}": round(float(v), 5) for (d, p), v in tv.items()},
    }
    with open(os.path.join(RES, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    md = ["# Results summary", "",
          "Generated by `code/make_tables.py` from the CSV files in `data/`.", "",
          "## Exact counts (`data/interval_counts.csv`)", "",
          f"* {summary['interval_counts']['rows']} triples (d, q, f) with d = {summary['interval_counts']['d_range'][0]}..{summary['interval_counts']['d_range'][1]}.",
          f"* Smallest N_irr observed: **{summary['interval_counts']['min_N_irr']}** (never zero).",
          f"* Theorem D closed forms checked on {summary['interval_counts']['closed_form_checked']} rows, mismatches: **{summary['interval_counts']['closed_form_mismatches']}**.",
          "", "| d | max abs (2d N - q^(d+1)) / q^d, p > 2d |", "|---|---|"]
    for d, v in summary["interval_counts"]["max_abs_err_over_q^d_by_d (p>2d)"].items():
        md.append(f"| {d} | {v:.4f} |")
    md += ["", "## Twisted identity (`data/twisted_identity.csv`)", "",
           f"* {summary['twisted_identity']['rows']} cases; #W(F_q) = 2d N_irr + E in all: **{summary['twisted_identity']['identity_holds_all']}**.",
           f"* E within the bound of Theorem B: **{summary['twisted_identity']['E_within_bound_all']}**; |#W - q^(d+1)| within Cafure-Matera bound: **{summary['twisted_identity']['within_CM_bound_all']}**.",
           "", "## Exhaustive search over f (`data/exhaustive_min_counts.csv`)", "",
           "| d | q | p > 2d | classes | min N_irr | max N_irr | q^(d+1)/2d |", "|---|---|---|---|---|---|---|"]
    for _, r in e.iterrows():
        md.append(f"| {r.d} | {r.p} | {'yes' if r.p_gt_2d else 'no'} | {r.classes_checked} | {r.min_N_irr} | {r.max_N_irr} | {r.main_term} |")
    md += ["", "## Factorization types vs S_2d (`data/factorization_types.csv`)", "",
           "Total variation distance between observed type frequencies and 1/z_lambda:", "",
           "| case | TV distance |", "|---|---|"]
    for k, v in summary["factorization_types_total_variation_distance"].items():
        md.append(f"| {k} | {v} |")
    md += ["", "## Thresholds (`data/thresholds.csv`)", "",
           "log10 of the threshold q_0(d) beyond which N_irr >= 1. For d <= 3 no threshold is needed (Theorem D).", "",
           "| d | v5 claim (unproved) | Katz + Deligne | Cafure-Matera | Sawin 2021 |", "|---|---|---|---|---|"]
    for _, r in t[t.d <= 12].iterrows():
        sw = "" if str(r.log10_sawin) in ("", "nan") else r.log10_sawin
        md.append(f"| {r.d} | {r.log10_v5_claimed} | {r.log10_katz_betti} | {r.log10_revised_rigorous} | {sw} |")
    with open(os.path.join(RES, "summary.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(md) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
