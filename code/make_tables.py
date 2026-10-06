"""
Build the LaTeX tables in paper/tables/ and the summary in results/ from data/*.csv.
"""
import json
import math
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
    lines = [r"\begin{tabular}{rrrrrr}", r"\toprule",
             r"$d$ & classical & v5 (unproved) & Katz + Deligne & Cafure--Matera & via Sawin \\",
             r"\midrule"]
    for _, r in t[t.d.isin([1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40, 50])].iterrows():
        vals = {"katz": r.log10_katz_betti, "cm": r.log10_revised_rigorous}
        if str(r.log10_sawin) not in ("", "nan"):
            vals["sawin"] = float(r.log10_sawin)
        cell = {k: f"{v:.1f}" for k, v in vals.items()}
        classical = f"\\textbf{{{math.log10(r.d - 2):.2f}}}" if r.d >= 3 else r"\textbf{--}"
        lines.append(f"{int(r.d)} & {classical} & {r.log10_v5_claimed:.1f} & {cell['katz']} & {cell['cm']} & "
                     f"{cell.get('sawin', '--')} \\\\")
    lines += [r"\bottomrule",
              r"\multicolumn{6}{l}{\footnotesize Classical: $N_{\rm irr}\ge1$ for $q>d-2$ (Theorem~1); for $d\le2$ for every $q$.}",
              r"\end{tabular}"]
    return "\n".join(lines)


def open_table(s):
    lines = [r"\begin{tabular}{rrrrrrrr}", r"\toprule",
             r"$q$ & $d$ & intervals & orbits & $\min N_{\rm irr}$ & ratio & $\mathrm{Var}\,\Psi/q^{d+1}$ & $d-2$ \\",
             r"\midrule"]
    last = None
    for _, r in s.sort_values(["q", "d"]).iterrows():
        if r.q > r.d - 2 and r.q != 2:
            continue  # covered by Theorem 1 (keep q = 2 rows for completeness)
        if last is not None and r.q != last:
            lines.append(r"\addlinespace")
        last = r.q
        mark = r"$^\dagger$" if r["mode"] != "all" else ""
        var = f"{r['var_Lambda_over_q^(d+1)']:.2f}" if r["mode"] == "all" else "--"
        lines.append(f"{int(r.q)} & {int(r.d)}{mark} & {int(r.intervals_total)} & {int(r.orbits_examined)} & "
                     f"{int(r.min_N_irr)} & {r.min_ratio:.3f} & {var} & {int(r.d) - 2 if r.d >= 4 else '--'} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines)


def prime_powers(upto):
    out = []
    for n in range(2, upto + 1):
        p = next(k for k in range(2, n + 1) if n % k == 0)
        m = n
        while m % p == 0:
            m //= p
        if m == 1:
            out.append(n)
    return out


def verified_macros(s):
    full = s[s["mode"] == "all"]
    ok = {(int(r.q), int(r.d)) for _, r in full.iterrows() if r.min_N_irr >= 1}

    def maxd(q):
        """largest D such that every d with q <= d - 2 <= D - 2 is fully computed (contiguous)."""
        d = q + 1  # d <= q + 1 is covered by Theorem 1
        while (q, d + 1) in ok:
            d += 1
        return d

    ver = {q: maxd(q) for q in (2, 3, 4, 5)}
    # all-q bound: largest D such that every (q, d) with d <= D and q <= d - 2 is fully computed
    D = 2
    while True:
        dn = D + 1
        need = [q for q in prime_powers(max(dn - 2, 1)) if q <= dn - 2]
        if all((q, dn) in ok for q in need):
            D = dn
        else:
            break
        if D > 60:
            break
    open_rows = full[(full.q <= full.d - 2)]
    worst = open_rows.loc[open_rows.min_ratio.idxmin()]
    # heuristic maximum deficit over q >= 2, 3 <= d <= 200
    best = (0, None)
    for q in prime_powers(64):
        for d in range(3, 200):
            v = math.sqrt(2 * (d - 1) * (d - 2) * math.log(q) * math.exp(-(d + 1) * math.log(q)))
            if v > best[0]:
                best = (v, (q, d))
    lines = [
        f"\\newcommand{{\\VerTwo}}{{{ver[2]}}}",
        f"\\newcommand{{\\VerThree}}{{{ver[3]}}}",
        f"\\newcommand{{\\VerFour}}{{{ver[4]}}}",
        f"\\newcommand{{\\VerFive}}{{{ver[5]}}}",
        f"\\newcommand{{\\AllQd}}{{{D}}}",
        f"\\newcommand{{\\AllQdNext}}{{{D + 1}}}",
        f"\\newcommand{{\\MinRatioOpen}}{{{worst.min_ratio:.3f}}}",
        f"\\newcommand{{\\MinRatioCase}}{{({int(worst.q)},{int(worst.d)})}}",
        f"\\newcommand{{\\MaxDeficit}}{{{best[0]:.2f}}}",
        f"\\newcommand{{\\MaxDeficitCase}}{{({best[1][0]},{best[1][1]})}}",
    ]
    return "\n".join(lines) + "\n", {"ver": ver, "all_q_d": D, "min_ratio_open": float(worst.min_ratio),
                                      "min_ratio_case": [int(worst.q), int(worst.d)],
                                      "max_deficit": best[0], "max_deficit_case": list(best[1])}


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
    from open_data import load_summary
    open_summary = load_summary()
    if open_summary is not None:
        with open(os.path.join(TAB, "open_regime.tex"), "w") as fh:
            fh.write(open_table(open_summary))
        macros, open_info = verified_macros(open_summary)
        with open(os.path.join(TAB, "verified.tex"), "w") as fh:
            fh.write(macros)

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
    # classical Hayes-Weil bounds (Theorem 1) on the q > d rows
    cq = c.assign(lo=(c.p ** (c.d + 1) - c.d * c.p ** c.d) / (2 * c.d),
                  hi=(c.p ** (c.d + 1) + (c.d - 2).clip(lower=0) * c.p ** c.d) / (2 * c.d))
    m = cq.p > cq.d
    summary["classical_bound_check"] = {
        "rows_q_gt_d": int(m.sum()),
        "all_within_bounds": bool(((cq.N_irr >= cq.lo) & (cq.N_irr <= cq.hi))[m].all()),
    }
    if open_summary is not None:
        summary["open_range"] = open_info
    with open(os.path.join(RES, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    md = ["# Results summary", "",
          "Generated by `code/make_tables.py` from the CSV files in `data/`.", ""]
    if open_summary is not None:
        md += ["## Open range q <= d - 2 (`data/open_regime_summary.csv`)", "",
               f"* Every Legendre interval contains an irreducible polynomial for q = 2, d <= {open_info['ver'][2]}; "
               f"q = 3, d <= {open_info['ver'][3]}; q = 4, d <= {open_info['ver'][4]}; q = 5, d <= {open_info['ver'][5]}.",
               f"* With Theorem 1 (q >= d - 1), the analogue of Legendre's conjecture holds for **every q when d <= {open_info['all_q_d']}**.",
               f"* Smallest normalized count 2d N_irr / q^(d+1) in the open range: {open_info['min_ratio_open']:.3f} at (q, d) = {tuple(open_info['min_ratio_case'])}.",
               "", "| q | d | mode | intervals | orbits | min N_irr | min ratio | Var(Psi)/q^(d+1) | KR limit d-2 |",
               "|---|---|---|---|---|---|---|---|---|"]
        for _, r in open_summary.sort_values(["q", "d"]).iterrows():
            md.append(f"| {r.q} | {r.d} | {r['mode']} | {r.intervals_total} | {r.orbits_examined} | {r.min_N_irr} | "
                      f"{r.min_ratio:.4f} | {r['var_Lambda_over_q^(d+1)']:.3f} | {r.d - 2} |")
        md += [""]
    md += ["## Classical range check", "",
           f"* {summary['classical_bound_check']['rows_q_gt_d']} exact counts with q > d; all within the Hayes-Weil bounds of Theorem 1: "
           f"**{summary['classical_bound_check']['all_within_bounds']}**.", "",
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
           "log10 of the threshold q_0(d) beyond which N_irr >= 1 is proved. The classical Hayes-Weil bound "
           "(Theorem 1 of v7) gives q_0 = d - 2 and supersedes the other columns, which are kept from v6.", "",
           "| d | classical (v7) | v5 claim (unproved) | Katz + Deligne | Cafure-Matera | Sawin 2021 |", "|---|---|---|---|---|---|"]
    for _, r in t[t.d <= 12].iterrows():
        sw = "" if str(r.log10_sawin) in ("", "nan") else r.log10_sawin
        cl = f"{math.log10(r.d - 2):.3f}" if r.d >= 3 else "-"
        md.append(f"| {r.d} | {cl} | {r.log10_v5_claimed} | {r.log10_katz_betti} | {r.log10_revised_rigorous} | {sw} |")
    with open(os.path.join(RES, "summary.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(md) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
