"""
Build every figure in figures/ from the CSV files in data/.

fig1_ratio.pdf/png                 2d*N_irr / q^(d+1) -> 1  (all d on one axis)
fig2_normalized_error.pdf/png      (2d*N_irr - q^(d+1)) / q^d, one panel per d (p > 2d)
fig3_factorization_types.pdf/png   observed factorisation-type frequency vs 1/z_lambda
fig4_thresholds.pdf/png            log10 q_0(d): v5 claim vs the two rigorous bounds
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
FIG = os.path.join(ROOT, "figures")

# validated categorical palette (fixed order) + chrome
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
INK, INK2, MUTED, GRID, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Segoe UI", "DejaVu Sans", "Arial"],
    "font.size": 9,
    "axes.edgecolor": AXIS,
    "axes.labelcolor": INK2,
    "axes.titlecolor": INK,
    "axes.titlesize": 10,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "legend.frameon": False,
    "legend.labelcolor": INK2,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
})


def save(fig, name):
    os.makedirs(FIG, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"))
    plt.close(fig)


def fig1(counts):
    c = counts[counts.p_gt_2d == 1]
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    for i, d in enumerate(sorted(c.d.unique())):
        sub = c[c.d == d].groupby("p")["ratio_2dN_over_q^(d+1)"].mean().reset_index()
        ax.plot(sub.p, sub.iloc[:, 1], "-o", color=SERIES[i], lw=2, ms=4, label=f"d = {d}")
    ax.axhline(1, color=AXIS, lw=1, ls="--")
    ax.set_xscale("log")
    ax.set_xlabel("q = p  (prime, p > 2d)")
    ax.set_ylabel(r"$2d\,N_{\rm irr}\,/\,q^{d+1}$")
    ax.set_title("Irreducibles in the Legendre interval approach the Chebotarev density 1/(2d)", loc="left")
    ax.legend(ncol=3, fontsize=8, loc="lower right")
    save(fig, "fig1_ratio")


def fig2(counts):
    c = counts[counts.p_gt_2d == 1]
    ds = sorted(c.d.unique())
    fig, axes = plt.subplots(2, 3, figsize=(8.4, 4.8), sharex=False)
    for ax, d in zip(axes.flat, ds):
        sub = c[c.d == d]
        ax.scatter(sub.p, sub["err_over_q^d"], s=22, color=SERIES[0], edgecolor=SURFACE, linewidth=1, zorder=3)
        ax.axhline(0, color=AXIS, lw=1)
        ax.set_title(f"d = {d}", loc="left")
        ax.set_xlabel("q")
        lim = max(0.25, 1.15 * np.abs(sub["err_over_q^d"]).max())
        ax.set_ylim(-lim, lim)
    for ax in axes.flat[len(ds):]:
        ax.axis("off")
    axes[0, 0].set_ylabel(r"$(2d\,N_{\rm irr}-q^{d+1})/q^{d}$")
    axes[1, 0].set_ylabel(r"$(2d\,N_{\rm irr}-q^{d+1})/q^{d}$")
    fig.suptitle("Normalized error stays bounded (each point is one f)", x=0.01, ha="left", color=INK, fontsize=10)
    fig.tight_layout()
    save(fig, "fig2_normalized_error")


def fig3(ft):
    groups = list(ft.groupby(["d", "p"]))
    ncol = 3
    nrow = (len(groups) + ncol - 1) // ncol
    fig, axs = plt.subplots(nrow, ncol, figsize=(8.4, 2.9 * nrow), squeeze=False)
    axes = list(axs.flat)
    for ax in axes[len(groups):]:
        ax.axis("off")
    for ax, ((d, p), sub) in zip(axes, groups):
        x = sub.S_2d_probability.values
        y = sub.frequency.values
        mask = y > 0
        lo = min(x.min(), y[mask].min()) / 2
        ax.plot([lo, 1], [lo, 1], color=AXIS, lw=1, ls="--", zorder=1)
        ax.scatter(x[mask], y[mask], s=14, color=SERIES[0], edgecolor=SURFACE, linewidth=0.8, zorder=3)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(f"d={d}, q={p}", loc="left")
        ax.set_xlabel(r"$1/z_\lambda$ in $S_{%d}$" % (2 * d))
        irr = sub[sub.partition == str(2 * d)]
        if len(irr):
            ax.scatter(irr.S_2d_probability, irr.frequency, s=46, facecolor="none",
                       edgecolor=SERIES[1], linewidth=1.5, zorder=4)
            ax.annotate("irreducible", (irr.S_2d_probability.iloc[0], irr.frequency.iloc[0]),
                        xytext=(7, -12), textcoords="offset points", ha="left", fontsize=7, color=INK2)
    for r in range(nrow):
        axes[r * ncol].set_ylabel("observed frequency")
    fig.suptitle("Factorization types in $\\mathcal{I}_{t^d}$ follow the cycle types of $S_{2d}$",
                 x=0.01, ha="left", color=INK, fontsize=10)
    fig.tight_layout()
    save(fig, "fig3_factorization_types")


def fig4(th):
    th = th[th.d <= 40]
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    ax.plot(th.d, th.log10_v5_claimed, "--", color=SERIES[1], lw=2, label="v5 claim $(8d)^{2d+6}$ (unproved)")
    ax.plot(th.d, th.log10_katz_betti, "-", color=SERIES[2], lw=2, label="Katz Betti bound + Deligne (Prop. 5.4)")
    ax.plot(th.d, th.log10_revised_rigorous, "-", color=SERIES[0], lw=2, label="Cafureâ€“Matera (Cor. 5.3)")
    ax.axvspan(0.5, 3.5, color=GRID, alpha=0.6, lw=0)
    ax.text(2, ax.get_ylim()[1] * 0.45, "exact\n(Thm D)", ha="center", va="top", fontsize=8, color=INK2)
    ax.set_xlabel("d")
    ax.set_ylabel(r"$\log_{10} q_0(d)$")
    ax.set_title("Thresholds beyond which every Legendre interval contains an irreducible", loc="left")
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(0.13, 1.0))
    save(fig, "fig4_thresholds")


def main():
    counts = pd.read_csv(os.path.join(DATA, "interval_counts.csv"))
    ft = pd.read_csv(os.path.join(DATA, "factorization_types.csv"), dtype={"partition": str})
    th = pd.read_csv(os.path.join(DATA, "thresholds.csv"))
    fig1(counts)
    fig2(counts)
    fig3(ft)
    fig4(th)
    print("figures written to", FIG)


if __name__ == "__main__":
    main()
