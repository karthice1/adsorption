"""Supplementary Fig. S1: profile likelihoods of the reductant capacity (v0.2 designs D1, D3, D4).

    python scripts/make_supp_figs.py   -> manuscript/figures/figS1_profiles.png
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
d = pd.read_csv(ROOT / "results/v2/ident_profiles.csv")
studies = sorted(d.study_id.unique())
fig, axes = plt.subplots(1, len(studies), figsize=(4.2 * len(studies), 3.4), squeeze=False)
colors = {"D1": "#4C72B0", "D3": "#DD8452", "D4": "#55A868"}
for ax, sid in zip(axes[0], studies):
    s = d[d.study_id == sid]
    for des, g in s.groupby("design"):
        g = g.sort_values("log10_R0")
        ax.plot(g.log10_R0 - g.true_log10_R0, g.delta, "-o", ms=3, color=colors.get(des), label=f"{des} (n = {g.n_obs.iloc[0]})")
    ax.axhline(3.84, ls=":", c="0.4", lw=1)
    ax.axvline(0, ls="--", c="0.6", lw=0.8)
    ax.set_xlabel("log10 R0 − true value")
    ax.set_ylabel("Δ(−2 log L)")
    ax.set_title(f"{sid}: true reduction share {s.true_share.iloc[0]:.2f}", fontsize=10)
    ax.set_yscale("symlog", linthresh=10)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=8, frameon=False)
fig.tight_layout()
out = ROOT / "manuscript/figures/figS1_profiles.png"
fig.savefig(out, dpi=300)
print("written", out)
