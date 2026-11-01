"""True partition of Cr(VI) disappearance by scenario, and Monte Carlo precision of the benchmark.

    python scripts/bench_partition.py  -> results/v3/T_bench_partition.csv, T_bench_mcse.csv
Uses results/v2/bench_studies.csv.gz (5 scenarios x 100 replicate corpora x 30 studies).
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
d = pd.read_csv(ROOT / "results/v2/bench_studies.csv.gz")
cols = {"share_ads": "adsorbed Cr(VI)", "share_red": "reduced (total)", "share_cr3_aq": "dissolved Cr(III)",
        "share_cr3_surface": "Cr(III) on solid", "share_precip": "precipitated Cr(III)"}
rows = []
for sc, g in d.groupby("scenario"):
    for c, lab in cols.items():
        x = g[c].dropna()
        rows.append({"scenario": sc, "component": lab, "median": x.median(), "p05": x.quantile(0.05), "p95": x.quantile(0.95), "mean": x.mean()})
P = pd.DataFrame(rows).round(3)
P.to_csv(ROOT / "results/v3/T_bench_partition.csv", index=False)

# Monte Carlo precision: per-replicate estimates -> between-replicate SD -> MCSE at R replicates,
# and replicates needed for a 95% half-width of 0.01 (proportions) or 0.05 (correlation).
d["fit95"] = d.R2_best >= 0.95
rows = []
for sc, g in d.groupby("scenario"):
    per = g.groupby("rep").agg(fit95=("fit95", "mean"), share_red=("share_red", "mean"))
    hi = g[g.share_red > 0.2].groupby("rep").fit95.mean()
    rho = g.groupby("rep").apply(lambda x: x[["share_red", "qmax_ratio"]].corr(method="spearman").iloc[0, 1], include_groups=False)
    for name, v, hw in (("P(R2>=0.95)", per.fit95, 0.01), ("P(R2>=0.95 | share_red>0.2)", hi, 0.02),
                        ("mean reduction share", per.share_red, 0.01), ("Spearman(share_red, qmax ratio)", rho, 0.05)):
        v = v.dropna(); R = len(v)
        if R < 2 or not np.isfinite(v.std(ddof=1)) or v.std(ddof=1) == 0:
            continue
        sd = v.std(ddof=1)
        rows.append({"scenario": sc, "quantity": name, "estimate": v.mean(), "between_rep_SD": sd, "R": R,
                     "MCSE": sd / np.sqrt(R), "target_95_halfwidth": hw, "R_needed": int(np.ceil((1.96 * sd / hw) ** 2))})
M = pd.DataFrame(rows).round(4)
M.to_csv(ROOT / "results/v3/T_bench_mcse.csv", index=False)
print(P.to_string()); print(M.to_string())
