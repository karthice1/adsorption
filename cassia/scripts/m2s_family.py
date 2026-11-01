"""M2s (exact study-level split conformal) by material family at nominal 90% (Table S5b).

    python scripts/m2s_family.py  -> results/v3/T_m2s_family_full.csv
Source: results/v3/m2s_per_study.csv (run_v3.py stage m2s; 10 corpora, outer leave-one-study-out).
Unbounded intervals count as covering every observation; width is averaged over bounded intervals.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm


def wilson(k, n, a=0.05):
    z = norm.ppf(1 - a / 2); p = k / n; c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return max(0.0, c - h), min(1.0, c + h)

ROOT = Path(__file__).resolve().parents[1]
d = pd.read_csv(ROOT / "results/v3/m2s_per_study.csv")
d = d[np.isclose(d.alpha, 0.10)]
rows = []
for (tgt, fam), g in d.groupby(["target", "family"]):
    b = g[~g.unbounded]
    k = int((g.coverage >= 0.9 - 1e-9).sum()); lo, hi = wilson(k, len(g))
    rows.append({"target": tgt, "family": fam, "n_studies": len(g), "mean_coverage": g.coverage.mean(),
                 "frac_studies_cov_ge_0.9": k / len(g), "wilson95_lo": lo, "wilson95_hi": hi, "frac_unbounded": g.unbounded.mean(),
                 "mean_width_log10_bounded": b.width_log10.mean() if len(b) else np.nan})
T = pd.DataFrame(rows).round(3)
T.to_csv(ROOT / "results/v3/T_m2s_family_full.csv", index=False)
print(T.to_string())
