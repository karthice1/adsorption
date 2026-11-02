"""M2s (exact study-level split conformal) by material family at nominal 90% (Table S5b).

    python scripts/m2s_family.py  -> results/v3/T_m2s_family_full.csv
Source: results/v3/m2s_per_study.csv (run_v3.py stage m2s; 10 corpora, outer leave-one-study-out).
Unbounded intervals count as covering every observation; width is averaged over bounded intervals.

Uncertainty of the fraction of studies with coverage >= 0.9, per family:
  Wilson 95 %   : descriptive, treats the studies of a family as independent Bernoulli observations;
  cluster 95 %  : primary; corpus-cluster percentile bootstrap (the 10 corpora resampled with replacement,
                  B = 4000, seed 11), because studies in one corpus share the redrawn family hyper-priors;
  envelope 95 % : the wider of the two (lower = min, upper = max), reported in Table S5b, because a percentile
                  bootstrap over only 10 clusters can itself be too narrow and collapses to a point when every
                  study in the family covers.
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
B = 4000
rng = np.random.default_rng(11)
corpora = np.sort(d.corpus.unique())
boot_idx = rng.integers(0, len(corpora), (B, len(corpora)))
rows = []
for (tgt, fam), g in d.groupby(["target", "family"]):
    b = g[~g.unbounded]
    k = int((g.coverage >= 0.9 - 1e-9).sum()); lo, hi = wilson(k, len(g))
    hit = (g.coverage >= 0.9 - 1e-9).groupby(g.corpus).agg(["sum", "size"]).reindex(corpora, fill_value=0)
    ks, ns = hit["sum"].values[boot_idx].sum(1), hit["size"].values[boot_idx].sum(1)
    fr = ks[ns > 0] / ns[ns > 0]
    rows.append({"target": tgt, "family": fam, "n_studies": len(g), "mean_coverage": g.coverage.mean(),
                 "frac_studies_cov_ge_0.9": k / len(g), "wilson95_lo": lo, "wilson95_hi": hi,
                 "cluster95_lo": float(np.quantile(fr, 0.025)), "cluster95_hi": float(np.quantile(fr, 0.975)),
                 "envelope95_lo": min(lo, float(np.quantile(fr, 0.025))), "envelope95_hi": max(hi, float(np.quantile(fr, 0.975))),
                 "n_corpora": int(g.corpus.nunique()), "frac_unbounded": g.unbounded.mean(),
                 "mean_width_log10_bounded": b.width_log10.mean() if len(b) else np.nan})
T = pd.DataFrame(rows).round(3)
T.to_csv(ROOT / "results/v3/T_m2s_family_full.csv", index=False)
print(T.to_string())
