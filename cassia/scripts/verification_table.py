"""Consolidated accepted-design counts and unsafe rates for both verification rules (Table 8c).

    python scripts/verification_table.py  -> results/v3/T_verification_denominators.csv

Two sources, same 300 new materials per target (20 cohorts x 15):
  * design stage (run_v2.py stage design, rule nominal_verify = simple rule): one realised verification
    sequence per material -> observed counts, Wilson 95% upper bounds with two denominators
    (all materials; accepted designs only).
  * sequential-verification stage (run_v3.py stage seqver): 1,000 simulated sequences per material and rule
    -> expected numbers of accepted designs and of unsafe acceptances, and the conditional unsafe rate.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]


def wilson_hi(k, n, a=0.05):
    z = norm.ppf(1 - a / 2); p = k / n
    return (p + z * z / (2 * n) + z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / (1 + z * z / n)


O = pd.read_csv(ROOT / "results/v2/design_outcomes.csv")
O = O[O.rule == "nominal_verify"]
rows = []
for t, g in O.groupby("target"):
    acc = g.feasible.astype(bool); uns = acc & ~g.met_both.astype(bool)
    n, a, u = len(g), int(acc.sum()), int(uns.sum())
    rows.append({"source": "design stage (1 realisation)", "target": t, "rule": "simple (n = 2, mean <= 0.9 limit)",
                 "bias_model": "shared within study", "materials": n, "accepted": a, "unsafe_accepted": u,
                 "unsafe_rate_all_materials": u / n, "wilson95_upper_all": wilson_hi(u, n),
                 "unsafe_rate_among_accepted": u / a, "wilson95_upper_among_accepted": wilson_hi(u, a)})
R = pd.read_csv(ROOT / "results/v3/seqver.csv")
lab = {"v2_ratio_n2": "simple (n = 2, mean <= 0.9 limit)", "ucb_n3_bonf": "error-controlled (n = 3, UCB <= limit)"}
for (t, rule, bm), g in R.groupby(["target", "rule", "bias_model"]):
    ea, eu = g.p_accept.sum(), g.p_false_accept.sum()
    rows.append({"source": "seqver (1,000 sequences per material)", "target": t, "rule": lab[rule],
                 "bias_model": bm.replace("bias_", ""), "materials": len(g), "accepted": ea, "unsafe_accepted": eu,
                 "unsafe_rate_all_materials": eu / len(g), "wilson95_upper_all": np.nan,
                 "unsafe_rate_among_accepted": eu / ea, "wilson95_upper_among_accepted": np.nan,
                 "max_p_false_accept_material": g.p_false_accept.max(), "mean_batches": g.mean_batches.mean()})
T = pd.DataFrame(rows)
T.to_csv(ROOT / "results/v3/T_verification_denominators.csv", index=False)
pd.set_option("display.width", 250)
print(T.round(4).to_string())
