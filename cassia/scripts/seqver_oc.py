"""Operating characteristic of the sequential verification rules (Table 8b inset, Table S9).

    PYTHONPATH=. python scripts/seqver_oc.py   -> results/v3/T_seqver_oc.csv   (~10 min)

Both endpoints are held at a fixed ratio to their limits (S1: 0.5 and 2.0 mg/L) in all four rounds;
20,000 simulated sequences per ratio, rule and bias model (seed 1).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_v3 as r3  # noqa: E402

rng = np.random.default_rng(1)
rows = []
for ratio in [0.7, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.3, 1.5]:
    truths = [(ratio * 0.5, ratio * 2.0)] * 4
    for rule, (n, kind) in r3.RULES_SEQ.items():
        for bm, (rho, bc) in {"bias_shared": (0.5, 1.0), "bias_independent": (0.0, 0.0)}.items():
            N = 20000
            acc = sum(r3._sequence(rng, truths, (0.5, 2.0), n, kind, rho, bc)[0] is not None for _ in range(N))
            rows.append({"true_over_limit": ratio, "rule": rule, "bias_model": bm, "p_accept_within_4_rounds": acc / N})
T = pd.DataFrame(rows)
T.to_csv(r3.RES / "T_seqver_oc.csv", index=False)
print(T.pivot_table(index="true_over_limit", columns=["rule", "bias_model"], values="p_accept_within_4_rounds").round(4))
