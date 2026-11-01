"""Out-of-study accuracy of the four endpoint models used by the design layer (RF, speciation features).

    PYTHONPATH=. python scripts/endpoint_validation.py  -> results/v3/endpoint_validation.csv, T_endpoint_validation.csv
Same 10 corpora as the conformal analysis (run_v2.py stage uq). Total-Cr models are trained only on
rows with a total-Cr measurement (Class A/B studies).
"""
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_v2 as r2  # noqa: E402
from cassia import validation as val  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TARGETS = ("q_app_VI", "q_tot_aq", "Ce_VI", "Ce_total")


def work(k):
    d, _ = r2.corpus_features(r2.SEEDS["ml"] + k, 30)
    rows = []
    for tgt in TARGETS:
        for scheme in ("random", "LOSO"):
            res, pred_log, dd, _ = val.evaluate(d, "RF", val.CHEM, scheme, target=tgt)
            ylog = val.to_target(dd[tgt].values, tgt)
            e = pred_log - ylog
            r2log = 1 - np.sum(e ** 2) / np.sum((ylog - ylog.mean()) ** 2)
            rows.append({"corpus": k, "target": tgt, "scheme": scheme, "n_obs": res["n_obs"], "n_studies": res["n_studies"],
                         "R2_linear": res["R2"], "R2_log10": r2log, "RMSE_log10": float(np.sqrt(np.mean(e ** 2))),
                         "MedAE_log10": float(np.median(np.abs(e))), "bias_log10": float(np.mean(e))})
    return rows


if __name__ == "__main__":
    with Pool(4) as p:
        R = pd.DataFrame([x for rows in p.map(work, range(10)) for x in rows])
    R.to_csv(ROOT / "results/v3/endpoint_validation.csv", index=False)
    S = R.groupby(["target", "scheme"]).agg(n_studies=("n_studies", "mean"), R2_log10=("R2_log10", "mean"), R2_log10_sd=("R2_log10", "std"),
                                            R2_linear=("R2_linear", "mean"), RMSE_log10=("RMSE_log10", "mean"),
                                            MedAE_log10=("MedAE_log10", "mean"), bias_log10=("bias_log10", "mean")).round(3)
    S.to_csv(ROOT / "results/v3/T_endpoint_validation.csv")
    print(S.to_string())
