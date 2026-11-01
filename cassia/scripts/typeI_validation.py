"""Type-I error of the adsorption-only rejection test on an independent cohort of negative controls.

    PYTHONPATH=. CTRL_PART=1 python scripts/typeI_validation.py   # controls, corpora 0-5   (~45 min, 4 cores)
    PYTHONPATH=. CTRL_PART=2 python scripts/typeI_validation.py   # controls, corpora 6-11  (~45 min)
    PYTHONPATH=. python scripts/typeI_validation.py               # summary -> results/v3/T_typeI_validation.csv

Calibration cohort: the 120 paired controls of run_v3.py stage equiv3. Validation cohort: 240 new controls
(12 new corpora x 20 studies; adsorption-only versions of scenario-E materials), fitted with exactly the same
procedure (equivalence_test_v3, 3 starts, warm start across designs). Thresholds compared:
  empirical : 95th percentile of calibration chi2 (method of the previous version)
  conformal : the ceil(0.95 (n + 1))-th smallest calibration chi2 (finite-sample type-I <= 0.05 for an
              exchangeable null study; Vovk et al., 2005)
Type-I errors on the validation cohort are reported with Wilson 95% intervals. Power is recomputed at the
conformal threshold on the study arm of stage equiv3.
"""
import os
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_v3 as r3  # noqa: E402
from cassia import identifiability as idf, simulate as sim  # noqa: E402

RES = r3.RES
SEED0 = 900000
N_CORPORA, N_ST = 12, 20


def work(args):
    rep, i0, i1 = args
    seed = SEED0 + rep
    _, mats = sim.virtual_corpus(N_ST, seed=seed, scenario="E")
    out = []
    for k, (_, meta) in enumerate(mats.iterrows()):
        if not (i0 <= k < i1):
            continue
        meta = meta.to_dict()
        mm = sim.apply_scenario(r3.v2.material_from(meta), "A")
        for sd in (0.03, 0.10):
            x_prev = None
            for d in idf.DESIGNS_V3:
                try:
                    r = idf.equivalence_test_v3(mm, meta, d, sd_rel=sd, seed=seed + k, n_starts=3, x_init=x_prev)
                except (RuntimeError, ValueError):
                    continue
                x_prev = r.pop("_x")
                out.append({"rep": rep, "study_id": meta["study_id"], "arm": "control_validation", **r})
    return out


def wilson(k, n, a=0.05):
    z = stats.norm.ppf(1 - a / 2); p = k / n; c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return max(0.0, c - h), min(1.0, c + h)


def conformal_threshold(x, alpha=0.05):
    x = np.sort(np.asarray(x)); k = int(np.ceil((1 - alpha) * (len(x) + 1)))
    return float(x[k - 1]) if k <= len(x) else np.inf


def summary():
    V = pd.concat([pd.read_csv(f) for f in sorted(RES.glob("typeI_ctrl_part*.csv.gz"))], ignore_index=True)
    Q = pd.read_csv(RES / "equiv3.csv.gz")
    rows = []
    for (sd, d), g in Q.groupby(["sd_rel", "design"]):
        cal = g[g.arm == "control_A"].chi2.values
        val_ = V[(V.sd_rel == sd) & (V.design == d)].chi2.values
        tru = g[g.arm == "true"]
        dof = int(g.dof.iloc[0])
        rec = {"sd_rel": sd, "design": d, "dof": dof, "n_calibration": len(cal), "n_validation": len(val_)}
        for name, thr in (("nominal", stats.chi2.ppf(0.95, dof)), ("empirical", np.quantile(cal, 0.95)),
                          ("conformal", conformal_threshold(cal))):
            k = int((val_ > thr).sum()); lo, hi = wilson(k, len(val_))
            rec |= {f"thr_{name}": thr, f"typeI_{name}": k / len(val_), f"typeI_{name}_lo": lo, f"typeI_{name}_hi": hi}
        # final threshold for reporting power: conformal on the pooled 360 controls
        thr_all = conformal_threshold(np.concatenate([cal, val_]))
        rec["thr_conformal_pooled"] = thr_all
        for lo_, hi_, lab in ((0.05, 0.3, "5-30%"), (0.3, 1.01, ">30%")):
            gg = tru[(tru.share_red > lo_) & (tru.share_red <= hi_)]
            rec[f"n|{lab}"] = len(gg)
            rec[f"power_conformal|{lab}"] = float((gg.chi2 > rec["thr_conformal"]).mean())
            rec[f"power_conformal_pooled|{lab}"] = float((gg.chi2 > thr_all).mean())
            rej = gg.chi2 > thr_all
            for eps in (0.05, 0.1):
                eq = gg[f"equiv_{eps}"].astype(bool) & ~rej
                rec[f"rejected|{lab}|eps={eps}"] = float(rej.mean())
                rec[f"equivalent|{lab}|eps={eps}"] = float(eq.mean())
                rec[f"inconclusive|{lab}|eps={eps}"] = float((~rej & ~eq).mean())
        rows.append(rec)
    T = pd.DataFrame(rows)
    T.to_csv(RES / "T_typeI_validation.csv", index=False)
    pd.set_option("display.width", 250)
    print(T[[c for c in T.columns if c.startswith(("sd", "design", "n_", "typeI", "power_conformal|"))]].round(3).to_string())


if __name__ == "__main__":
    part = int(os.environ.get("CTRL_PART", "0"))
    if part:
        reps = range(0, 6) if part == 1 else range(6, 12)
        jobs = [(r, i, i + 5) for r in reps for i in range(0, N_ST, 5)]
        with Pool(r3.NPROC) as p:
            V = pd.DataFrame([x for c in p.imap_unordered(work, jobs, chunksize=1) for x in c])
        V.to_csv(RES / f"typeI_ctrl_part{part}.csv.gz", index=False)
        print(f"typeI part {part} written", len(V))
    else:
        summary()
