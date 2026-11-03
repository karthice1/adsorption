"""Sensitivity of the A1 rejection test to misspecified XPS/XANES error models (Table S16).

    PYTHONPATH=. python scripts/meas_sensitivity.py            (about 1-2 h on 4 cores; resumable per setting)
      -> results/v3/meas_sens_warmstarts.json   start values of the nominal D1->D4x chain (as in stage equiv3)
         results/v3/meas_sens_<setting>.csv.gz   one row per study, arm (true / control_A) and setting
         results/v3/T_meas_sensitivity.csv       Table S16
         results/v3/meas_sensitivity.json

Design. The 120 studies and paired adsorption-only controls of stage equiv3 (scenario E replicates 0-2, B/C/D
replicate 0; 3% analytical noise) are re-simulated for designs D4x (XPS) and D4n (XANES). The oxidation-state
data are generated with a perturbed technique model (true_tech), while the analyst's adsorption-only fit keeps
the nominal model of Table S1 (nominal SDs and beam-reduction bound): the true analytical error is larger than
assumed. For each setting we report
  * type-I error of the published study-pooled threshold (360 nominal controls, Table S14b), applied unchanged;
  * a threshold recalibrated on the 120 perturbed controls (95th percentile), its split-half cross-validated
    type-I error, and the power to reject adsorption-only for reduction shares 5-30% and > 30%.
The 'nominal' settings repeat the procedure with the nominal model (same code path) as the reference; with the
warm starts of stage equiv3 they reproduce its D4x/D4n fits.
"""
from __future__ import annotations

import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
import run_v2 as v2  # noqa: E402
from cassia import identifiability as idf, simulate as sim  # noqa: E402

RES = ROOT / "results" / "v3"
SEED0 = 500000                      # = run_v3.SEEDS3["equiv"]: same studies as stage equiv3
JOBS = {"E": [0, 1, 2], "B": [0], "C": [0], "D": [0]}
N_ST, SD_REL = 20, 0.03
SETTINGS = {                        # name: (design, overrides of cassia.measurement.TECHNIQUES[tech] for data generation)
    "XPS_nominal":           ("D4x", {}),
    "XPS_noise_x2":          ("D4x", {"noise": 0.10}),
    "XPS_enrichment_sd_x2":  ("D4x", {"enrich_logsd": 0.8}),
    "XPS_beam_x2":           ("D4x", {"beam_max": 0.30}),
    "XANES_nominal":         ("D4n", {}),
    "XANES_beam_0.10":       ("D4n", {"beam_max": 0.10}),
    "XANES_noise_x2":        ("D4n", {"noise": 0.06}),
}


WARM = RES / "meas_sens_warmstarts.json"


def _studies(scen, rep, i0, i1):
    seed = SEED0 + 1000 * sim.SCENARIOS.index(scen) + rep
    _, mats = sim.virtual_corpus(N_ST, seed=seed, scenario=scen)
    for k, (_, meta) in enumerate(mats.iterrows()):
        if i0 <= k < i1:
            meta = meta.to_dict()
            m = v2.material_from(meta)
            for arm, mm in (("true", m), ("control_A", sim.apply_scenario(m, "A"))):
                yield seed + k, meta, arm, mm


def warm_worker(args):
    """Nominal warm-start chain D1 -> D2 -> D3 -> D4x exactly as in stage equiv3 (run_v3._equiv3_worker, 3% noise):
    the start values for D4x (fit of D3) and D4n (fit of D4x). D1-D3 contain no oxidation-state data, so their fits
    do not depend on the technique model; using the nominal D4x fit as the D4n start is an optimiser choice only."""
    scen, rep, i0, i1 = args
    out = {}
    for seed, meta, arm, mm in _studies(scen, rep, i0, i1):
        x_prev, inits = None, {}
        for d in ("D1", "D2", "D3", "D4x"):
            if d == "D4x":
                inits["D4x"] = None if x_prev is None else [float(v) for v in x_prev]
            try:
                r = idf.equivalence_test_v3(mm, meta, d, sd_rel=SD_REL, seed=seed, n_starts=3, x_init=x_prev)
            except (RuntimeError, ValueError):
                continue
            x_prev = r.pop("_x")
        inits["D4n"] = None if x_prev is None else [float(v) for v in x_prev]
        out[f"{scen}|{rep}|{meta['study_id']}|{arm}"] = inits
    return out


def worker(args):
    setting, scen, rep, i0, i1 = args
    design, over = SETTINGS[setting]
    W = json.loads(WARM.read_text())
    tech = idf.TECH_OF[design]
    out = []
    for seed, meta, arm, mm in _studies(scen, rep, i0, i1):
        x0 = W[f"{scen}|{rep}|{meta['study_id']}|{arm}"][design]
        try:
            r = idf.equivalence_test_v3(mm, meta, design, sd_rel=SD_REL, seed=seed, n_starts=3,
                                        x_init=None if x0 is None else np.array(x0),
                                        true_tech=dict(over) if over else {"noise": idf.ms.TECHNIQUES[tech]["noise"]})
        except (RuntimeError, ValueError):
            continue
        r.pop("_x")
        out.append({"setting": setting, "scenario": scen, "rep": rep, "study_id": meta["study_id"], "arm": arm, **r})
    return out


def summarise():
    T = pd.read_csv(RES / "T_typeI_validation.csv")
    pub = {d: float(T[(T.sd_rel == SD_REL) & (T.design == d)].thr_conformal_pooled.iloc[0]) for d in ("D4x", "D4n")}
    rng = np.random.default_rng(1)
    rows = []
    for s, (design, over) in SETTINGS.items():
        f = RES / f"meas_sens_{s}.csv.gz"
        if not f.exists():
            continue
        g = pd.read_csv(f)
        ctrl, tru = g[g.arm == "control_A"], g[g.arm == "true"]
        idx = rng.permutation(len(ctrl)); h = len(idx) // 2
        c1, c2 = ctrl.chi2.values[idx[:h]], ctrl.chi2.values[idx[h:]]
        cv = 0.5 * ((c2 > np.quantile(c1, 0.95)).mean() + (c1 > np.quantile(c2, 0.95)).mean())
        thr = float(np.quantile(ctrl.chi2, 0.95))
        rec = {"setting": s, "design": design, "perturbation": json.dumps(over) if over else "none",
               "n_controls": len(ctrl), "published_threshold": pub[design],
               "typeI_published_threshold": float((ctrl.chi2 > pub[design]).mean()),
               "recalibrated_threshold": thr, "typeI_recalibrated_cv": float(cv)}
        for lo, hi, lab in ((0.05, 0.3, "5-30%"), (0.3, 1.01, ">30%")):
            gg = tru[(tru.share_red > lo) & (tru.share_red <= hi)]
            rec[f"n|{lab}"] = len(gg)
            rec[f"power_published|{lab}"] = float((gg.chi2 > pub[design]).mean())
            rec[f"power_recalibrated|{lab}"] = float((gg.chi2 > thr).mean())
        rows.append(rec)
    S = pd.DataFrame(rows)
    S.to_csv(RES / "T_meas_sensitivity.csv", index=False)
    (RES / "meas_sensitivity.json").write_text(json.dumps({"settings": {k: [d, o] for k, (d, o) in SETTINGS.items()},
                                                           "sd_rel": SD_REL, "table": S.round(4).to_dict("records")}, indent=1))
    pd.set_option("display.width", 250)
    print(S.round(3).to_string())


def main():
    if not WARM.exists():
        jobs = [(sc, r, i, i + 5) for sc, rs in JOBS.items() for r in rs for i in range(0, N_ST, 5)]
        W = {}
        with Pool(4) as p:
            for c in p.imap_unordered(warm_worker, jobs, chunksize=1):
                W.update(c)
        WARM.write_text(json.dumps(dict(sorted(W.items()))))
        print("written", WARM.name, len(W), flush=True)
    for s in SETTINGS:
        f = RES / f"meas_sens_{s}.csv.gz"
        if f.exists():
            continue
        jobs = [(s, sc, r, i, i + 5) for sc, rs in JOBS.items() for r in rs for i in range(0, N_ST, 5)]
        with Pool(4) as p:
            D = pd.DataFrame([x for c in p.imap_unordered(worker, jobs, chunksize=1) for x in c])
        D.to_csv(f, index=False)
        print("written", f.name, len(D), flush=True)
    summarise()


if __name__ == "__main__":
    summarise() if "--summarise" in sys.argv else main()
