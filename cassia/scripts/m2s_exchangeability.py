"""M2s evaluated by the random role assignment of Appendix A (main M2s evaluation: Table 7b, Table S5c, Fig. 6).

    PYTHONPATH=. python scripts/m2s_exchangeability.py   (~15 min, 4 cores)
      -> results/v3/m2s_roles_per_study.csv.gz   one row per corpus x proper-training set x test study x alpha
         results/v3/T_m2s_roles_summary.csv      Table 7b (role-averaged), with the theoretical value k/(J+1)
         results/v3/T_m2s_exchangeability.csv    Table S5c (family-conditional and new-corpus coverage)
         results/v3/m2s_exchangeability.json

Design facts (cassia.simulate.virtual_corpus): 30 studies per corpus; material families are assigned in a fixed
rotation (study i has family i mod 6, 5 studies per family); family hyper-priors are drawn once per corpus and
shared by its studies; everything else is drawn independently per study. Studies are therefore NOT i.i.d. draws
from one population, and families are not randomised.

The M2s guarantee rests instead on the random assignment of roles: in the leave-one-study-out evaluation the test
study is each study in turn (equal weight) and the remaining studies are split at random into proper-training and
calibration studies. Over the test index and the split, (test, proper-training, calibration) is a uniformly random
partition of the corpus's studies, so given the corpus and the proper-training set, the test study is a uniformly
random member of the J + 1 non-training studies. This is the exchangeability the theorem needs, and it holds
conditionally on the realised corpus, whatever the family allocation and shared random effects.

This script IS the role-randomised evaluation used for the reported M2s results (same 10 corpora, model and
beta = 0.9 as run_v3.py stage m2s, whose single fixed-seed leave-one-study-out run is kept as one realisation):
  * R random proper-training sets per corpus and target, each drawn uniformly and independently (one generator per
    corpus, seeded 5000 + corpus); for each, every one of the J + 1 non-training studies is taken as test in turn,
    calibrated on the other J. Over the R draws and the enumeration this is a Monte Carlo average over uniformly
    random partitions (test, proper-training, calibration), i.e. over the role assignment of Appendix A;
  * within-corpus rate of studies with >= 90% of observations covered at alpha = 0.10 (marginal: guaranteed >= 0.9
    in expectation over role assignment, for every corpus);
  * the same rate conditional on the test study's family (not guaranteed);
  * transfer to a new corpus: model and calibration from corpus k, test studies from corpus k+1 (different
    hyper-prior draw; not covered by the guarantee).
"""
import json
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_v2 as v2  # noqa: E402
from cassia import validation as val  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/v3"
R, BETA, ALPHA, NCORP = 50, 0.9, 0.10, 10
ALPHAS = (0.5, 0.3, 0.2, 0.1, 0.05)
TARGETS = ("q_app_VI", "Ce_VI", "q_tot_aq", "Ce_total")


def study_score(r, beta=BETA):
    return np.sort(r)[int(np.ceil(beta * len(r))) - 1]


def work(k):
    d, _ = v2.corpus_features(v2.SEEDS["ml"] + k, 30)
    d2, _ = v2.corpus_features(v2.SEEDS["ml"] + (k + 1) % NCORP, 30)
    rng = np.random.default_rng(5000 + k)
    out = []
    for tgt in TARGETS:
        dd, ee = val.prepare(d, tgt), val.prepare(d2, tgt)
        y, ye = val.to_target(dd[tgt].values, tgt), val.to_target(ee[tgt].values, tgt)
        g, ge = dd.study_id.values, ee.study_id.values
        fam = dd.groupby("study_id").material_family.first()
        studies = np.unique(g); N = len(studies)
        J = int(round(0.5 * (N - 1))); P = N - 1 - J          # as in StudySplitConformal within LOSO
        kq = int(np.ceil((1 - ALPHA) * (J + 1)))
        for r in range(R):
            prop = rng.choice(studies, P, replace=False)
            pool = np.setdiff1d(studies, prop)
            m = val.make_model("RF", k).fit(dd[val.CHEM][np.isin(g, prop)], y[np.isin(g, prop)])
            res = np.abs(y - m.predict(dd[val.CHEM]))
            sc = {s: study_score(res[g == s]) for s in pool}
            for t in pool:                                         # each non-training study as test in turn
                others = np.sort([sc[s] for s in pool if s != t])
                rt = res[g == t]
                for a in ALPHAS:
                    ka = int(np.ceil((1 - a) * (J + 1)))
                    q = others[ka - 1] if ka <= J else np.inf
                    cov = float((rt <= q).mean())
                    out.append({"corpus": k, "target": tgt, "rep": r, "kind": "within", "study_id": t, "family": fam[t],
                                "alpha": a, "N": N, "J": J, "coverage": cov, "ok": cov >= BETA - 1e-9,
                                "unbounded": bool(np.isinf(q)), "width_log10": float(min(2 * q, 10)), "n": int(rt.size),
                                "theory": min(ka, J + 1) / (J + 1)})
            # transfer: calibration on a random J of the pool, tests from the next corpus
            cal = rng.choice(pool, J, replace=False)
            q = np.sort([sc[s] for s in cal])[kq - 1] if kq <= J else np.inf
            rese = np.abs(ye - m.predict(ee[val.CHEM]))
            for t in np.unique(ge):
                cov = float((rese[ge == t] <= q).mean())
                out.append({"corpus": k, "target": tgt, "rep": r, "kind": "new_corpus", "study_id": t, "alpha": ALPHA,
                            "family": ee[ge == t].material_family.iloc[0], "N": N, "J": J, "coverage": cov,
                            "ok": cov >= BETA - 1e-9, "unbounded": bool(np.isinf(q))})
    return out


def cluster_ci(per_corpus, B=4000, seed=11):
    v = np.asarray(per_corpus, float); rng = np.random.default_rng(seed)
    bs = v[rng.integers(0, len(v), (B, len(v)))].mean(1)
    return float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))


def main():
    with Pool(4) as p:
        U = pd.DataFrame([x for c in p.map(work, range(NCORP)) for x in c])
    U.to_csv(RES / "m2s_roles_per_study.csv.gz", index=False)
    W = U[U.kind == "within"]
    # Table 7b (role-averaged): rates over corpora x proper-training sets x enumerated test roles
    rows = []
    for (tgt, a), g in W.groupby(["target", "alpha"]):
        b = g[~g.unbounded]
        pc = g.groupby("corpus").ok.mean()
        rows.append({"target": tgt, "nominal_1-alpha": 1 - a, "n_role_assignments": len(g), "J_mean": g.groupby("corpus").J.first().mean(),
                     "frac_unbounded": g.unbounded.mean(), "frac_studies_cov_ge_0.9": g.ok.mean(), "theory_k_over_J+1": g.theory.mean(),
                     "min_corpus": pc.min(), "mean_cov": g.coverage.mean(),
                     "mean_width_log10_bounded": b.width_log10.mean() if len(b) else np.nan})
    S = pd.DataFrame(rows)
    S.to_csv(RES / "T_m2s_roles_summary.csv", index=False)
    # Table S5c at alpha = 0.10: family-conditional and new-corpus coverage (not covered by the theorem)
    rows = []
    U10 = U[np.isclose(U.alpha, ALPHA)]
    for tgt, g in U10.groupby("target"):
        w, n = g[g.kind == "within"], g[g.kind == "new_corpus"]
        pc = w.groupby("corpus").ok.mean()
        rec = {"target": tgt, "J_mean": w.groupby("corpus").J.first().mean(), "frac_unbounded": w.unbounded.mean(),
               "within_rate": w.ok.mean(), "within_min_corpus": pc.min(), "within_max_corpus": pc.max(),
               "new_corpus_rate": n.ok.mean(), "new_corpus_min_corpus": n.groupby("corpus").ok.mean().min()}
        for f, h in w.groupby("family"):
            lo, hi = cluster_ci(h.groupby("corpus").ok.mean().values)
            rec[f"family_{f}"] = h.ok.mean(); rec[f"family_{f}_lo"] = lo; rec[f"family_{f}_hi"] = hi
        fr = w.groupby("family").ok.mean()
        rec["family_min"], rec["family_max"] = fr.min(), fr.max()
        rows.append(rec)
    T = pd.DataFrame(rows)
    T.to_csv(RES / "T_m2s_exchangeability.csv", index=False)
    (RES / "m2s_exchangeability.json").write_text(json.dumps({"R": R, "alphas": ALPHAS, "beta": BETA, "corpora": NCORP,
                                                              "table7b": S.round(4).to_dict("records"),
                                                              "tableS5c": T.round(4).to_dict("records")}, indent=1))
    pd.set_option("display.width", 250)
    print(S.round(3).to_string()); print(T.round(3).T.to_string())


if __name__ == "__main__":
    main()
