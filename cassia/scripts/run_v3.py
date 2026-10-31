"""CASSIA v0.3 experiment driver (second revision). Adds to run_v2.py:

  chem3    background-electrolyte sensitivity; undersaturation thresholds by database and I
  ident3   technique-aware Fisher analysis (XPS vs XANES vs inventory/phase) and detection tests
  equiv3   calibrated goodness-of-fit rejection + TOST-type equivalence test with margins,
           operating characteristics against the noise-free truth
  m2s      exact study-level split conformal (M2s) for the four endpoints (10 corpora)
  mlscen   ML validation and intervals for scenarios A-E (5 corpora each) + cross-mechanism transfer
  censor   sensitivity of ML results to censoring / LOD treatment
  lc3      learning curve, 8 replicates per corpus size
  seqver   error-controlled sequential verification of treatment designs
  figures  updated figures 3-6 and new figure 8

    python scripts/run_v3.py --stage all       # ~2-3 h on 4 cores
Outputs: results/v3/. Seeds: SEEDS3 below and run_v2.SEEDS.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import json
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_v2 as v2  # noqa: E402
from cassia import chemistry as ch, identifiability as idf, simulate as sim, validation as val  # noqa: E402

ROOT = v2.ROOT
RES = ROOT / "results" / "v3"
RES.mkdir(parents=True, exist_ok=True)
FIG = v2.FIG
NPROC = 4
SEEDS3 = {"ident": 400000, "equiv": 500000, "mlscen": 600000, "lc": 700000, "seq": 800000}
C, INK, MUTED = v2.C, v2.INK, v2.MUTED
panel, wilson = v2.panel, v2.wilson


def dump(name, obj):
    (RES / f"{name}.json").write_text(json.dumps(v2._clean(obj), indent=2, default=str))


# ============================================================================= chem3
def stage_chem3(quick=False):
    out = {}
    rows = []
    for e in ("NaCl", "NaNO3", "Na2SO4"):
        for I in (0.01, 0.1):
            for pH in (3.0, 5.0, 6.5):
                fr = ch.speciate_cr6(pH, 10.0, ionic_strength=I, electrolyte=e).fractions
                rows.append({"electrolyte": e, "I": I, "pH": pH, "f_HCrO4": fr.get("HCrO4-", 0), "f_CrO4": fr.get("CrO4-2", 0),
                             "f_Cr2O7": fr.get("Cr2O7-2", 0), "f_CrO3SO4": fr.get("CrO3SO4-2", 0)})
            p = np.round(np.arange(5.5, 7.21, 0.02), 3)
            d = []
            for x in p:
                fr = ch.speciate_cr6(float(x), 10.0, ionic_strength=I, electrolyte=e).fractions
                d.append(fr.get("HCrO4-", 0) - fr.get("CrO4-2", 0))
            out[f"crossover|{e}|I={I}"] = float(np.interp(0, -np.array(d), p))
    T = pd.DataFrame(rows)
    T.to_csv(RES / "T_electrolyte_sensitivity.csv", index=False)
    out["speciation"] = T.round(4).to_dict("records")
    S = pd.read_csv(v2.RES / "T_cr3_solubility_SI.csv")
    th = S[S.status == "undersaturated_to_cmax"].groupby(["database", "I"]).pH.max()
    out["undersat_max_pH_by_db_I"] = {f"{a}|I={b}": float(v) for (a, b), v in th.items()}
    dump("chem3", out)
    print("chem3 done")


# ============================================================================= ident3
def _ident3_worker(args):
    scen, rep, n_st = args
    seed = SEEDS3["ident"] + 1000 * sim.SCENARIOS.index(scen) + rep
    _, mats = sim.virtual_corpus(n_st, seed=seed, scenario=scen)
    rng = np.random.default_rng(seed)
    out = []
    for _, meta in mats.iterrows():
        meta = meta.to_dict()
        m = v2.material_from(meta)
        try:
            fs = idf.fisher_study_v3(m, meta)
        except RuntimeError:
            continue
        rec = {"scenario": scen, "rep": rep, "study_id": meta["study_id"], "share_red": fs["base"]["share_red"],
               "share_precip": fs["base"]["share_precip"]}
        for d in idf.DESIGNS_V3:
            r = idf.posterior_sd_v3(fs, d)
            rec.update({f"{k}_{d}": v for k, v in r.items() if k != "design"})
        rec.update(idf.detection_tests_v3(rng, fs, meta, n_rep=20))
        mA = sim.apply_scenario(m, "A")
        fsA = idf.fisher_study_v3(mA, meta, designs=("D1",))
        rec.update({f"ctrl_{k}": v for k, v in idf.detection_tests_v3(rng, fsA, meta, n_rep=20).items()})
        out.append(rec)
    return out


def stage_ident3(quick=False):
    n_st = 6 if quick else 25
    reps = {"E": 2 if quick else 8, "B": 1 if quick else 2, "C": 1 if quick else 2, "D": 1 if quick else 2}
    with Pool(NPROC) as pool:
        F = pd.DataFrame([x for c in pool.imap_unordered(_ident3_worker, [(s, r, n_st) for s, n in reps.items() for r in range(n)]) for x in c])
    F.to_csv(RES / "ident3_fim_detection.csv", index=False)
    bins = [0, 0.05, 0.2, 0.4, 1.0]
    F["red_bin"] = pd.cut(F.share_red, bins, include_lowest=True)
    summ = {"n": F.groupby("scenario").size().to_dict(), "red_bins": [str(b) for b in F.red_bin.cat.categories]}
    for d in idf.DESIGNS_V3:
        summ[f"post_sd_red_by_bin|{d}"] = F.groupby("red_bin", observed=True)[f"post_sd_red_{d}"].median().round(4).tolist()
        pr = F[F.share_precip > 0.02]
        summ[f"post_sd_precip_median_precip>0.02|{d}"] = float(pr[f"post_sd_precip_{d}"].median()) if len(pr) else None
    summ["n_precip>0.02"] = int((F.share_precip > 0.02).sum())
    for k in ("dissolved_cr3", "xps_cr3", "xanes_cr3"):
        summ[f"detect_by_bin|{k}"] = F.groupby("red_bin", observed=True)[k].mean().round(3).tolist()
        summ[f"false_alarm_controls|{k}"] = float(F[f"ctrl_{k}"].mean())
    dump("ident3", summ)
    print("ident3 done", len(F))


# ============================================================================= equiv3
def _equiv3_worker(args):
    scen, rep, n_st = args
    seed = SEEDS3["equiv"] + 1000 * sim.SCENARIOS.index(scen) + rep
    _, mats = sim.virtual_corpus(n_st, seed=seed, scenario=scen)
    out = []
    for k, (_, meta) in enumerate(mats.iterrows()):
        meta = meta.to_dict()
        m = v2.material_from(meta)
        for arm, mm in (("true", m), ("control_A", sim.apply_scenario(m, "A"))):
            for sd in (0.03, 0.10):
                for d in idf.DESIGNS_V3:
                    try:
                        r = idf.equivalence_test_v3(mm, meta, d, sd_rel=sd, seed=seed + k, n_starts=2)
                    except (RuntimeError, ValueError):
                        continue
                    out.append({"scenario": scen, "rep": rep, "study_id": meta["study_id"], "arm": arm, **r})
    return out


def stage_equiv3(quick=False):
    n_st = 3 if quick else 20
    reps = {"E": 1 if quick else 3, "B": 0 if quick else 1, "C": 0 if quick else 1, "D": 0 if quick else 1}
    with Pool(NPROC) as pool:
        Q = pd.DataFrame([x for c in pool.imap_unordered(_equiv3_worker, [(s, r, n_st) for s, n in reps.items() for r in range(n)]) for x in c])
    Q.to_csv(RES / "equiv3.csv.gz", index=False)
    rng = np.random.default_rng(1)
    summ = {"n_studies": int(Q[Q.arm == "true"].groupby(["scenario", "rep", "study_id"]).ngroups)}
    rows = []
    for (sd, d), g in Q.groupby(["sd_rel", "design"]):
        ctrl = g[g.arm == "control_A"]
        dof = int(g.dof.iloc[0])
        # calibrated rejection threshold: 95th percentile of control chi2 (split-half cross-validated type-I)
        idx = rng.permutation(len(ctrl)); h = len(idx) // 2
        c1, c2 = ctrl.chi2.values[idx[:h]], ctrl.chi2.values[idx[h:]]
        t1, t2 = np.quantile(c1, 0.95), np.quantile(c2, 0.95)
        typeI_cv = 0.5 * ((c2 > t1).mean() + (c1 > t2).mean())
        thr = float(np.quantile(ctrl.chi2, 0.95))
        nominal = float(stats.chi2.ppf(0.95, dof))
        tru = g[g.arm == "true"].copy()
        rec = {"sd_rel": sd, "design": d, "dof": dof, "n_controls": len(ctrl), "nominal_threshold": nominal,
               "calibrated_threshold": thr, "typeI_nominal": float((ctrl.chi2 > nominal).mean()),
               "typeI_calibrated_cv": float(typeI_cv)}
        for lo_, hi_, lab in ((0.05, 0.3, "5-30%"), (0.3, 1.01, ">30%")):
            gg = tru[(tru.share_red > lo_) & (tru.share_red <= hi_)]
            rec[f"n|{lab}"] = len(gg)
            rec[f"power_calibrated|{lab}"] = float((gg.chi2 > thr).mean()) if len(gg) else None
            for eps in (0.05, 0.1, 0.2):
                rej = gg.chi2 > thr
                eq = gg[f"equiv_{eps}"].astype(bool) & ~rej
                rec[f"rejected|{lab}|eps={eps}"] = float(rej.mean()) if len(gg) else None
                rec[f"equivalent|{lab}|eps={eps}"] = float(eq.mean()) if len(gg) else None
                rec[f"inconclusive|{lab}|eps={eps}"] = float((~rej & ~eq).mean()) if len(gg) else None
        for eps in (0.05, 0.1, 0.2):
            allr = pd.concat([tru, ctrl])
            non = allr[~allr[f"truly_within_{eps}"].astype(bool)]
            win = allr[allr[f"truly_within_{eps}"].astype(bool)]
            rec[f"false_equivalence|eps={eps}"] = float(non[f"equiv_{eps}"].astype(bool).mean()) if len(non) else None
            rec[f"n_truly_outside|eps={eps}"] = len(non)
            rec[f"equivalence_power|eps={eps}"] = float(win[f"equiv_{eps}"].astype(bool).mean()) if len(win) else None
            rec[f"controls_equivalent|eps={eps}"] = float(ctrl[f"equiv_{eps}"].astype(bool).mean())
            rec[f"high_red_truly_within|eps={eps}"] = float(tru[tru.share_red > 0.3][f"truly_within_{eps}"].astype(bool).mean())
        rows.append(rec)
    S = pd.DataFrame(rows)
    S.to_csv(RES / "T_equiv3_summary.csv", index=False)
    summ["summary"] = S.round(4).to_dict("records")
    dump("equiv3", summ)
    print("equiv3 done", len(Q))


# ============================================================================= m2s
def _m2s_worker(args):
    k, quick = args
    d, _ = v2.corpus_features(v2.SEEDS["ml"] + k, 12 if quick else 30)
    rows = []
    for tgt in ("q_app_VI", "q_tot_aq", "Ce_VI", "Ce_total"):
        dd = val.prepare(d, tgt)
        y = val.to_target(dd[tgt].values, tgt)
        g = dd["study_id"].values
        for test in np.unique(g):
            tr, te = g != test, g == test
            sp = val.StudySplitConformal("RF", val.CHEM, 0.5, 0.9, seed=k).fit(dd[tr], tgt)
            for a in (0.5, 0.3, 0.2, 0.1, 0.05):
                lo, hi, _ = sp.interval(dd[te], a)
                cov = (y[te] >= lo) & (y[te] <= hi)
                rows.append({"corpus": k, "target": tgt, "study_id": test, "family": dd["material_family"].values[te][0],
                             "alpha": a, "coverage": float(cov.mean()), "n": int(te.sum()),
                             "width_log10": float(np.minimum(hi - lo, 10).mean()), "n_cal": sp.J,
                             "unbounded": bool(np.isinf(hi).any())})
    return rows


def stage_m2s(quick=False):
    with Pool(NPROC) as pool:
        U = pd.DataFrame([x for c in pool.map(_m2s_worker, [(k, quick) for k in range(2 if quick else 10)]) for x in c])
    U.to_csv(RES / "m2s_per_study.csv", index=False)
    rows = []
    for (tgt, a), g in U.groupby(["target", "alpha"]):
        rows.append({"target": tgt, "nominal_1-alpha": 1 - a, "n_test_studies": len(g), "n_cal_mean": g.n_cal.mean(),
                     "frac_unbounded": g.unbounded.mean(), "frac_studies_cov_ge_0.9": float((g.coverage >= 0.9 - 1e-9).mean()),
                     "mean_cov": g.coverage.mean(), "worst_cov": g.coverage.min(),
                     "mean_width_log10_bounded": g[~g.unbounded].width_log10.mean() if (~g.unbounded).any() else np.nan})
    S = pd.DataFrame(rows)
    S.to_csv(RES / "T_m2s_summary.csv", index=False)
    fam = U[U.alpha == 0.1].groupby(["target", "family"]).apply(lambda x: (x.coverage >= 0.9 - 1e-9).mean(), include_groups=False).unstack().round(3)
    fam.to_csv(RES / "T_m2s_family.csv")
    dump("m2s", {"summary": S.round(4).to_dict("records"), "family_frac_ge_0.9": fam.to_dict()})
    print("m2s done")


# ============================================================================= mlscen
def _mlscen_worker(args):
    scen, k, quick = args
    seed = SEEDS3["mlscen"] + 1000 * sim.SCENARIOS.index(scen) + k
    n = 12 if quick else 30
    df, _ = sim.virtual_corpus(n, seed=seed, scenario=scen)
    d = v2.add_features(df)
    rows, uq = [], []
    for mname in ("RF", "XGB"):
        for fs_name, fs in (("base", val.BASE), ("chem", val.CHEM)):
            for sc in ("random", "LOSO", "LOFO"):
                r, *_ = val.evaluate(d, mname, fs, sc)
                rows.append(r | {"scenario": scen, "corpus": k})
    # cross-mechanism transfer: train on a scenario-E corpus (same seed index), test on this scenario
    if scen != "E":
        dfe, _ = sim.virtual_corpus(n, seed=SEEDS3["mlscen"] + 1000 * sim.SCENARIOS.index("E") + k, scenario="E")
        de = v2.add_features(dfe)
        for mname in ("RF", "XGB"):
            mdl = val.make_model(mname).fit(de[val.CHEM], val.to_target(de.q_app_VI.values))
            sm = val.duan_smear(val.to_target(de.q_app_VI.values) - mdl.predict(de[val.CHEM]))
            qp = val.from_target(mdl.predict(d[val.CHEM]), sm)
            rows.append({"model": mname, "features": "chem", "scheme": "transfer_from_E", "scenario": scen, "corpus": k,
                         **val.metrics(np.clip(d.q_app_VI.values, 0, None), qp)})
    cf = val.conformal_loso(d, "q_app_VI", "RF", val.CHEM, K=10, alphas=(0.1,), beta=0.9, seed=k)
    cf["scenario"], cf["corpus"] = scen, k
    return rows, cf


def stage_mlscen(quick=False):
    jobs = [(s, k, quick) for s in sim.SCENARIOS for k in range(1 if quick else 5)]
    with Pool(NPROC) as pool:
        res = pool.map(_mlscen_worker, jobs, chunksize=1)
    M = pd.DataFrame([r for rows, _ in res for r in rows])
    U = pd.concat([u for _, u in res], ignore_index=True)
    M.to_csv(RES / "mlscen_metrics.csv", index=False)
    U.to_csv(RES / "mlscen_uq.csv", index=False)
    A = M.groupby(["scenario", "model", "features", "scheme"]).R2.agg(["mean", "std", "count"]).round(4).reset_index()
    A.to_csv(RES / "T_mlscen_summary.csv", index=False)
    UQ = U.groupby(["scenario", "method"]).apply(lambda x: pd.Series({"mean_cov": x.coverage.mean(),
                                                                     "frac_ge_0.9": (x.coverage >= 0.9 - 1e-9).mean(),
                                                                     "width": x.width_log10.mean()}), include_groups=False).round(3).reset_index()
    UQ.to_csv(RES / "T_mlscen_uq.csv", index=False)
    piv = M[M.scheme == "LOSO"].pivot_table(index=["scenario", "corpus"], columns=["model", "features"], values="R2")
    chem_gain = {s: float((piv.loc[s][("RF", "chem")] - piv.loc[s][("RF", "base")]).mean()) for s in sim.SCENARIOS}
    dump("mlscen", {"summary": A.to_dict("records"), "uq": UQ.to_dict("records"), "RF_chem_minus_base_LOSO_by_scenario": chem_gain})
    print("mlscen done")


# ============================================================================= censor
def _censor_worker(args):
    k, quick = args
    df, _ = sim.virtual_corpus(12 if quick else 30, seed=v2.SEEDS["ml"] + k, scenario="E")
    out = []
    variants = {"log_censor_LODhalf": ("log", None), "asinh_keep_negative": ("asinh", None),
                "log_drop_censored": ("log", "drop"), "log_LOD_full": ("log", "lodfull")}
    for name, (mode, opt) in variants.items():
        dd = df.copy()
        if opt == "lodfull":
            dd.loc[dd.Ce_VI_below_LOD, "Ce_VI"] = 0.01
        d = v2.add_features(dd)
        val.TRANSFORM["mode"] = mode
        try:
            if opt == "drop":
                keep = d.q_app_VI >= 0
                res, pl, dsub, qp = val.evaluate(d[keep].reset_index(drop=True), "RF", val.CHEM, "LOSO")
            else:
                res, pl, dsub, qp = val.evaluate(d, "RF", val.CHEM, "LOSO")
            cf = val.conformal_loso(d if opt != "drop" else d[d.q_app_VI >= 0].reset_index(drop=True), "q_app_VI", "RF",
                                    val.CHEM, K=5, alphas=(0.1,), with_split=False)
            m1 = cf[cf.method == "M1"]
            out.append({"corpus": k, "variant": name, "R2": res["R2"], "RMSE": res["RMSE"], "bias": res["bias"],
                        "n_obs": res["n_obs"], "n_censored": int((d.q_app_VI < 0).sum()),
                        "M1_mean_cov": float(m1.coverage.mean()), "M1_frac_ge_0.9": float((m1.coverage >= 0.9 - 1e-9).mean())})
        finally:
            val.TRANSFORM["mode"] = "log"
    return out


def stage_censor(quick=False):
    with Pool(NPROC) as pool:
        R = pd.DataFrame([x for c in pool.map(_censor_worker, [(k, quick) for k in range(2 if quick else 5)]) for x in c])
    R.to_csv(RES / "censor_sensitivity.csv", index=False)
    S = R.groupby("variant")[["R2", "RMSE", "bias", "M1_mean_cov", "M1_frac_ge_0.9", "n_censored"]].mean().round(4)
    S.to_csv(RES / "T_censor_summary.csv")
    dump("censor", {"summary": S.reset_index().to_dict("records")})
    print("censor done")


# ============================================================================= lc3
def _lc3_worker(args):
    n, rep = args
    d, _ = v2.corpus_features(SEEDS3["lc"] + 100 * n + rep, n)
    r, *_ = val.evaluate(d, "RF", val.CHEM, "LOSO")
    r2, *_ = val.evaluate(d, "RF", val.CHEM, "random")
    return {"n_studies": n, "rep": rep, "R2_LOSO": r["R2"], "RMSE_LOSO": r["RMSE"], "R2_random": r2["R2"]}


def stage_lc3(quick=False):
    ns = (12, 24) if quick else (12, 24, 48, 96)
    with Pool(NPROC) as pool:
        L = pd.DataFrame(pool.map(_lc3_worker, [(n, r) for n in ns for r in range(2 if quick else 8)], chunksize=1))
    L.to_csv(RES / "lc3.csv", index=False)
    S = L.groupby("n_studies").agg(R2_LOSO_mean=("R2_LOSO", "mean"), R2_LOSO_sd=("R2_LOSO", "std"),
                                   R2_LOSO_se=("R2_LOSO", lambda x: x.std() / np.sqrt(len(x))),
                                   R2_random_mean=("R2_random", "mean"), n=("R2_LOSO", "size")).round(4)
    rho = stats.spearmanr(L.n_studies, L.R2_LOSO)
    S.to_csv(RES / "T_lc3.csv")
    dump("lc3", {"summary": S.reset_index().to_dict("records"), "spearman_n_vs_R2": [float(rho.statistic), float(rho.pvalue)]})
    print("lc3 done")


# ============================================================================= seqver
RULES_SEQ = {
    # name: (replicates per round, acceptance: 'ratio' (mean <= 0.9 limit) or 'ucb' (upper confidence bound <= limit))
    "v2_ratio_n2": (2, "ratio"),
    "ucb_n3_bonf": (3, "ucb"),
}
SIG_R, SIG_B = 0.03, 0.05      # analytical random SD (relative), calibration-bias SD (relative, shared)


def _sequence(rng, truths, lims, n_rep, rule, rho_e=0.5, bias_corr=1.0, alpha=0.05, rounds=4):
    """One verification sequence. truths: list of (Ce_VI, Ce_tot) along the dose ladder."""
    b6 = rng.normal(0, SIG_B)
    bt = bias_corr * b6 + np.sqrt(1 - bias_corr ** 2) * rng.normal(0, SIG_B)
    z = stats.norm.ppf(1 - alpha / (2 * rounds))           # Bonferroni over 2 endpoints x rounds
    for r, (c6, ct) in enumerate(truths[:rounds], start=1):
        e = rng.multivariate_normal([0, 0], [[SIG_R ** 2, rho_e * SIG_R ** 2], [rho_e * SIG_R ** 2, SIG_R ** 2]], size=n_rep)
        m6 = np.mean(c6 * np.exp(b6 + e[:, 0])); mt = np.mean(ct * np.exp(bt + e[:, 1]))
        if rule == "ratio":
            ok = (m6 <= 0.9 * lims[0]) and (mt <= 0.9 * lims[1])
        else:
            f = np.exp(z * np.sqrt(SIG_R ** 2 / n_rep + SIG_B ** 2))   # one-sided bound on the log scale
            ok = (m6 * f <= lims[0]) and (mt * f <= lims[1])
        if ok:
            return r, (c6 > lims[0]) or (ct > lims[1]), r * n_rep
    return None, False, min(len(truths), rounds) * n_rep


def _seqver_worker(args):
    k, quick, n_new, n_mc = args
    seed = v2.SEEDS["ml"] + k
    d, _ = v2.corpus_features(seed, 12 if quick else 30)
    cps = {t: val.StudyCVPlus("RF", val.CHEM, K=10).fit(d, t) for t in ("Ce_VI", "Ce_total")}
    hyper = sim.corpus_hyperpriors(np.random.default_rng(seed))
    rng2 = np.random.default_rng(v2.SEEDS["design"] + k)
    rng = np.random.default_rng(SEEDS3["seq"] + k)
    fams = list(sim.FAMILIES)
    out = []
    for j in range(n_new):
        m = sim.random_material(rng2, fams[j % len(fams)], 9000 + j, hyper, "E")
        for tname, (C0, l6, lt) in v2.TARGETS.items():
            G = v2._design_frame(m, C0)
            best, Gp = v2._choose(G, cps, C0, l6, lt, "nominal")
            if best is None:
                Gp = Gp.assign(score=np.maximum(Gp.Ce6_pt / l6, Gp.Cet_pt / lt))
                best = Gp.sort_values(["score", "D"]).iloc[0]
            ladder = [D_ for D_ in v2.GRID_D if D_ >= best.D]
            truths = [v2._truth(m, C0, best.pH, D_) for D_ in ladder[:4]]
            noncomp = [(a > l6) or (b > lt) for a, b in truths]
            for rule, (n_rep, kind) in RULES_SEQ.items():
                for scen_name, (rho_e, bc) in {"bias_shared": (0.5, 1.0), "bias_independent": (0.0, 0.0)}.items():
                    acc, fa, reps_used = 0, 0, []
                    for _ in range(n_mc):
                        r, false_acc, nr = _sequence(rng, truths, (l6, lt), n_rep, kind, rho_e, bc)
                        acc += r is not None; fa += false_acc; reps_used.append(nr)
                    out.append({"cohort": k, "material": m.material_id, "target": tname, "rule": rule, "bias_model": scen_name,
                                "pH": best.pH, "start_D": best.D, "any_noncompliant_on_ladder": any(noncomp),
                                "first_compliant_round": next((i + 1 for i, nc in enumerate(noncomp) if not nc), None),
                                "p_accept": acc / n_mc, "p_false_accept": fa / n_mc, "mean_batches": float(np.mean(reps_used))})
    return out


def stage_seqver(quick=False):
    n_c, n_new, n_mc = (2, 4, 100) if quick else (20, 15, 1000)
    with Pool(NPROC) as pool:
        R = pd.DataFrame([x for c in pool.map(_seqver_worker, [(k, quick, n_new, n_mc) for k in range(n_c)]) for x in c])
    R.to_csv(RES / "seqver.csv", index=False)
    rows = []
    for (t, rule, bm), g in R.groupby(["target", "rule", "bias_model"]):
        risky = g[g.any_noncompliant_on_ladder]
        rows.append({"target": t, "rule": rule, "bias_model": bm, "n_materials": len(g),
                     "mean_p_false_accept_all": g.p_false_accept.mean(),
                     "mean_p_false_accept_risky": risky.p_false_accept.mean() if len(risky) else np.nan,
                     "max_p_false_accept": g.p_false_accept.max(), "n_risky": len(risky),
                     "frac_materials_pfa_gt_0.05": float((g.p_false_accept > 0.05).mean()),
                     "mean_p_accept": g.p_accept.mean(), "mean_batches": g.mean_batches.mean()})
    S = pd.DataFrame(rows)
    S.to_csv(RES / "T_seqver_summary.csv", index=False)
    dump("seqver", {"n_mc_per_material": n_mc, "summary": S.round(5).to_dict("records")})
    print("seqver done")


# ============================================================================= figures
def figures():
    # Fig 3 (v3): identifiability by design
    F = pd.read_csv(RES / "ident3_fim_detection.csv")
    Q = pd.read_csv(RES / "equiv3.csv.gz")
    S = pd.read_csv(RES / "T_equiv3_summary.csv")
    fig, ax = plt.subplots(2, 2, figsize=(7.4, 6.0), constrained_layout=True)
    bins = np.array([0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0]); mids = 0.5 * (bins[1:] + bins[:-1])
    cut = pd.cut(F.share_red, bins)
    lab = {"D1": "D1 Cr(VI) eq. (C)", "D2": "D2 +kinetics (C)", "D3": "D3 +total Cr (B)", "D4x": "D4x +XPS (A1)",
           "D4n": "D4n +XANES (A1)", "D5": "D5 +inventory/phase (A2)"}
    for i, d in enumerate(idf.DESIGNS_V3):
        ax[0, 0].plot(mids, F.groupby(cut, observed=False)[f"post_sd_red_{d}"].median().values, "-o", ms=3, color=C[i % 6], label=lab[d])
    ax[0, 0].plot(mids, F.groupby(cut, observed=False)["prior_sd_red_D1"].median().values, ":", color=MUTED, label="prior")
    ax[0, 0].set(yscale="log", xlabel="True reduction share", ylabel="Posterior SD of reduction share", title="Fisher analysis (model structure known)")
    ax[0, 0].legend(fontsize=6.2, loc="lower right")
    panel(ax[0, 0], "a")
    for i, (k, l) in enumerate((("dissolved_cr3", "dissolved Cr(III) z-test (B)"), ("xps_cr3", "XPS Cr(III) > 0.20 (A1)"), ("xanes_cr3", "XANES Cr(III) > 0.10 (A1)"))):
        ax[0, 1].plot(mids, F.groupby(cut, observed=False)[k].mean().values, "-o", ms=3, color=C[i], label=f"{l}; controls {F[f'ctrl_{k}'].mean():.2f}")
    ax[0, 1].set(xlabel="True reduction share", ylabel="Detection rate", ylim=(0, 1.02), title="Model-free diagnostics (false alarms in legend)")
    ax[0, 1].legend(fontsize=6.2)
    panel(ax[0, 1], "b")
    s10 = S[S.sd_rel == 0.10].set_index("design").loc[list(idf.DESIGNS_V3)]
    x = np.arange(len(idf.DESIGNS_V3)); w = 0.38
    for j, grp in enumerate((">30%", "5-30%")):
        bottom = np.zeros(len(x))
        for kk, (outc, col) in enumerate((("rejected", C[1]), ("inconclusive", MUTED), ("equivalent", C[2]))):
            v = s10[f"{outc}|{grp}|eps=0.1"].values.astype(float)
            ax[1, 0].bar(x + (j - 0.5) * w, v, w * 0.92, bottom=bottom, color=col, alpha=1 if j == 0 else 0.55,
                         hatch=None if j == 0 else "//", edgecolor="white", label=f"{outc}" if j == 0 else None)
            bottom += v
    ax[1, 0].set_xticks(x, list(idf.DESIGNS_V3))
    ax[1, 0].set(ylabel="Fraction of studies", ylim=(0, 1.0), title="Adsorption-only: 10% noise, ε = 10% (solid >30%, hatched 5–30%)")
    ax[1, 0].title.set_fontsize(8)
    ax[1, 0].legend(fontsize=6.5, loc="upper right")
    panel(ax[1, 0], "c")
    tru = Q[(Q.arm == "true") & (Q.sd_rel == 0.10)]
    thr = S[S.sd_rel == 0.10].set_index("design").calibrated_threshold
    b2 = np.array([0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0]); m2 = 0.5 * (b2[1:] + b2[:-1])
    for i, d in enumerate(idf.DESIGNS_V3):
        g = tru[tru.design == d]
        pw = (g.chi2 > thr[d]).groupby(pd.cut(g.share_red, b2), observed=False).mean()
        ax[1, 1].plot(m2, pw.values, "-o", ms=3, color=C[i % 6], label=d)
    ax[1, 1].axhline(0.05, color=MUTED, ls=":", lw=1)
    ax[1, 1].set(xlabel="True reduction share", ylabel="Rejection rate (calibrated, α = 0.05)", ylim=(0, 1.02),
                 title="Power of the calibrated misfit test (10% noise)")
    ax[1, 1].legend(fontsize=6.5, ncol=2)
    panel(ax[1, 1], "d")
    fig.savefig(FIG / "fig3_identifiability.png"); plt.close(fig)

    # Fig 8: ML across mechanisms + learning curve + censoring
    M = pd.read_csv(RES / "mlscen_metrics.csv")
    U = pd.read_csv(RES / "mlscen_uq.csv")
    L = pd.read_csv(RES / "lc3.csv")
    Cn = pd.read_csv(RES / "T_censor_summary.csv")
    fig, ax = plt.subplots(1, 3, figsize=(7.6, 3.0), constrained_layout=True)
    sc = list(sim.SCENARIOS); x = np.arange(len(sc)); w = 0.2
    for i, (scheme, col) in enumerate((("random", MUTED), ("LOSO", C[0]), ("LOFO", C[2]), ("transfer_from_E", C[1]))):
        g = M[(M.model == "RF") & (M.features == "chem") & (M.scheme == scheme)].groupby("scenario").R2
        mu = g.mean().reindex(sc); sdv = g.std().reindex(sc)
        ax[0].bar(x + (i - 1.5) * w, mu.values, w * 0.92, color=col, label={"transfer_from_E": "trained on E"}.get(scheme, scheme))
        ax[0].errorbar(x + (i - 1.5) * w, mu.values, yerr=sdv.values, fmt="none", ecolor=INK, lw=0.7, capsize=1.5)
    ax[0].set_xticks(x, sc)
    ax[0].set(xlabel="Mechanism scenario", ylabel="R$^2$ (RF, +spec.)", ylim=(-0.2, 1), title="Across mechanisms (5 corpora each)")
    ax[0].legend(fontsize=6.3, ncol=2, loc="lower left")
    panel(ax[0], "a")
    g = L.groupby("n_studies")
    ax[1].errorbar(g.R2_LOSO.mean().index, g.R2_LOSO.mean(), yerr=g.R2_LOSO.std() / np.sqrt(g.size()), fmt="-o", ms=3.5, color=C[0], capsize=2, label="LOSO (mean ± SE)")
    ax[1].scatter(L.n_studies, L.R2_LOSO, s=7, color=C[0], alpha=0.3)
    ax[1].errorbar(g.R2_random.mean().index, g.R2_random.mean(), fmt="--s", ms=3.5, color=MUTED, label="random")
    ax[1].set(xscale="log", xlabel="Number of studies", ylabel="R$^2$ (RF, +spec.)", title="Learning curve (8 replicates)")
    ax[1].xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax[1].set_xticks([12, 24, 48, 96], ["12", "24", "48", "96"])
    ax[1].legend(fontsize=6.5)
    panel(ax[1], "b")
    uq = U.groupby(["scenario", "method"]).coverage.apply(lambda c: (c >= 0.9 - 1e-9).mean()).unstack().reindex(sc)
    for i, meth in enumerate(("M0", "M1", "M2", "M2s")):
        if meth in uq:
            ax[2].plot(x, uq[meth].values, "-o", ms=3.5, color=[C[0], C[1], C[2], C[5]][i], label=meth)
    ax[2].axhline(0.9, color=MUTED, ls=":", lw=1)
    ax[2].set_xticks(x, sc)
    ax[2].set(xlabel="Mechanism scenario", ylabel="Studies with ≥ 90% coverage", ylim=(0.4, 1.02), title="Study-level intervals (90%)")
    ax[2].legend(fontsize=6.5)
    panel(ax[2], "c")
    fig.savefig(FIG / "fig8_mechanism_transfer.png"); plt.close(fig)

    # Fig 5 (v3): conformal incl. M2s
    U2 = pd.read_csv(v2.RES / "uq_per_study.csv")
    S2 = pd.read_csv(RES / "m2s_per_study.csv")
    fig, ax = plt.subplots(1, 3, figsize=(7.6, 2.9), constrained_layout=True)
    u = U2[U2.target == "q_app_VI"]
    names = {"M0": C[0], "M1": C[1], "M2": C[2]}
    for meth, col in names.items():
        gg = u[u.method == meth].groupby("alpha").coverage.apply(lambda c: (c >= 0.9 - 1e-9).mean())
        ax[0].plot(1 - gg.index.values, gg.values, "-o", ms=3.5, color=col, label=meth)
    s = S2[S2.target == "q_app_VI"]
    gg = s.groupby("alpha").coverage.apply(lambda c: (c >= 0.9 - 1e-9).mean())
    ax[0].plot(1 - gg.index.values, gg.values, "-s", ms=3.5, color=C[5], label="M2s (split, exact)")
    ax[0].plot([0.45, 1], [0.45, 1], color=MUTED, lw=1, ls=":")
    ax[0].set(xlabel="Nominal 1 − α", ylabel="Studies with ≥ 90% coverage", title="Study-level guarantee, q$_{app,VI}$")
    ax[0].legend(fontsize=6.5)
    panel(ax[0], "a")
    g90 = u[u.alpha == 0.1]
    data = [g90[g90.method == m_].coverage.values for m_ in names] + [s[s.alpha == 0.1].coverage.values]
    ax[1].boxplot(data, tick_labels=["M0", "M1", "M2", "M2s"], widths=0.55, medianprops={"color": INK}, flierprops={"markersize": 2.5})
    ax[1].axhline(0.9, color=MUTED, ls=":", lw=1)
    ax[1].set(ylabel="Per-study coverage at 90%", title="Spread across studies")
    panel(ax[1], "b")
    tl = {"q_app_VI": "q$_{app,VI}$", "q_tot_aq": "q$_{tot,aq}$", "Ce_VI": "C$_{e,VI}$", "Ce_total": "C$_{e,tot}$"}
    for i, tgt in enumerate(tl):
        w_ = [U2[(U2.target == tgt) & (U2.alpha == 0.1) & (U2.method == m_)].width_log10.mean() for m_ in names]
        ss = S2[(S2.target == tgt) & (S2.alpha == 0.1)]
        w_.append(ss[~ss.unbounded].width_log10.mean() if (~ss.unbounded).any() else np.nan)
        ax[2].bar(np.arange(4) + (i - 1.5) * 0.2, w_, 0.19, color=C[i], label=tl[tgt])
    ax[2].set_xticks(range(4), ["M0", "M1", "M2", "M2s"])
    ax[2].set(ylabel="Mean 90% width (log$_{10}$)", title="Width (bounded intervals)")
    ax[2].legend(fontsize=6.3, ncol=2)
    panel(ax[2], "c")
    fig.savefig(FIG / "fig5_conformal.png"); plt.close(fig)

    # Fig 6 (v3): add sequential verification
    O = pd.read_csv(v2.RES / "design_outcomes.csv")
    O["met_both"] = O["met_both"].fillna(False).astype(bool)
    R = pd.read_csv(RES / "seqver.csv")
    fig, ax = plt.subplots(1, 3, figsize=(7.6, 3.1), constrained_layout=True, gridspec_kw={"width_ratios": [1.7, 1, 1]})
    rules = ["nominal", "marginal90", "joint90_bonf", "joint90_adapt_k4", "joint80_adapt_k4", "nominal_verify"]
    rl = {"nominal": "point", "marginal90": "marg.\n90%", "joint90_bonf": "joint\n90%", "joint90_adapt_k4": "J90\n+4b",
          "joint80_adapt_k4": "J80\n+4b", "nominal_verify": "point\n+verify"}
    xx = np.arange(len(rules)); w = 0.2
    for j, (t, hatch) in enumerate((("S1", None), ("S2", "//"))):
        uns, found = [], []
        for r_ in rules:
            g = O[(O.target == t) & (O.rule == r_)]
            uns.append((g.feasible & ~g.met_both).mean()); att = g[g.attainable]; found.append((att.feasible & att.met_both).mean())
        off = -1.5 * w if j == 0 else 0.5 * w
        ax[0].bar(xx + off, uns, w * 0.95, color=C[1], hatch=hatch, alpha=1 if j == 0 else 0.55, edgecolor="white", label=f"unsafe, {t}")
        ax[0].bar(xx + off + w, found, w * 0.95, color=C[2], hatch=hatch, alpha=1 if j == 0 else 0.55, edgecolor="white", label=f"found, {t}")
    ax[0].set_xticks(xx, [rl[r_] for r_ in rules], fontsize=7)
    ax[0].set(ylabel="Fraction of materials", ylim=(0, 1), title="Design rules (simulated cohorts)")
    ax[0].legend(fontsize=6.3, ncol=2, loc="upper left")
    panel(ax[0], "a")
    for i, (rule, col) in enumerate((("v2_ratio_n2", C[1]), ("ucb_n3_bonf", C[0]))):
        g = R[(R.rule == rule) & (R.bias_model == "bias_shared") & R.any_noncompliant_on_ladder]
        ax[1].hist(g.p_false_accept, bins=np.linspace(0, 0.5, 26), color=col, alpha=0.6,
                   label={"v2_ratio_n2": "mean ≤ 0.9·limit, n = 2", "ucb_n3_bonf": "UCB ≤ limit, n = 3"}[rule])
    ax[1].axvline(0.05, color=INK, ls=":", lw=1)
    ax[1].set(xlabel="P(false acceptance) per material", ylabel="Materials (ladder with a non-compliant dose)", title="Sequential verification")
    ax[1].legend(fontsize=6.3)
    panel(ax[1], "b")
    Sg = R.groupby(["rule", "target"]).agg(pa=("p_accept", "mean"), nb=("mean_batches", "mean")).reset_index()
    for i, (rule, col) in enumerate((("v2_ratio_n2", C[1]), ("ucb_n3_bonf", C[0]))):
        g = Sg[Sg.rule == rule]
        ax[2].bar(np.arange(2) + (i - 0.5) * 0.36, g.pa.values, 0.34, color=col, label={"v2_ratio_n2": "ratio, n = 2", "ucb_n3_bonf": "UCB, n = 3"}[rule])
        for x_, nb_ in zip(np.arange(2) + (i - 0.5) * 0.36, g.nb.values):
            ax[2].text(x_, 0.03, f"{nb_:.1f}", ha="center", fontsize=6.5, color="white")
    ax[2].set_xticks(range(2), ["S1", "S2"])
    ax[2].set(ylabel="P(acceptance within 4 rounds)", ylim=(0, 1), title="Cost (batches in bar)")
    ax[2].legend(fontsize=6.3)
    panel(ax[2], "c")
    fig.savefig(FIG / "fig6_design.png"); plt.close(fig)
    print("figures v3 done")


STAGES = {"chem3": stage_chem3, "ident3": stage_ident3, "equiv3": stage_equiv3, "m2s": stage_m2s,
          "mlscen": stage_mlscen, "censor": stage_censor, "lc3": stage_lc3, "seqver": stage_seqver}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all", choices=list(STAGES) + ["all", "figures"])
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    todo = list(STAGES) if a.stage == "all" else ([] if a.stage == "figures" else [a.stage])
    for s in todo:
        t = time.time()
        STAGES[s](a.quick)
        print(f"[{s}] {time.time() - t:.1f} s", flush=True)
    if a.stage in ("all", "figures"):
        figures()
