"""CASSIA v0.2 experiment driver: regenerates every table, figure and number of the revised paper.

    python scripts/run_v2.py --stage all            # full run (about 1.5-2 h on 4 cores)
    python scripts/run_v2.py --stage bench --quick  # one stage, reduced size (smoke test)

Stages: chem, bench, ident, ml, uq, design, missing, figures. Outputs: results/v2/.
All randomness is seeded; stage seeds are listed in SEEDS.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):  # one BLAS thread per worker
    os.environ.setdefault(_v, "1")

import argparse
import json
import time
from multiprocessing import Pool
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from cassia import chemistry as ch, dataset as ds, fitting as fit, identifiability as idf, simulate as sim
from cassia import validation as val

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results" / "v2"
FIG = ROOT / "manuscript" / "figures"
RES.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
NPROC = 4
SEEDS = {"bench": 100000, "ident": 200000, "equiv": 300000, "ml": 5000, "lc": 9000, "ext": 5000, "design": 7000}

C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"]
INK, MUTED, GRID = "#0b0b0b", "#6e6c66", "#e1e0d9"
plt.rcParams.update({"font.size": 10, "axes.titlesize": 10, "axes.labelsize": 10, "xtick.labelsize": 9,
                     "ytick.labelsize": 9, "legend.fontsize": 8.5, "axes.edgecolor": "#a8a69c",
                     "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True,
                     "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "lines.linewidth": 1.8, "savefig.dpi": 300,
                     "figure.facecolor": "white", "legend.frameon": False, "font.family": "DejaVu Sans"})
SC_COL = {"A": MUTED, "B": C[0], "C": C[2], "D": C[3], "E": C[1]}


def panel(ax, letter):
    ax.set_title(f"({letter})  " + ax.get_title(), loc="left", fontsize=10)
    ax.set_title("")


def _clean(o):
    if isinstance(o, dict):
        return {("|".join(map(str, k)) if isinstance(k, tuple) else str(k)): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def dump(name, obj):
    (RES / f"{name}.json").write_text(json.dumps(_clean(obj), indent=2, default=str))


def material_from(meta):
    return sim.Material(**{k: meta[k] for k in sim.Material.__dataclass_fields__})


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (c - h, c + h)


# =============================================================================
# Stage: chemistry
# =============================================================================
def stage_chem(quick=False):
    out = {}
    prov = ch.database_provenance()
    prov.to_csv(RES / "T_db_provenance.csv", index=False)
    (ch.DB_DIR / "SHA256SUMS").write_text("".join(f"{r.sha256}  {r.database}\n" for r in prov.itertuples()))
    out["provenance"] = prov.to_dict("records")
    lk = ch.logk_table()
    lk.to_csv(RES / "T_logK.csv")
    out["logK"] = lk.round(4).to_dict()
    pH = np.round(np.arange(1.5, 10.01, 0.1 if not quick else 0.5), 2)
    tabs = {}
    for c0 in (10.0, 500.0):
        t = ch.speciation_table(pH, cr_mgL=c0, databases=ch.DATABASES)
        t.to_csv(RES / f"speciation_C0_{int(c0)}.csv", index=False)
        tabs[c0] = t
    disc = {}
    for c0, t in tabs.items():
        ref = t[t.database == "minteq.v4.dat"].set_index("pH")
        for db in ("minteq.dat", "llnl.dat"):
            o = t[t.database == db].set_index("pH")
            cols = [c for c in ("f_HCrO4-", "f_CrO4-2", "f_Cr2O7-2", "f_H2CrO4") if c in ref and c in o]
            d = (ref[cols] - o[cols]).abs()
            disc[f"{db}_vs_minteq.v4_C0_{int(c0)}"] = {c: [float(d[c].max()), float(d[c].idxmax())] for c in cols}
    out["speciation_max_abs_diff"] = disc
    cross = {}
    t10 = tabs[10.0]
    for db in ch.DATABASES:
        s = t10[t10.database == db]
        diff = (s["f_HCrO4-"] - s["f_CrO4-2"]).values
        cross[db] = float(np.interp(0, -diff, s["pH"].values))
    out["crossover_pH_10mgL"] = cross
    g = ch.grid()
    p = np.linspace(5.5, 7.2, 200)
    cI = {}
    for I in (0.001, 0.01, 0.1):
        f = g(p, 10.0, I, 298.15)
        cI[str(I)] = float(np.interp(0, -(f[:, 0] - f[:, 1]), p))
    cT = {}
    for T in (288.15, 298.15, 318.15):
        f = g(p, 10.0, 0.01, T)
        cT[str(T)] = float(np.interp(0, -(f[:, 0] - f[:, 1]), p))
    out["crossover_vs_I"], out["crossover_vs_T"] = cI, cT
    out["dimer_fraction_pH3"] = {db: {str(c): float(ch.speciate_cr6(3.0, c, database=db).fractions.get("Cr2O7-2", 0))
                                      for c in (1, 10, 50, 100, 300, 700)} for db in ch.DATABASES}
    out["grid_validation"] = g.validate(300 if not quick else 50)
    # redox boundary
    pe_grid = np.arange(-4, 21.01, 0.25)
    thr = []
    for db in ch.DATABASES:
        for ph in np.arange(2, 10.01, 0.5 if not quick else 2):
            pes = pe_grid[pe_grid <= 20.75 - ph]
            f6 = np.array([ch.redox_scenario(ph, 50.0, pe, database=db)["frac_CrVI"] for pe in pes])
            crossed = f6.max() >= 0.5
            pe50 = float(np.interp(0.5, f6, pes)) if crossed else np.nan
            thr.append({"database": db, "pH": ph, "pe50": pe50, "Eh50_V": pe50 * 0.05916,
                        "frac_CrVI_at_water_limit": float(f6[-1])})
    thr = pd.DataFrame(thr)
    thr.to_csv(RES / "T_redox_boundary.csv", index=False)
    piv = thr.pivot(index="pH", columns="database", values="Eh50_V")
    out["Eh50_V"] = piv.round(3).to_dict()
    d1 = (piv["minteq.v4.dat"] - piv["llnl.dat"]).dropna() * 1000
    d2 = (piv["minteq.v4.dat"] - piv["minteq.dat"]).dropna() * 1000
    out["Eh50_diff_mV"] = {"minteq.v4-llnl": [float(d1.min()), float(d1.max())],
                           "minteq.v4-minteq": [float(d2.min()), float(d2.max())]}
    # Cr(III) solubility by direct SI bisection
    rows = []
    for db in ("minteq.v4.dat", "minteq.dat"):
        for I in (0.001, 0.01, 0.1):
            for ph in np.arange(3.0, 10.01, 0.25 if not quick else 1.0):
                rows.append(ch.cr3_solubility_si(float(ph), db, I))
    sol = pd.DataFrame(rows)
    sol.to_csv(RES / "T_cr3_solubility_SI.csv", index=False)
    out["cr3_solubility_status_counts"] = sol.groupby(["database", "status"]).size().to_dict()
    s01 = sol[sol.I == 0.01].pivot(index="pH", columns="database", values="S_mgL")
    out["cr3_S_mgL_I0.01"] = {str(k): v for k, v in s01.loc[[x for x in (3.0, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0) if x in s01.index]].to_dict().items()}
    out["llnl_cr3_solubility"] = "not computable: llnl.dat has no Cr(III) hydroxide; eskolaite is written via CrO4-2 + O2"
    # SI map
    if not quick:
        phs = np.arange(3.0, 9.01, 0.25)
        cs = np.logspace(-3, 3, 25)
        SI = np.full((len(phs), len(cs)), np.nan)
        conv = 0
        for i, ph in enumerate(phs):
            for j, c in enumerate(cs):
                v, ok = ch.si_cr3(float(ph), float(c))
                SI[i, j] = v if ok else np.nan
                conv += ok
        np.savez_compressed(RES / "SI_map_minteq_v4.npz", SI=SI, pH=phs, C=cs)
        out["SI_map_converged_fraction"] = conv / SI.size
    dump("chem", out)
    print("chem done")


# =============================================================================
# Stage: isotherm-fit benchmark across mechanisms (>= 100 replicate corpora)
# =============================================================================
def _bench_worker(args):
    scen, rep, n_studies, boot = args
    seed = SEEDS["bench"] + 1000 * sim.SCENARIOS.index(scen) + rep
    df, mats = sim.virtual_corpus(n_studies, seed=seed, scenario=scen)
    mats = mats.set_index("study_id")
    out = []
    for sid, gdf in df.groupby("study_id"):
        m = mats.loc[sid]
        iso = ds.endpoints(gdf[gdf.series == "isotherm"]).sort_values("Ce_VI")
        fits = fit.fit_isotherms(iso["Ce_VI"].values, iso["q_app_VI"].values, T=float(m["T"]))
        if "Langmuir" not in fits:
            continue
        sel = fit.selection_table(fits)
        lang = fits["Langmuir"]
        rec = {"scenario": scen, "rep": rep, "study_id": sid, "family": m["family"], "evidence": m["study_evidence"],
               "mu_R_family": m["mu_R_family"],
               "best_iso": sel.loc[0, "model"], "R2_best": float(sel["R2"].max()), "R2_L": float(lang.r2),
               "runs_p_L": float(sel.set_index("model").loc["Langmuir", "runs_p"]),
               "qmax_L": float(lang.params[1]), "qmax_true": float(m["qmax"]),
               "q_ads_true_maxC0": float(iso["q_adsVI_true"].iloc[-1]),
               "q_app_maxC0": float(iso["q_app_VI"].iloc[-1])}
        rem = (iso["C0"] - iso["Ce_VI_true"]).clip(lower=1e-12)
        rec["share_red"] = float(iso["reduced_mgL"].sum() / rem.sum())
        rec["share_precip"] = float(iso["Cr3_precip_true"].sum() / rem.sum())
        rec["share_cr3_aq"] = float(iso["Cr3_aq_true"].sum() / rem.sum())
        rec["share_cr3_surface"] = float(iso["Cr3_surface_true"].sum() / rem.sum())
        rec["share_ads"] = 1 - rec["share_red"]
        if "Ce_total" in iso and iso["Ce_total"].notna().all():
            rec["aq_cr3_index"] = float((iso["Ce_total"] - iso["Ce_VI"]).clip(lower=0).sum() / (iso["C0"] - iso["Ce_VI"]).sum())
        if boot:
            ci, _ = fit.bootstrap_ci(lang, iso["Ce_VI"].values, T=float(m["T"]), B=200, seed=rep)
            rec["qmax_ci_lo"], rec["qmax_ci_hi"] = float(ci[0, 1]), float(ci[1, 1])
            rec["qmax_true_in_ci"] = bool(ci[0, 1] <= m["qmax"] <= ci[1, 1])
        kin = ds.endpoints(gdf[gdf.series == "kinetic"]).sort_values("t")
        kt = fit.selection_table(fit.fit_kinetics(kin["t"].values, kin["q_app_VI"].values))
        rec["best_kin"] = kt.loc[0, "model"] if len(kt) else None
        out.append(rec)
    return out


def stage_bench(quick=False):
    n_rep, n_st, n_boot = (4, 12, 2) if quick else (100, 30, 20)
    jobs = [(s, r, n_st, r < n_boot) for s in sim.SCENARIOS for r in range(n_rep)]
    with Pool(NPROC) as pool:
        res = [x for chunk in pool.imap_unordered(_bench_worker, jobs, chunksize=2) for x in chunk]
    B = pd.DataFrame(res)
    B["qmax_ratio"] = B["qmax_L"] / B["qmax_true"]
    B.to_csv(RES / "bench_studies.csv.gz", index=False)
    summ = {"n_replicates_per_scenario": n_rep, "n_studies_per_corpus": n_st, "n_bootstrap_corpora": n_boot,
            "n_study_fits": len(B)}
    for s, g in B.groupby("scenario"):
        pub = g[g.R2_best >= 0.95]
        per_rep = g.groupby("rep").apply(lambda x: pd.Series({
            "share_R2ge95": (x.R2_best >= 0.95).mean(),
            "rho": stats.spearmanr(x.share_red, x.qmax_ratio).statistic if x.share_red.std() > 0 else np.nan}),
            include_groups=False)
        hi = g[g.share_red > 0.2]
        d = {"n": len(g),
             "share_red_quantiles": g.share_red.quantile([0.05, 0.25, 0.5, 0.75, 0.95]).round(4).tolist(),
             "share_precip_quantiles": g.share_precip.quantile([0.05, 0.5, 0.95]).round(4).tolist(),
             "share_cr3_aq_quantiles": g.share_cr3_aq.quantile([0.05, 0.5, 0.95]).round(4).tolist(),
             "share_cr3_surface_quantiles": g.share_cr3_surface.quantile([0.05, 0.5, 0.95]).round(4).tolist(),
             "frac_R2best_ge_0.95": float((g.R2_best >= 0.95).mean()),
             "frac_R2best_ge_0.95_mcse": float(per_rep.share_R2ge95.std() / np.sqrt(len(per_rep))),
             "n_red_gt_0.2": len(hi),
             "frac_of_red_gt_0.2_with_R2_ge_0.95": float((hi.R2_best >= 0.95).mean()) if len(hi) else np.nan,
             "pub_median_qmax_ratio": float(pub.qmax_ratio.median()),
             "pub_qmax_ratio_IQR": pub.qmax_ratio.quantile([0.25, 0.75]).tolist(),
             "pub_median_qmax_ratio_red_gt_0.3": float(pub[pub.share_red > 0.3].qmax_ratio.median()) if (pub.share_red > 0.3).any() else np.nan,
             "pub_median_qmax_ratio_red_lt_0.05": float(pub[pub.share_red < 0.05].qmax_ratio.median()) if (pub.share_red < 0.05).any() else np.nan,
             "spearman_red_vs_qmax_ratio_mean_sd": [float(per_rep.rho.mean()), float(per_rep.rho.std())],
             "runs_flag_rate": float((g.runs_p_L < 0.05).mean()),
             "best_kin_PSO_rate": float((g.best_kin == "PSO").mean())}
        if "qmax_true_in_ci" in g:
            bb = g[g.qmax_true_in_ci.notna()]
            bp = bb[bb.R2_best >= 0.95]
            d["bootstrap_ci_coverage_of_true_qmax_pub"] = float(bp.qmax_true_in_ci.astype(bool).mean()) if len(bp) else np.nan
            d["n_bootstrap_pub"] = len(bp)
        if g["aq_cr3_index"].notna().any():
            a = g[g.aq_cr3_index.notna() & (g.share_red > 0.02)]
            d["aq_index_spearman"] = float(stats.spearmanr(a.aq_cr3_index, a.share_red).statistic) if len(a) > 3 else np.nan
            d["aq_index_over_true_median"] = float((a.aq_cr3_index / a.share_red).median()) if len(a) else np.nan
        summ[s] = d
    # circularity check: is the most reductive family the same across replicates?
    E = B[B.scenario == "E"]
    top = E.groupby(["rep", "family"]).share_red.mean().reset_index().sort_values("share_red").groupby("rep").tail(1)
    summ["E_most_reductive_family_counts"] = top.family.value_counts().to_dict()
    dump("bench", summ)
    print("bench done", len(B))


# =============================================================================
# Stage: identifiability (FIM, detection, metadata bias, equivalence, profiles)
# =============================================================================
def _ident_worker(args):
    scen, rep, n_studies = args
    seed = SEEDS["ident"] + 1000 * sim.SCENARIOS.index(scen) + rep
    _, mats = sim.virtual_corpus(n_studies, seed=seed, scenario=scen)
    rng = np.random.default_rng(seed)
    out = []
    for _, meta in mats.iterrows():
        meta = meta.to_dict()
        m = material_from(meta)
        try:
            fs = idf.fisher_study(m, meta)
        except RuntimeError:
            continue
        rec = {"scenario": scen, "rep": rep, "study_id": meta["study_id"], "family": meta["family"],
               "share_red": fs["base"]["share_red"], "share_precip": fs["base"]["share_precip"]}
        for d in idf.DESIGNS:
            r = idf.posterior_share_sd(fs, d)
            rec[f"post_sd_{d}"], rec[f"prior_sd_{d}"], rec[f"cond_{d}"] = r["post_sd_share"], r["prior_sd_share"], r["fim_cond"]
        rec.update(idf.detection_tests(rng, fs, meta, n_rep=20))
        for d in ("D1", "D3"):
            for k, v in idf.metadata_bias(m, meta, fs, design=d).items():
                rec[f"{k}_{d}"] = v
        out.append(rec)
    return out


def _equiv_worker(args):
    scen, rep, n_studies = args
    seed = SEEDS["equiv"] + 1000 * sim.SCENARIOS.index(scen) + rep
    _, mats = sim.virtual_corpus(n_studies, seed=seed, scenario=scen)
    out = []
    for k, (_, meta) in enumerate(mats.iterrows()):
        meta = meta.to_dict()
        m = material_from(meta)
        for lab, mm in (("true", m), ("control_A", sim.apply_scenario(m, "A"))):
            for sd in (0.03, 0.10):
                for d in idf.DESIGNS:
                    try:
                        r = idf.equivalence_test(mm, meta, d, sd_rel=sd, seed=seed + k, n_starts=3)
                    except (RuntimeError, ValueError):
                        continue
                    out.append({"scenario": scen, "rep": rep, "study_id": meta["study_id"], "arm": lab, **r})
    return out


def _profile_worker(args):
    meta, design, grid = args
    m = material_from(meta)
    df = idf.profile_likelihood(m, meta, design, grid, seed=1, max_nfev=PROFILE_NFEV)
    df["study_id"] = meta["study_id"]
    df["true_log10_R0"] = np.log10(meta["R0"])
    df["true_share"] = idf.observe(m, meta)["share_red"]
    return df


PROFILE_NFEV = 150


def stage_profile(quick=False):
    """Profiles only (re-run with a larger optimisation budget and two starts per point)."""
    F = pd.read_csv(RES / "ident_fim_detection.csv")
    n_st = 6 if quick else 25
    E = F[F.scenario == "E"]
    picks = [E.iloc[(E.share_red - 0.15).abs().argmin()], E.iloc[(E.share_red - 0.6).abs().argmin()]]
    pj = []
    for pk in picks:
        _, mats = sim.virtual_corpus(n_st, seed=SEEDS["ident"] + 1000 * sim.SCENARIOS.index("E") + int(pk.rep), scenario="E")
        meta = mats.set_index("study_id").loc[pk.study_id].to_dict() | {"study_id": pk.study_id}
        lr0 = np.log10(meta["R0"])
        for d in ("D1", "D3", "D4"):
            pj.append((meta, d, np.linspace(lr0 - 1.2, lr0 + 1.2, 5 if quick else 13)))
    with Pool(2) as pool:
        P = pd.concat(pool.map(_profile_worker, pj), ignore_index=True)
    P.to_csv(RES / "ident_profiles.csv", index=False)
    summ = json.loads((RES / "ident.json").read_text())
    summ["profiles"] = {}
    for (s_, d), g in P.groupby(["study_id", "design"]):
        ok = g[g.delta <= 3.84]
        summ["profiles"][f"{s_}|{d}"] = {
            "true_share": float(g.true_share.iloc[0]), "n_obs": int(g.n_obs.iloc[0]),
            "min_m2logL": float(g.m2logL.min()), "share_at_min": float(g.loc[g.delta.idxmin(), "share_red"]),
            "share_range_within_95": [float(ok.share_red.min()), float(ok.share_red.max())],
            "log10R0_range_within_95_rel_true": [float((ok.log10_R0 - ok.true_log10_R0).min()), float((ok.log10_R0 - ok.true_log10_R0).max())],
            "frac_points_converged": float(g.converged.mean())}
    dump("ident", summ)
    print("profile done")


def stage_ident(quick=False):
    n_st = 6 if quick else 25
    reps = {"E": 2 if quick else 8, "B": 1 if quick else 2, "C": 1 if quick else 2, "D": 1 if quick else 2}
    jobs = [(s, r, n_st) for s, n in reps.items() for r in range(n)]
    with Pool(NPROC) as pool:
        res = [x for chunk in pool.imap_unordered(_ident_worker, jobs) for x in chunk]
    F = pd.DataFrame(res)
    F.to_csv(RES / "ident_fim_detection.csv", index=False)
    # equivalence
    ereps = {"E": 1 if quick else 4, "B": 0 if quick else 1, "C": 0 if quick else 1, "D": 0 if quick else 1}
    n_eq = 3 if quick else 25
    ejobs = [(s, r, n_eq) for s, n in ereps.items() for r in range(n)]
    with Pool(NPROC) as pool:
        eq = [x for chunk in pool.imap_unordered(_equiv_worker, ejobs) for x in chunk]
    Q = pd.DataFrame(eq)
    Q.to_csv(RES / "ident_equivalence.csv", index=False)
    # profiles: one low- and one high-reduction scenario-E study
    E = F[F.scenario == "E"]
    picks = [E.iloc[(E.share_red - 0.15).abs().argmin()], E.iloc[(E.share_red - 0.6).abs().argmin()]]
    pj = []
    for pk in picks:
        _, mats = sim.virtual_corpus(n_st, seed=SEEDS["ident"] + 1000 * sim.SCENARIOS.index("E") + int(pk.rep), scenario="E")
        meta = mats.set_index("study_id").loc[pk.study_id].to_dict() | {"study_id": pk.study_id}
        lr0 = np.log10(meta["R0"])
        grid = np.linspace(lr0 - 1.2, lr0 + 1.2, 5 if quick else 11)
        for d in ("D1", "D3", "D4"):
            pj.append((meta, d, grid))
    with Pool(NPROC) as pool:
        P = pd.concat(pool.map(_profile_worker, pj), ignore_index=True)
    P.to_csv(RES / "ident_profiles.csv", index=False)
    # summaries
    summ = {"n_fim_studies": F.groupby("scenario").size().to_dict()}
    bins = [0, 0.05, 0.2, 0.4, 1.0]
    F["red_bin"] = pd.cut(F.share_red, bins, include_lowest=True)
    for d in idf.DESIGNS:
        F[f"ratio_{d}"] = F[f"post_sd_{d}"] / F[f"prior_sd_{d}"]
    summ["post_sd_median_by_design"] = {d: float(F[f"post_sd_{d}"].median()) for d in idf.DESIGNS}
    summ["post_sd_p90_by_design"] = {d: float(F[f"post_sd_{d}"].quantile(0.9)) for d in idf.DESIGNS}
    summ["sd_ratio_median_by_design"] = {d: float(F[f"ratio_{d}"].median()) for d in idf.DESIGNS}
    summ["post_sd_by_redbin"] = {d: F.groupby("red_bin", observed=True)[f"post_sd_{d}"].median().round(4).tolist() for d in idf.DESIGNS}
    summ["detection_power_by_redbin"] = {k: F.groupby("red_bin", observed=True)[k].mean().round(3).tolist()
                                         for k in ("lack_of_fit", "dissolved_cr3", "xps_cr3")}
    summ["red_bins"] = [str(b) for b in F.red_bin.cat.categories]
    summ["metadata_bias_pH0.5"] = {d: F[f"bias_share_pH_{d}"].abs().quantile([0.5, 0.9]).round(4).tolist() for d in ("D1", "D3")}
    summ["metadata_bias_Ix10"] = {d: F[f"bias_share_I_{d}"].abs().quantile([0.5, 0.9]).round(4).tolist() for d in ("D1", "D3")}
    eqs = {}
    for (arm, sd, d), g in Q.groupby(["arm", "sd_rel", "design"]):
        gg = g if arm == "control_A" else g[g.share_red > 0.05]
        eqs[f"{arm}|sd={sd}|{d}"] = {"n": len(gg), "frac_indistinguishable": float(gg.indistinguishable.mean()) if len(gg) else np.nan}
    summ["equivalence"] = eqs
    hi = Q[(Q.arm == "true") & (Q.share_red > 0.3)]
    summ["equivalence_red_gt_0.3"] = hi.groupby(["sd_rel", "design"]).indistinguishable.mean().round(3).to_dict()
    summ["equivalence_red_gt_0.3"] = {f"{k[0]}|{k[1]}": v for k, v in summ["equivalence_red_gt_0.3"].items()}
    summ["profiles"] = {f"{s}|{d}": {"true_share": float(g.true_share.iloc[0]),
                                     "share_range_within_delta_3.84": [float(g[g.delta <= 3.84].share_red.min()), float(g[g.delta <= 3.84].share_red.max())]}
                        for (s, d), g in P.groupby(["study_id", "design"])}
    dump("ident", summ)
    print("ident done", len(F), len(Q))


# =============================================================================
# Stage: ML validation across replicate corpora
# =============================================================================
def add_features(df):
    d = ds.endpoints(df)
    d = ch.speciation_features(d, conc_col="C0")
    d["pH_minus_pzc"] = d["pH"] - d["pHpzc"]
    d["I_imputed"] = d["I_imputed"].astype(float)
    return d


def corpus_features(seed, n=30, prefix="VS"):
    df, mats = sim.virtual_corpus(n, seed=seed, scenario="E", prefix=prefix)
    return add_features(df), mats


def _ml_worker(args):
    k, quick, with_gpr = args
    d, _ = corpus_features(SEEDS["ml"] + k, 12 if quick else 30)
    rows, perstudy = [], []
    schemes = ("random", "grouped5", "LOSO", "LOFO")
    models = ["RF", "XGB"] + (["GPR"] if with_gpr else [])
    for mname in models:
        for fs_name, fs in (("base", val.BASE), ("chem", val.CHEM)):
            for sc in (schemes if mname != "GPR" else ("random", "LOSO")):
                r, pl, dd, qp = val.evaluate(d, mname, fs, sc)
                r["corpus"] = k
                rows.append(r)
                if sc == "LOSO" and mname == "RF":
                    ps = val.per_study(dd, np.clip(dd["q_app_VI"].values, 0, None), qp)
                    ps["features"], ps["corpus"] = fs_name, k
                    perstudy.append(ps)
                    if fs_name == "chem":
                        icc = val.between_study_variance(dd, pl)
                        r.update({f"icc_{a}": b for a, b in icc.items()})
                        # naive vs smeared retransformation bias
                        r2, _, _, qn = val.evaluate(d, mname, fs, sc, retransform="naive")
                        r["bias_naive"], r["bias_smear"] = r2["bias"], r["bias"]
    # external: train on this corpus, test on the next corpus (different hyper-prior draw)
    dx, _ = corpus_features(SEEDS["ext"] + 1000 + k, 12 if quick else 30, prefix="EX")
    for mname in ("RF", "XGB"):
        for fs_name, fs in (("base", val.BASE), ("chem", val.CHEM)):
            mdl = val.make_model(mname).fit(d[fs], val.to_target(d["q_app_VI"].values))
            sm = val.duan_smear(val.to_target(d["q_app_VI"].values) - mdl.predict(d[fs]))
            qp = val.from_target(mdl.predict(dx[fs]), sm)
            qt = np.clip(dx["q_app_VI"].values, 0, None)
            rows.append({"model": mname, "features": fs_name, "scheme": "external", "corpus": k,
                         **val.metrics(qt, qp), **val.cluster_bootstrap(dx, qt, qp)})
    return rows, perstudy


def _tune_worker(args):
    k, quick = args
    d, _ = corpus_features(SEEDS["ml"] + k, 12 if quick else 30)
    out = []
    for mname in ("RF", "XGB"):
        for tune in (False, True):
            r, *_ = val.evaluate(d, mname, val.CHEM, "LOSO", tune=tune)
            r["corpus"] = k
            out.append(r)
    return out


def _lc_worker(args):
    n, rep = args
    d, _ = corpus_features(SEEDS["lc"] + 100 * n + rep, n)
    r, *_ = val.evaluate(d, "RF", val.CHEM, "LOSO")
    r2, *_ = val.evaluate(d, "RF", val.CHEM, "random")
    return {"n_studies": n, "rep": rep, "R2_LOSO": r["R2"], "RMSE_LOSO": r["RMSE"], "R2_random": r2["R2"]}


def stage_ml(quick=False):
    n_c = 3 if quick else 20
    n_gpr = 1 if quick else 3
    jobs = [(k, quick, k < n_gpr) for k in range(n_c)]
    with Pool(NPROC) as pool:  # GPR jobs first, one job per dispatch
        res = pool.map(_ml_worker, jobs, chunksize=1)
    M = pd.DataFrame([r for rows, _ in res for r in rows])
    PS = pd.concat([p for _, ps in res for p in ps], ignore_index=True)
    M.to_csv(RES / "ml_metrics.csv", index=False)
    PS.to_csv(RES / "ml_per_study_LOSO_RF.csv", index=False)
    with Pool(NPROC) as pool:
        T = pd.DataFrame([x for c in pool.map(_tune_worker, [(k, quick) for k in range(1 if quick else 3)]) for x in c])
    T.to_csv(RES / "ml_tuning.csv", index=False)
    lc_n = (12, 24) if quick else (12, 24, 48, 96)
    with Pool(NPROC) as pool:
        LC = pd.DataFrame(pool.map(_lc_worker, [(n, r) for n in lc_n for r in range(2 if quick else 4)]))
    LC.to_csv(RES / "ml_learning_curve.csv", index=False)
    summ = {"n_corpora": n_c, "n_corpora_gpr": n_gpr}
    agg = M.groupby(["model", "features", "scheme"]).agg(R2_mean=("R2", "mean"), R2_sd=("R2", "std"),
                                                          R2_q025=("R2", lambda x: x.quantile(0.025)),
                                                          R2_q975=("R2", lambda x: x.quantile(0.975)),
                                                          RMSE_mean=("RMSE", "mean"), RMSE_sd=("RMSE", "std"),
                                                          R2_cb_lo=("R2_lo", "mean"), R2_cb_hi=("R2_hi", "mean"),
                                                          n=("R2", "size")).reset_index()
    agg.to_csv(RES / "T_ml_summary.csv", index=False)
    summ["summary"] = agg.round(4).to_dict("records")
    piv = M[M.scheme == "LOSO"].pivot_table(index="corpus", columns=["model", "features"], values="R2")
    if ("RF", "chem") in piv and ("XGB", "chem") in piv:
        dd = (piv[("RF", "chem")] - piv[("XGB", "chem")]).dropna()
        summ["RF_minus_XGB_LOSO_R2_chem"] = {"mean": float(dd.mean()), "sd": float(dd.std()),
                                             "frac_RF_better": float((dd > 0).mean()), "n": len(dd),
                                             "sign_test_p": float(stats.binomtest(int((dd > 0).sum()), len(dd)).pvalue)}
    for mname in ("RF", "XGB", "GPR"):
        if (mname, "chem") in piv and (mname, "base") in piv:
            dd = (piv[(mname, "chem")] - piv[(mname, "base")]).dropna()
            summ[f"{mname}_chem_minus_base_LOSO_R2"] = {"mean": float(dd.mean()), "sd": float(dd.std()),
                                                        "frac_positive": float((dd > 0).mean()), "n": len(dd)}
    rf = M[(M.model == "RF") & (M.features == "chem") & (M.scheme == "LOSO")]
    summ["ICC_RF_chem_LOSO"] = {"mean": float(rf.icc_ICC.mean()), "sd": float(rf.icc_ICC.std())}
    summ["retransform_bias_mg_g"] = {"naive_mean": float(rf.bias_naive.mean()), "smear_mean": float(rf.bias_smear.mean())}
    summ["n_censored_mean"] = float(M[M.scheme == "LOSO"].n_censored.mean())
    summ["per_study_RMSE_RF_chem"] = PS[PS.features == "chem"].RMSE.quantile([0.1, 0.5, 0.9]).round(3).tolist()
    summ["per_study_RMSE_by_family"] = PS[PS.features == "chem"].groupby("family").RMSE.median().round(3).to_dict()
    tt = T.pivot_table(index="corpus", columns=["model", "tuned"], values="R2")
    summ["tuning_delta_R2"] = {m: float((tt[(m, True)] - tt[(m, False)]).mean()) for m in ("RF", "XGB")}
    summ["learning_curve"] = LC.groupby("n_studies")[["R2_LOSO", "R2_random"]].agg(["mean", "std"]).round(4).to_dict()
    summ["learning_curve"] = {str(k): v for k, v in summ["learning_curve"].items()}
    dump("ml", summ)
    print("ml done")


# =============================================================================
# Stage: cluster-level conformal prediction
# =============================================================================
def _uq_worker(args):
    k, quick = args
    d, _ = corpus_features(SEEDS["ml"] + k, 12 if quick else 30)
    alphas = (0.5, 0.4, 0.3, 0.2, 0.1, 0.05)
    out = []
    for tgt in ("q_app_VI", "q_tot_aq", "Ce_VI", "Ce_total"):
        r = val.conformal_loso(d, target=tgt, model_name="RF", features=val.CHEM, K=10, alphas=alphas, beta=0.9)
        r["target"], r["corpus"] = tgt, k
        out.append(r)
    return pd.concat(out)


def stage_uq(quick=False):
    n_c = 2 if quick else 10
    with Pool(NPROC) as pool:
        U = pd.concat(pool.map(_uq_worker, [(k, quick) for k in range(n_c)]), ignore_index=True)
    U.to_csv(RES / "uq_per_study.csv", index=False)
    rows = []
    for (tgt, meth, a), g in U.groupby(["target", "method", "alpha"]):
        nominal = 1 - a
        rows.append({"target": tgt, "method": meth, "nominal": nominal, "n_test_studies": len(g),
                     "mean_cov": g.coverage.mean(), "median_cov": g.coverage.median(),
                     "p10_cov": g.coverage.quantile(0.1), "worst_cov": g.coverage.min(),
                     "frac_studies_cov_ge_nominal": (g.coverage >= nominal - 1e-9).mean(),
                     "frac_studies_cov_ge_0.9": (g.coverage >= 0.9 - 1e-9).mean(),
                     "obs_weighted_cov": np.average(g.coverage, weights=g.n),
                     "mean_width_log10": g.width_log10.mean(), "n_cal_studies": g.n_cal_studies.mean()})
    S = pd.DataFrame(rows)
    S.to_csv(RES / "T_uq_summary.csv", index=False)
    dump("uq", {"n_corpora": n_c, "summary": S.round(4).to_dict("records")})
    print("uq done")


# =============================================================================
# Stage: treatment design with two endpoints and joint constraints
# =============================================================================
TARGETS = {"S1": (20.0, 0.5, 2.0), "S2": (5.0, 0.1, 0.5)}
RULES = {"nominal": None, "marginal90": 0.10, "joint90_bonf": 0.05, "joint90_adapt_k4": 0.05,
         "joint80_adapt_k4": 0.10, "nominal_verify": None}
ADAPT_DOSES = {4: (0.5, 1.0, 2.0, 4.0)}  # confirmatory batches at pH 3.0 and the target C0
MAX_ROUNDS = 4  # verification rounds for the model-guided, test-confirmed protocol
GRID_PH = np.round(np.arange(2.0, 7.01, 0.5), 2)
GRID_D = np.array([0.5, 1, 2, 3, 4, 6, 8], float)  # inside the training dose range (dose series)


def _design_frame(m, C0, I=0.01, T=298.15):
    G = pd.DataFrame([{"C0": C0, "pH": p, "D": D, "t": 1440.0, "T": T, "I": I, "SBET": m.SBET, "Vp": m.Vp,
                       "dp": m.dp, "pHpzc": m.pHpzc} for p in GRID_PH for D in GRID_D])
    G = ch.speciation_features(G, conc_col="C0")
    G["pH_minus_pzc"] = G["pH"] - G["pHpzc"]
    G["I_imputed"] = 0.0
    return G


def _choose(G, cps, C0, lim6, limt, rule, adapt=None):
    """Pick the lowest-dose design whose (point or upper-bound) effluent meets both limits.

    The decision quantities are modelled directly: log10 residual Cr(VI) and log10 residual
    total Cr, each with its own study-level CV+ interval (M1).
    """
    c6 = val.from_target(cps["Ce_VI"].full.predict(G[val.CHEM]), target="Ce_VI")
    ct = val.from_target(cps["Ce_total"].full.predict(G[val.CHEM]), target="Ce_total")
    G = G.assign(Ce6_pt=c6, Cet_pt=np.maximum(ct, c6))
    if RULES[rule] is None:
        G["Ce6_ub"], G["Cet_ub"] = G.Ce6_pt, G.Cet_pt
    elif adapt is not None:
        # study-adapted: centre shifted by the offset estimated from confirmatory batches
        for tgt, col in (("Ce_VI", "Ce6"), ("Ce_total", "Cet")):
            off, hw = adapt[tgt]
            c = cps[tgt].full.predict(G[val.CHEM]) + off
            G[f"{col}_pt"] = val.from_target(c, target=tgt)
            G[f"{col}_ub"] = val.from_target(c + hw, target=tgt)
        G["Cet_ub"] = np.maximum(G.Cet_ub, G.Ce6_ub)
    else:
        a = RULES[rule]
        _, hi6, _ = cps["Ce_VI"].interval(G, a, "M1")
        _, hit, _ = cps["Ce_total"].interval(G, a, "M1")
        G["Ce6_ub"] = val.from_target(hi6, target="Ce_VI")
        G["Cet_ub"] = val.from_target(hit, target="Ce_total")
    ok = G[(G.Ce6_ub <= lim6) & (G.Cet_ub <= limt)]
    if ok.empty:
        return None, G
    ok = ok.assign(margin=np.maximum(ok.Ce6_ub / lim6, ok.Cet_ub / limt))
    return ok.sort_values(["D", "margin"]).iloc[0], G


def _verify(m, G, cps, C0, lim6, limt, rng, max_rounds=MAX_ROUNDS, n_rep=2):
    """Model-guided, test-confirmed design.

    Candidates are ranked by predicted effluent (point models) at the model-preferred pH, lowest
    feasible dose first. Each round runs n_rep confirmatory batches at the candidate (true
    response + 3% analytical noise + study bias); the design is accepted when the mean measured
    Cr(VI) and total Cr are both below 90% of their limits. Otherwise the dose is escalated to the
    next grid value. Returns (accepted design row or None, rounds used).
    """
    best, Gp = _choose(G, cps, C0, lim6, limt, "nominal")
    if best is None:  # start from the dose the model ranks best at the most favourable pH
        Gp = Gp.assign(score=np.maximum(Gp.Ce6_pt / lim6, Gp.Cet_pt / limt))
        best = Gp.sort_values(["score", "D"]).iloc[0]
    pH0 = best.pH
    doses = [D_ for D_ in GRID_D if D_ >= best.D]
    bias = np.exp(rng.normal(0, 0.05))
    for r_, D_ in enumerate(doses[:max_rounds], start=1):
        c6, ct = _truth(m, C0, pH0, D_)
        m6 = np.mean([c6 * bias * (1 + rng.normal(0, 0.03)) for _ in range(n_rep)])
        mt = np.mean([ct * bias * (1 + rng.normal(0, 0.03)) for _ in range(n_rep)])
        if m6 <= 0.9 * lim6 and mt <= 0.9 * limt:
            row = Gp[(Gp.pH == pH0) & (Gp.D == D_)].iloc[0].copy()
            return row, r_
    return None, min(len(doses), max_rounds)


def _truth(m, C0, pH, D, I=0.01, T=298.15):
    tr = sim.run_batch(m, C0, pH, D, [1440.0], T, I)
    return float(tr["Ce_VI_true"][0]), float(tr["Ce_VI_true"][0] + tr["Cr3_aq_true"][0])


def _design_worker(args):
    k, quick, n_new = args
    seed = SEEDS["ml"] + k
    d, _ = corpus_features(seed, 12 if quick else 30)
    cps = {t: val.StudyCVPlus("RF", val.CHEM, K=10).fit(d, t) for t in ("Ce_VI", "Ce_total")}
    rng = np.random.default_rng(seed)
    hyper = sim.corpus_hyperpriors(rng)  # identical to the training corpus hyper-priors
    rng2 = np.random.default_rng(SEEDS["design"] + k)
    out, sens, pred = [], [], []
    fams = list(sim.FAMILIES)
    for j in range(n_new):
        m = sim.random_material(rng2, fams[j % len(fams)], 9000 + j, hyper, "E")
        for tname, (C0, l6, lt) in TARGETS.items():
            G = _design_frame(m, C0)
            # true attainability: lowest dose (any grid pH) that meets both limits in the simulator
            min_true_D = np.nan
            for Dg in GRID_D:
                if any(all(v <= lim for v, lim in zip(_truth(m, C0, pg, Dg), (l6, lt))) for pg in GRID_PH):
                    min_true_D = float(Dg)
                    break
            adapt_by_k = {}
            for kk, doses in ADAPT_DOSES.items():
                sub = G[(G.pH == 3.0) & G.D.isin(doses)]
                bias = np.exp(rng2.normal(0, 0.05))
                obs6, obst = [], []
                for _, rr in sub.iterrows():
                    c6, ct = _truth(m, C0, rr.pH, rr.D)
                    obs6.append(max(c6 * bias * (1 + rng2.normal(0, 0.03)), 0.005))
                    obst.append(max(ct * bias * (1 + rng2.normal(0, 0.03)), 0.01))
                offs = {tgt: float(np.mean(val.to_target(np.array(obs), tgt) - cps[tgt].full.predict(sub[val.CHEM])))
                        for tgt, obs in (("Ce_VI", obs6), ("Ce_total", obst))}
                for a_ in (0.05, 0.10):  # per-endpoint alpha of the joint-90 and joint-80 rules
                    adapt_by_k[(kk, a_)] = {tgt: (offs[tgt], cps[tgt].adapted_halfwidth(kk, a_, seed=k)) for tgt in offs}
            for rule in RULES:
                ad = adapt_by_k[(int(rule[-1]), RULES[rule])] if "_adapt_k" in rule else None
                if rule == "nominal_verify":
                    best, rounds = _verify(m, G, cps, C0, l6, lt, rng2)
                else:
                    best, Gp = _choose(G, cps, C0, l6, lt, rule, ad)
                    rounds = 0
                rec = {"cohort": k, "material": m.material_id, "family": m.family, "target": tname, "rule": rule,
                       "feasible": best is not None, "attainable": bool(np.isfinite(min_true_D)),
                       "min_true_D": min_true_D, "rounds": rounds}
                if best is not None:
                    c6, ct = _truth(m, C0, best.pH, best.D)
                    rec.update(pH=best.pH, D=best.D, Ce6_pred=best.Ce6_pt, Cet_pred=best.Cet_pt,
                               Ce6_ub=best.Ce6_ub, Cet_ub=best.Cet_ub, Ce6_true=c6, Cet_true=ct,
                               met_VI=c6 <= l6, met_tot=ct <= lt, met_both=(c6 <= l6) and (ct <= lt),
                               bound_holds=(c6 <= best.Ce6_ub) and (ct <= best.Cet_ub))
                    if rule == "nominal_verify":
                        for dp_, fD in ((-0.5, 1), (-0.25, 1), (0.25, 1), (0.5, 1), (0, 0.8), (0, 1.2)):
                            c6p, ctp = _truth(m, C0, best.pH + dp_, best.D * fD)
                            sens.append({"cohort": k, "material": m.material_id, "target": tname, "dpH": dp_,
                                         "fD": fD, "met_both": (c6p <= l6) and (ctp <= lt)})
                out.append(rec)
            # predicted vs true over a sub-grid (for calibration of predicted effluent)
            if tname == "S1" and j < 5:
                sub = G.iloc[::5]
                p6 = val.from_target(cps["Ce_VI"].full.predict(sub[val.CHEM]), target="Ce_VI")
                pt = val.from_target(cps["Ce_total"].full.predict(sub[val.CHEM]), target="Ce_total")
                _, h6, _ = cps["Ce_VI"].interval(sub, 0.05, "M1")
                _, ht, _ = cps["Ce_total"].interval(sub, 0.05, "M1")
                for (_, rr), a6, at, u6, ut in zip(sub.iterrows(), p6, pt, h6, ht):
                    c6, ct = _truth(m, C0, rr.pH, rr.D)
                    pred.append({"cohort": k, "material": m.material_id, "pH": rr.pH, "D": rr.D,
                                 "Ce6_pred": a6, "Cet_pred": max(at, a6), "Ce6_true": c6, "Cet_true": ct,
                                 "Ce6_ub95": val.from_target(u6, target="Ce_VI"),
                                 "Cet_ub95": val.from_target(ut, target="Ce_total")})
    return out, sens, pred


def stage_design(quick=False):
    n_c, n_new = (2, 4) if quick else (20, 15)
    with Pool(NPROC) as pool:
        res = pool.map(_design_worker, [(k, quick, n_new) for k in range(n_c)])
    O = pd.DataFrame([r for o, _, _ in res for r in o])
    SN = pd.DataFrame([r for _, s, _ in res for r in s])
    PR = pd.DataFrame([r for _, _, p in res for r in p])
    O.to_csv(RES / "design_outcomes.csv", index=False)
    SN.to_csv(RES / "design_sensitivity.csv", index=False)
    PR.to_csv(RES / "design_pred_vs_true.csv", index=False)
    summ = {"n_cohorts": n_c, "n_materials_per_cohort": n_new}
    for c in ("met_both", "met_VI", "met_tot", "bound_holds"):
        if c in O:
            O[c] = O[c].fillna(False).astype(bool)
    for (t, rule), g in O.groupby(["target", "rule"]):
        f = g[g.feasible]
        k, n = int(f.met_both.sum()), len(f)
        per = f.groupby("cohort").met_both.mean()
        N = len(g)
        unsafe = int((g.feasible & ~g.met_both).sum())
        att = g[g.attainable]
        summ[f"{t}|{rule}"] = {"n_materials": N, "n_attainable": int(g.attainable.sum()),
                               "unsafe_recommendations": unsafe, "unsafe_rate": unsafe / N, "unsafe_wilson95": wilson(unsafe, N),
                               "compliant_found": k, "found_of_attainable": float((att.feasible & att.met_both).mean()) if len(att) else np.nan,
                               "declared_infeasible": int((~g.feasible).sum()),
                               "declared_infeasible_but_attainable": int((~g.feasible & g.attainable).sum()),
                               "dose_overhead_median": float((f[f.met_both].D / f[f.met_both].min_true_D).median()) if k else np.nan,
                               "n_feasible": n, "success": k / n if n else np.nan,
                               "wilson95": wilson(k, n), "cohort_mean": float(per.mean()), "cohort_sd": float(per.std()),
                               "median_dose": float(f.D.median()) if n else np.nan,
                               "fail_VI": int((~f.met_VI.astype(bool)).sum()), "fail_tot": int((~f.met_tot.astype(bool)).sum()),
                               "joint_bound_holds": float(f.bound_holds.mean()) if n else np.nan,
                               "median_pH": float(f.pH.median()) if n else np.nan,
                               "mean_rounds": float(g.rounds.mean()) if "rounds" in g else np.nan,
                               "halfwidth_note": "total-Cr study-level 95% bound needs >= 19 calibration studies"}
    if len(SN):
        summ["sensitivity_joint"] = SN.groupby(["target", "dpH", "fD"]).met_both.mean().round(3).reset_index().to_dict("records")
    if len(PR):
        for e in ("Ce6", "Cet"):
            lt, lp = np.log10(PR[f"{e}_true"] + 1e-3), np.log10(PR[f"{e}_pred"] + 1e-3)
            summ[f"pred_vs_true_{e}_log10"] = {"spearman": float(stats.spearmanr(lt, lp).statistic),
                                               "median_abs_err_log10": float(np.median(np.abs(lt - lp))),
                                               "upper95_covers_truth": float((PR[f"{e}_true"] <= PR[f"{e}_ub95"]).mean())}
        both = (PR.Ce6_true <= PR.Ce6_ub95) & (PR.Cet_true <= PR.Cet_ub95)
        summ["pred_vs_true_joint_upper_bounds_hold"] = float(both.mean())
    dump("design", summ)
    print("design done")


# =============================================================================
# Stage: propagation of unreported chemistry through the workflow
# =============================================================================
def _missing_worker(args):
    k, quick, n_draw, n_draw_design = args
    seed = SEEDS["ml"] + k
    d, _ = corpus_features(seed, 12 if quick else 30)
    rng = np.random.default_rng(42 + k)
    rows = []
    # prediction level: outer grouped 5-fold; studies that did not report I and/or final pH
    for tr, te in val.splitter("grouped5", d):
        cp = val.StudyCVPlus("RF", val.CHEM, K=10).fit(d.iloc[tr], "q_app_VI")
        dte = d.iloc[te]
        for sid, g in dte.groupby("study_id"):
            I_missing = g["I"].isna().all()
            pH_unknown = g["pH_final"].isna().all()
            if not (I_missing or pH_unknown):
                continue
            g = g[g.series == "pH_edge"]
            base_lo, base_hi, base_c = cp.interval(g, 0.1, "M1")
            n = len(g)
            Iv = 10 ** rng.uniform(-3, -1, (n_draw, n)) if I_missing else np.full((n_draw, n), 0.01)
            dpH = rng.uniform(-0.5, 0.5, (n_draw, 1)) if pH_unknown else np.zeros((n_draw, 1))
            big = pd.concat([g] * n_draw, ignore_index=True)
            big["I"] = Iv.ravel()
            big["pH"] = np.clip(big["pH"].values + np.repeat(dpH, n, 1).ravel(), 1.5, 10)
            big = ch.speciation_features(big.drop(columns=[c for c in ("f_HCrO4", "f_CrO4", "f_Cr2O7", "f_H2CrO4", "z_mean_CrVI", "log_I") if c in big]), conc_col="C0")
            big["pH_minus_pzc"] = big["pH"] - big["pHpzc"]
            big["I_imputed"] = 0.0
            lo, hi, c = cp.interval(big, 0.1, "M1")
            c, lo, hi = c.reshape(n_draw, n), lo.reshape(n_draw, n), hi.reshape(n_draw, n)
            fH = big["f_HCrO4"].values.reshape(n_draw, n)
            for i in range(n):
                rows.append({"corpus": k, "study_id": sid, "pH": g["pH"].values[i], "I_missing": I_missing,
                             "pH_unknown": pH_unknown,
                             "fH_p05": np.quantile(fH[:, i], 0.05), "fH_p95": np.quantile(fH[:, i], 0.95),
                             "pred_spread90_log10": np.quantile(c[:, i], 0.95) - np.quantile(c[:, i], 0.05),
                             "width_nominal_log10": base_hi[i] - base_lo[i],
                             "width_union_log10": np.quantile(hi[:, i], 0.95) - np.quantile(lo[:, i], 0.05)})
    # design level: designer does not know I (true I = 0.01) -> dose and outcome distributions
    cps = {t: val.StudyCVPlus("RF", val.CHEM, K=10).fit(d, t) for t in ("Ce_VI", "Ce_total")}
    hyper = sim.corpus_hyperpriors(np.random.default_rng(seed))
    rng2 = np.random.default_rng(SEEDS["design"] + k)
    drows = []
    fams = list(sim.FAMILIES)
    for j in range(3 if quick else 6):
        m = sim.random_material(rng2, fams[j % len(fams)], 9000 + j, hyper, "E")
        C0, l6, lt = TARGETS["S1"]
        cache = {}
        for r in range(n_draw_design):
            I = float(10 ** rng.uniform(-3, -1))
            G = _design_frame(m, C0, I=I)          # designer's features use the assumed I
            best, rounds = _verify(m, G, cps, C0, l6, lt, np.random.default_rng(1000 * j + r))
            if best is None:
                drows.append({"corpus": k, "material": m.material_id, "draw": r, "I_assumed": I, "feasible": False})
                continue
            key = (best.pH, best.D)
            if key not in cache:
                cache[key] = _truth(m, C0, best.pH, best.D)
            c6, ct = cache[key]
            drows.append({"corpus": k, "material": m.material_id, "draw": r, "I_assumed": I, "feasible": True,
                          "pH": best.pH, "D": best.D, "Ce6_true": c6, "Cet_true": ct, "met_both": c6 <= l6 and ct <= lt,
                          "rounds": rounds})
    return pd.DataFrame(rows), pd.DataFrame(drows)


def stage_missing(quick=False):
    n_c = 1 if quick else 4
    with Pool(NPROC) as pool:
        res = pool.map(_missing_worker, [(k, quick, 200 if quick else 1000, 30 if quick else 300) for k in range(n_c)])
    P = pd.concat([a for a, _ in res], ignore_index=True)
    Dd = pd.concat([b for _, b in res], ignore_index=True)
    P.to_csv(RES / "missing_prediction.csv", index=False)
    Dd.to_csv(RES / "missing_design.csv", index=False)
    P["pH_bin"] = pd.cut(P.pH, [1.5, 3.5, 5.0, 6.5, 8.5])
    summ = {"n_rows": len(P),
            "by_pH_bin": P.groupby("pH_bin", observed=True)[["pred_spread90_log10", "width_nominal_log10", "width_union_log10"]].median().round(3).reset_index().astype(str).to_dict("records"),
            "fH_width90_by_pH_bin": (P.fH_p95 - P.fH_p05).groupby(P.pH_bin, observed=True).median().round(3).astype(float).to_dict()}
    summ["fH_width90_by_pH_bin"] = {str(k): v for k, v in summ["fH_width90_by_pH_bin"].items()}
    f = Dd[Dd.feasible.astype(bool)] if "D" in Dd else Dd.iloc[:0]
    if f.empty:
        summ["design"] = {"feasible_frac": float(Dd.feasible.mean()) if len(Dd) else None}
        dump("missing", summ)
        print("missing done (no feasible designs)")
        return
    summ["design"] = {"n_materials": int(Dd.groupby(["corpus", "material"]).ngroups),
                      "draws_per_material": int(Dd.groupby(["corpus", "material"]).size().median()),
                      "feasible_frac": float(Dd.feasible.mean()),
                      "dose_range_per_material": f.groupby(["corpus", "material"])["D"].agg(["min", "median", "max"]).reset_index().to_dict("records"),
                      "compliance_prob_per_material": f.groupby(["corpus", "material"]).met_both.mean().round(3).tolist(),
                      "Cet_true_p95_per_material": f.groupby(["corpus", "material"]).Cet_true.quantile(0.95).round(3).tolist()}
    dump("missing", summ)
    print("missing done")


# =============================================================================
# Figures
# =============================================================================
def figures():
    chem = json.loads((RES / "chem.json").read_text())
    # Fig 1: chemistry
    fig, ax = plt.subplots(2, 2, figsize=(7.4, 5.8), constrained_layout=True)
    lss = {"minteq.v4.dat": "-", "minteq.dat": ":", "llnl.dat": "--"}
    for a, c0, L in ((ax[0, 0], 10, "a"), (ax[0, 1], 500, "b")):
        t = pd.read_csv(RES / f"speciation_C0_{c0}.csv")
        for k, (sp, lab) in enumerate([("f_HCrO4-", "HCrO$_4^-$"), ("f_CrO4-2", "CrO$_4^{2-}$"), ("f_Cr2O7-2", "Cr$_2$O$_7^{2-}$ (per Cr)")]):
            for db, ls in lss.items():
                s = t[t.database == db]
                a.plot(s.pH, s[sp], ls, color=C[k], lw=1.5, label=lab if db == "minteq.v4.dat" else None)
        a.set(xlabel="pH", ylabel="Fraction of Cr(VI)", title=f"Cr(VI) {c0} mg L$^{{-1}}$, I = 0.01 M")
        panel(a, L)
    ax[0, 0].legend(loc="center left")
    ax[0, 1].text(0.02, 0.97, "solid minteq.v4\ndotted minteq\ndashed llnl", transform=ax[0, 1].transAxes, va="top", fontsize=8, color=MUTED)
    rb = pd.read_csv(RES / "T_redox_boundary.csv")
    for db, ls in lss.items():
        s = rb[rb.database == db]
        ax[1, 0].plot(s.pH, s.Eh50_V, ls, color=C[0], label=db.replace(".dat", ""))
    pw = np.array([2, 10])
    ax[1, 0].plot(pw, (20.75 - pw) * 0.05916, color=MUTED, lw=1, ls="-.", label="O$_2$/H$_2$O")
    ax[1, 0].set(xlabel="pH", ylabel="Eh at Cr(VI) = Cr(III) (V)", title="Redox boundary, 50 mg L$^{-1}$ Cr")
    ax[1, 0].legend()
    panel(ax[1, 0], "c")
    so = pd.read_csv(RES / "T_cr3_solubility_SI.csv")
    for db, ls in (("minteq.v4.dat", "-"), ("minteq.dat", ":")):
        s = so[(so.database == db) & (so.I == 0.01)]
        b = s[s.status == "bracketed"]
        ax[1, 1].plot(b.pH, b.S_mgL, ls, color=C[2], label=db.replace(".dat", ""), marker="o", ms=3)
        u = s[s.status == "undersaturated_to_cmax"]
        ax[1, 1].scatter(u.pH, np.full(len(u), 5000), marker="^", color=C[2], s=22, facecolors="white" if db == "minteq.dat" else C[2])
    ax[1, 1].set(yscale="log", xlabel="pH", ylabel="Dissolved Cr(III) at SI = 0 (mg L$^{-1}$)",
                 title="Cr(OH)$_3$(am) solubility (I = 0.01 M)")
    ax[1, 1].text(3.05, 2500, "▲ undersaturated up to 5000 mg L$^{-1}$", fontsize=7.5, color=MUTED, va="top")
    ax[1, 1].legend(loc="lower left")
    panel(ax[1, 1], "d")
    fig.savefig(FIG / "fig1_chemistry.png"); plt.close(fig)

    # Fig 2: benchmark across mechanisms (synthetic)
    B = pd.read_csv(RES / "bench_studies.csv.gz")
    fig, ax = plt.subplots(2, 2, figsize=(7.4, 5.8), constrained_layout=True)
    data = [B[B.scenario == s].share_red.values for s in sim.SCENARIOS]
    bp = ax[0, 0].boxplot(data, tick_labels=list(sim.SCENARIOS), showfliers=False, patch_artist=True, widths=0.55,
                          medianprops={"color": INK})
    for patch, s in zip(bp["boxes"], sim.SCENARIOS):
        patch.set_facecolor(SC_COL[s]); patch.set_alpha(0.55); patch.set_edgecolor(SC_COL[s])
    ax[0, 0].set(xlabel="Synthetic scenario", ylabel="True reduction share", title="Mechanism mix (synthetic)")
    panel(ax[0, 0], "a")
    bins = np.array([0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0])
    mids = 0.5 * (bins[1:] + bins[:-1])
    for s in ("B", "C", "D", "E"):
        g = B[B.scenario == s]
        cut = pd.cut(g.share_red, bins)
        fr = g.groupby(cut, observed=False).R2_best.apply(lambda x: (x >= 0.95).mean() if len(x) >= 20 else np.nan)
        ax[0, 1].plot(mids, fr.values, "-o", ms=3.5, color=SC_COL[s], label=s)
        pub = g[(g.R2_best >= 0.95) & (g.qmax_L < 0.99e5)]
        cutp = pd.cut(pub.share_red, bins)
        med = pub.groupby(cutp, observed=False).qmax_ratio.median()
        q1 = pub.groupby(cutp, observed=False).qmax_ratio.quantile(0.25)
        q3 = pub.groupby(cutp, observed=False).qmax_ratio.quantile(0.75)
        n_ = pub.groupby(cutp, observed=False).size()
        ok = (n_ >= 20).values
        ax[1, 0].plot(mids[ok], med.values[ok], "-o", ms=3.5, color=SC_COL[s], label=s)
        ax[1, 0].fill_between(mids[ok], q1.values[ok], q3.values[ok], color=SC_COL[s], alpha=0.15, lw=0)
    gA = B[B.scenario == "A"]
    ax[1, 0].axhline(gA[gA.R2_best >= 0.95].qmax_ratio.median(), color=MUTED, ls=":", lw=1.2, label="A (control)")
    ax[0, 1].axhline((gA.R2_best >= 0.95).mean(), color=MUTED, ls=":", lw=1.2, label="A (control)")
    ax[0, 1].set(xlabel="True reduction share", ylabel="Share of fits with R$^2$ ≥ 0.95", ylim=(0, 1.02),
                 title="Fit quality vs mechanism")
    ax[0, 1].legend(title="Scenario", ncol=2)
    panel(ax[0, 1], "b")
    ax[1, 0].set(xlabel="True reduction share", ylabel="Fitted / true q$_{max}$", yscale="log", ylim=(0.7, 30),
                 title="q$_{max}$ bias, fits with R$^2$ ≥ 0.95")
    ax[1, 0].legend(ncol=2)
    panel(ax[1, 0], "c")
    E = B[B.scenario == "E"]
    parts = [("share_ads", "adsorbed\nCr(VI)"), ("share_cr3_surface", "Cr(III)\non solid"),
             ("share_cr3_aq", "dissolved\nCr(III)"), ("share_precip", "Cr(OH)$_3$\nprecip.")]
    ax[1, 1].boxplot([E[p].values for p, _ in parts], tick_labels=[l for _, l in parts], showfliers=False, widths=0.55,
                     medianprops={"color": INK})
    ax[1, 1].set(ylabel="Share of Cr(VI) loss", title="Fate of removed Cr(VI), scenario E")
    panel(ax[1, 1], "d")
    fig.savefig(FIG / "fig2_benchmark.png"); plt.close(fig)

    # Fig 3: identifiability
    F = pd.read_csv(RES / "ident_fim_detection.csv")
    Q = pd.read_csv(RES / "ident_equivalence.csv")
    P = pd.read_csv(RES / "ident_profiles.csv")
    fig, ax = plt.subplots(2, 2, figsize=(7.4, 5.8), constrained_layout=True)
    bins = np.array([0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0]); mids = 0.5 * (bins[1:] + bins[:-1])
    cut = pd.cut(F.share_red, bins)
    for i, dsn in enumerate(idf.DESIGNS):
        ax[0, 0].plot(mids, F.groupby(cut, observed=False)[f"post_sd_{dsn}"].median().values, "-o", ms=3.5,
                      color=C[i], label=f"{dsn}: {idf.DESIGN_LABEL[dsn]}")
    ax[0, 0].plot(mids, F.groupby(cut, observed=False)["prior_sd_D1"].median().values, ":", color=MUTED, label="prior")
    ax[0, 0].set(yscale="log", xlabel="True reduction share", ylabel="Posterior SD of reduction share",
                 title="Local (Fisher) identifiability, model known")
    ax[0, 0].legend(fontsize=6.5, loc="lower right")
    panel(ax[0, 0], "a")
    labs = {"lack_of_fit": "isotherm lack of fit (C)", "dissolved_cr3": "dissolved Cr(III) z-test (B)", "xps_cr3": "XPS Cr(III) > 0.15 (A)"}
    for i, (k, l) in enumerate(labs.items()):
        ax[0, 1].plot(mids, F.groupby(cut, observed=False)[k].mean().values, "-o", ms=3.5, color=C[i], label=l)
    ax[0, 1].set(xlabel="True reduction share", ylabel="Detection rate", ylim=(0, 1.02), title="Model-free diagnostics")
    ax[0, 1].legend(fontsize=7.5)
    panel(ax[0, 1], "b")
    qq = Q.copy()
    qq["grp"] = np.where(qq.arm == "control_A", "control (A)", np.where(qq.share_red > 0.3, "reduction > 30%", np.where(qq.share_red > 0.05, "reduction 5–30%", "<5%")))
    qq = qq[qq.grp != "<5%"]
    x = np.arange(len(idf.DESIGNS)); w = 0.13
    k = 0
    for sd, hatch in ((0.03, None), (0.10, "//")):
        for grp, col in (("control (A)", MUTED), ("reduction 5–30%", C[3]), ("reduction > 30%", C[1])):
            v = [qq[(qq.sd_rel == sd) & (qq.grp == grp) & (qq.design == d)].indistinguishable.mean() for d in idf.DESIGNS]
            ax[1, 0].bar(x + (k - 2.5) * w, v, w * 0.92, color=col, hatch=hatch, edgecolor="white" if hatch is None else col,
                         alpha=0.9 if hatch is None else 0.5, label=f"{grp}, {int(sd*100)}% noise")
            k += 1
    ax[1, 0].set_xticks(x, idf.DESIGNS)
    ax[1, 0].set(ylabel="Adsorption-only model not rejected", ylim=(0, 1.42), title="Equivalence test (p > 0.05)")
    ax[1, 0].set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax[1, 0].legend(fontsize=6.0, ncol=2, loc="upper center")
    panel(ax[1, 0], "c")
    for i, ((sid, dsn), g) in enumerate(P.groupby(["study_id", "design"])):
        col = C[["D1", "D3", "D4"].index(dsn)]
        hi_ = g.true_share.iloc[0] > 0.4
        g = g.sort_values("log10_R0")
        smin = g.loc[g.delta.idxmin(), "share_red"]
        ax[1, 1].plot(g.log10_R0 - g.true_log10_R0, g.delta, "-" if hi_ else "--", marker="o", ms=3, color=col,
                      label=f"{dsn}, true {g.true_share.iloc[0]:.2f}, at min {smin:.2f}")
    ax[1, 1].axhline(3.84, color=MUTED, lw=1, ls=":")
    ax[1, 1].axvline(0, color=MUTED, lw=0.8)
    ax[1, 1].set(xlabel="log$_{10}$(R$_0$ / true R$_0$)", ylabel="Δ(−2 log L)", ylim=(-1, 60),
                 title="Profile likelihood of reductant capacity")
    ax[1, 1].legend(fontsize=6.3, title="reduction share", title_fontsize=6.5)
    panel(ax[1, 1], "d")
    fig.savefig(FIG / "fig3_identifiability.png"); plt.close(fig)

    # Fig 4: ML validation
    M = pd.read_csv(RES / "ml_metrics.csv")
    LC = pd.read_csv(RES / "ml_learning_curve.csv")
    PS = pd.read_csv(RES / "ml_per_study_LOSO_RF.csv")
    fig, ax = plt.subplots(1, 3, figsize=(7.6, 3.0), constrained_layout=True, gridspec_kw={"width_ratios": [2.2, 1, 1.2]})
    schemes = ["random", "grouped5", "LOSO", "LOFO", "external"]
    combos = [(m, f) for m in ("RF", "XGB", "GPR") for f in ("base", "chem")]
    w = 0.13
    for i, (m, f) in enumerate(combos):
        g = M[(M.model == m) & (M.features == f)]
        mu = [g[g.scheme == s].R2.mean() if (g.scheme == s).any() else np.nan for s in schemes]
        lo = [g[g.scheme == s].R2.quantile(0.025) if (g.scheme == s).any() else np.nan for s in schemes]
        hi = [g[g.scheme == s].R2.quantile(0.975) if (g.scheme == s).any() else np.nan for s in schemes]
        xx = np.arange(len(schemes)) + (i - 2.5) * w
        col = C[("RF", "XGB", "GPR").index(m)]
        ax[0].bar(xx, mu, w * 0.92, color=col, alpha=1 if f == "chem" else 0.45, label=f"{m} {'+spec.' if f == 'chem' else 'base'}")
        ax[0].errorbar(xx, mu, yerr=[np.array(mu) - np.array(lo), np.array(hi) - np.array(mu)], fmt="none", ecolor=INK, lw=0.7, capsize=1.5)
    ax[0].set_xticks(range(len(schemes)), ["Random\n5-fold", "Grouped\n5-fold", "LOSO", "LOFO", "External\ncorpus"])
    ax[0].set(ylabel="R$^2$ (apparent uptake)", ylim=(min(-0.1, np.nanmin(M.R2) - 0.05), 1))
    ax[0].legend(fontsize=6.5, ncol=3, loc="lower left")
    panel(ax[0], "a")
    g = LC.groupby("n_studies")
    ax[1].errorbar(g.R2_LOSO.mean().index, g.R2_LOSO.mean(), yerr=g.R2_LOSO.std(), fmt="-o", ms=3.5, color=C[0], capsize=2, label="LOSO")
    ax[1].errorbar(g.R2_random.mean().index, g.R2_random.mean(), yerr=g.R2_random.std(), fmt="--s", ms=3.5, color=MUTED, capsize=2, label="random")
    ax[1].set(xscale="log", xlabel="Number of studies", ylabel="R$^2$ (RF, +spec.)")
    ax[1].xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax[1].set_xticks([12, 24, 48, 96], ["12", "24", "48", "96"])
    ax[1].legend()
    panel(ax[1], "b")
    pc = PS[PS.features == "chem"]
    fams = list(sim.FAMILIES)
    abbr = {"AC_chemical": "ACc", "AC_physical": "ACp", "biochar_wood": "BCw", "biochar_manure": "BCm",
            "hydrochar": "HC", "porous_carbon": "PC"}
    ax[2].boxplot([pc[pc.family == f].RMSE.values for f in fams], tick_labels=[abbr[f] for f in fams], showfliers=False,
                  widths=0.55, medianprops={"color": INK})
    ax[2].tick_params(axis="x", labelsize=7, rotation=45)
    ax[2].set(ylabel="Per-study RMSE, LOSO (mg g$^{-1}$)")
    panel(ax[2], "c")
    fig.savefig(FIG / "fig4_ml_validation.png"); plt.close(fig)

    # Fig 5: conformal
    U = pd.read_csv(RES / "uq_per_study.csv")
    fig, ax = plt.subplots(1, 3, figsize=(7.6, 2.9), constrained_layout=True)
    names = {"M0": "pooled residual (M0)", "M1": "study CV+, weighted (M1)", "M2": "study CV+, study score (M2)"}
    u = U[U.target == "q_app_VI"]
    for i, (meth, lab) in enumerate(names.items()):
        g = u[u.method == meth].groupby("alpha").coverage
        nom = 1 - g.mean().index.values
        ax[0].plot(nom, g.mean().values, "-o", ms=3.5, color=C[i], label=lab)
        ax[0].plot(nom, g.quantile(0.1).values, ":", color=C[i], lw=1)
    ax[0].plot([0.45, 1], [0.45, 1], color=MUTED, lw=1)
    ax[0].set(xlabel="Nominal level", ylabel="Coverage (mean; dotted: 10th pct study)", title="Calibration, unseen studies")
    ax[0].legend(fontsize=6.5)
    panel(ax[0], "a")
    g90 = u[u.alpha == 0.1]
    ax[1].boxplot([g90[g90.method == m].coverage.values for m in names], tick_labels=["M0", "M1", "M2"], widths=0.55, showfliers=True,
                  medianprops={"color": INK},
                  flierprops={"markersize": 2.5})
    ax[1].axhline(0.9, color=MUTED, ls=":", lw=1)
    ax[1].set(ylabel="Per-study coverage at 90%", title="Spread across studies")
    panel(ax[1], "b")
    tl = {"q_app_VI": "q$_{app,VI}$", "q_tot_aq": "q$_{tot,aq}$", "Ce_VI": "C$_{e,VI}$", "Ce_total": "C$_{e,tot}$"}
    for i, tgt in enumerate(tl):
        gg = U[(U.target == tgt) & (U.alpha == 0.1)].groupby("method").coverage.mean()
        ax[2].bar(np.arange(3) + (i - 1.5) * 0.2, gg.reindex(list(names)).values, 0.19, color=C[i], label=tl[tgt])
    ax[2].set_xticks(range(3), ["M0", "M1", "M2"])
    ax[2].axhline(0.9, color=MUTED, ls=":", lw=1)
    ax[2].set(ylabel="Mean per-study coverage (90%)", ylim=(0.5, 1.02), title="By endpoint")
    ax[2].legend(fontsize=6.5, ncol=2, loc="lower right")
    panel(ax[2], "c")
    fig.savefig(FIG / "fig5_conformal.png"); plt.close(fig)

    # Fig 6: design
    O = pd.read_csv(RES / "design_outcomes.csv")
    SN = pd.read_csv(RES / "design_sensitivity.csv") if (RES / "design_sensitivity.csv").stat().st_size > 2 else pd.DataFrame()
    for c in ("met_both",):
        O[c] = O[c].fillna(False).astype(bool)
    fig, ax = plt.subplots(1, 3, figsize=(7.6, 3.1), constrained_layout=True, gridspec_kw={"width_ratios": [1.9, 1, 1]})
    rules = ["nominal", "marginal90", "joint90_bonf", "joint90_adapt_k4", "joint80_adapt_k4", "nominal_verify"]
    rl = {"nominal": "point", "marginal90": "marg.\n90%", "joint90_bonf": "joint\n90%", "joint90_adapt_k4": "joint 90%\n+4 b.",
          "joint80_adapt_k4": "joint 80%\n+4 b.", "nominal_verify": "point\n+verify"}
    x = np.arange(len(rules)); w = 0.2
    for j, (t, hatch) in enumerate((("S1", None), ("S2", "//"))):
        uns, found, ue, fe = [], [], [[], []], [[], []]
        for r in rules:
            g = O[(O.target == t) & (O.rule == r)]
            ku, n = int((g.feasible & ~g.met_both).sum()), len(g)
            att = g[g.attainable]
            kf, na = int((att.feasible & att.met_both).sum()), len(att)
            pu, pf = ku / n, kf / na
            lu, hu = wilson(ku, n); lf, hf = wilson(kf, na)
            uns.append(pu); found.append(pf)
            ue[0].append(max(pu - lu, 0)); ue[1].append(max(hu - pu, 0)); fe[0].append(max(pf - lf, 0)); fe[1].append(max(hf - pf, 0))
        off = -1.5 * w if j == 0 else 0.5 * w
        ax[0].bar(x + off, uns, w * 0.95, color=C[1], hatch=hatch, edgecolor="white" if hatch is None else C[1],
                  alpha=1 if hatch is None else 0.55, label=f"unsafe, {t}")
        ax[0].bar(x + off + w, found, w * 0.95, color=C[2], hatch=hatch, edgecolor="white" if hatch is None else C[2],
                  alpha=1 if hatch is None else 0.55, label=f"found, {t}")
        ax[0].errorbar(x + off, uns, yerr=ue, fmt="none", ecolor=INK, lw=0.7, capsize=1.5)
        ax[0].errorbar(x + off + w, found, yerr=fe, fmt="none", ecolor=INK, lw=0.7, capsize=1.5)
    ax[0].set_xticks(x, [rl[r] for r in rules], fontsize=7.5)
    ax[0].set(ylabel="Fraction of materials", ylim=(0, 1.0), title="Design rules (S1 solid, S2 hatched)")
    ax[0].legend(fontsize=6.5, ncol=2, loc="upper left")
    panel(ax[0], "a")
    f = O[(O.target == "S1") & O.feasible & O.rule.isin(["nominal", "nominal_verify"])]
    for i, r in enumerate(("nominal", "nominal_verify")):
        g = f[f.rule == r]
        ax[1].scatter(g.Cet_pred + 1e-3, g.Cet_true + 1e-3, s=10, color=[C[1], C[2]][i], alpha=0.6, label=rl[r].replace("\n", " "),
                      edgecolor="none")
    ax[1].plot([1e-3, 30], [1e-3, 30], color=MUTED, lw=1)
    ax[1].axhline(2.0, color=C[1], lw=0.8, ls=":")
    ax[1].set(xscale="log", yscale="log", xlabel="Predicted total Cr (mg L$^{-1}$)", ylabel="True total Cr (mg L$^{-1}$)",
              title="S1 designs")
    ax[1].legend(fontsize=6.5)
    panel(ax[1], "b")
    if len(SN):
        s1 = SN[SN.target == "S1"].groupby(["dpH", "fD"]).met_both.mean().reset_index()
        labels = [f"pH{'+' if a_ > 0 else ''}{a_:g}" if b_ == 1 else f"dose ×{b_:g}" for a_, b_ in zip(s1.dpH, s1.fD)]
        ax[2].bar(range(len(s1)), s1.met_both, color=C[0])
        ax[2].axhline(1.0, color=INK, ls=":", lw=1)
        ax[2].set_xticks(range(len(s1)), labels, rotation=40, fontsize=7.5)
        ax[2].set(ylabel="Compliance", ylim=(0, 1.05), title="Sensitivity")
    panel(ax[2], "c")
    fig.savefig(FIG / "fig6_design.png"); plt.close(fig)

    # Fig 7: missing chemistry
    MP = pd.read_csv(RES / "missing_prediction.csv")
    MD = pd.read_csv(RES / "missing_design.csv")
    fig, ax = plt.subplots(1, 3, figsize=(7.6, 2.9), constrained_layout=True)
    g = MP.groupby("pH")
    ax[0].plot(g.pH.mean(), (g.fH_p95.median() - g.fH_p05.median()), "-o", ms=3.5, color=C[0])
    ax[0].set(xlabel="Reported pH", ylabel="90% range of f(HCrO$_4^-$)", title="Speciation feature")
    panel(ax[0], "a")
    ax[1].plot(g.pH.mean(), g.width_nominal_log10.median(), "-o", ms=3.5, color=C[0], label="I and pH as reported")
    ax[1].plot(g.pH.mean(), g.width_union_log10.median(), "-s", ms=3.5, color=C[1], label="union over 1000 scenarios")
    ax[1].plot(g.pH.mean(), g.pred_spread90_log10.median(), "--", color=MUTED, label="spread of point prediction")
    ax[1].set(xlabel="Reported pH", ylabel="log$_{10}$ width", title="Prediction (M1, 90%)")
    ax[1].legend(fontsize=6.5)
    panel(ax[1], "b")
    f = MD[MD.feasible]
    mats = list(f.groupby(["corpus", "material"]).groups)
    ax[2].boxplot([f[(f.corpus == c) & (f.material == m)].D.values for c, m in mats], showfliers=False, widths=0.6,
                  medianprops={"color": INK})
    for i, (c, m) in enumerate(mats):
        fr = MD[(MD.corpus == c) & (MD.material == m)].feasible.mean()
        ax[2].text(i + 1, 8.9, f"{fr:.2f}", ha="center", fontsize=6.5, rotation=90, va="bottom", color=MUTED)
    ax[2].set_xticks(range(1, len(mats) + 1), [str(i + 1) for i in range(len(mats))], fontsize=7)
    ax[2].set(xlabel="New material", ylabel="Accepted dose (g L$^{-1}$)", ylim=(0, 10.5), title="Unknown I")
    panel(ax[2], "c")
    fig.savefig(FIG / "fig7_missing_chemistry.png"); plt.close(fig)
    print("figures done")


STAGES = {"chem": stage_chem, "bench": stage_bench, "ident": stage_ident, "profile": stage_profile, "ml": stage_ml, "uq": stage_uq,
          "design": stage_design, "missing": stage_missing}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all", choices=list(STAGES) + ["all", "figures"])
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    todo = list(STAGES) if a.stage == "all" else ([] if a.stage == "figures" else [a.stage])
    timing = {}
    for s in todo:
        t = time.time()
        STAGES[s](a.quick)
        timing[s] = round(time.time() - t, 1)
        print(f"[{s}] {timing[s]} s", flush=True)
    if a.stage in ("all", "figures"):
        figures()
    if todo:
        f = RES / "timing.json"
        old = json.loads(f.read_text()) if f.exists() else {}
        f.write_text(json.dumps(old | timing, indent=2))
    print(f"total {time.time() - t0:.0f} s")
