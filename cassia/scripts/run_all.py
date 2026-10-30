"""Reproduce every number and figure in the CASSIA manuscript.

    python scripts/run_all.py            # full run (~30 min on 4 cores)
    python scripts/run_all.py --quick    # reduced grid for smoke testing
"""
from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from cassia import chemistry as ch, dataset as ds, fitting as fit, simulate as sim, validation as val

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = ROOT / "manuscript" / "figures"
RES.mkdir(exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

# Reference categorical palette (validated, fixed order), ink and grid tokens.
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"]
INK, MUTED, GRID = "#0b0b0b", "#898781", "#e1e0d9"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#c3c2b7", "axes.labelcolor": INK,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True,
                     "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "lines.linewidth": 2, "savefig.dpi": 300,
                     "figure.facecolor": "white", "legend.frameon": False})
OUT = {}


def logk_table():
    rx = {"HCrO4-": r"CrO4-2 \+ H\+ = HCrO4-|H\+ \+ CrO4-2 = HCrO4-",
          "H2CrO4": r"CrO4-2 \+ 2 H\+ = H2CrO4|2 H\+ \+ CrO4-2 = H2CrO4",
          "Cr2O7-2": r"2 CrO4-2 \+ 2 H\+ = Cr2O7-2 \+ H2O|2 H\+ \+ 2 CrO4-2 = Cr2O7-2 \+ H2O"}
    rows = []
    for db in ch.DATABASES:
        lines = (ch.DB_DIR / db).read_text(errors="ignore").splitlines()
        for sp, pat in rx.items():
            for i, l in enumerate(lines):
                if re.fullmatch(r"\s*(" + pat + r")\s*", l):
                    for l2 in lines[i + 1:i + 4]:
                        m = re.search(r"log_?k\s+(-?[\d.]+)", l2)
                        if m:
                            rows.append({"database": db, "species": sp, "logK": float(m.group(1))})
                            break
                    break
    return pd.DataFrame(rows).pivot(index="species", columns="database", values="logK")


def exp_chemistry():
    t = logk_table()
    t.to_csv(RES / "T2_logK_comparison.csv")
    OUT["logK"] = t.round(3).to_dict()
    pH = np.round(np.arange(1.5, 10.01, 0.1), 2)
    tabs = {c0: ch.speciation_table(pH, cr_mgL=c0) for c0 in (10.0, 500.0)}
    for c0, tb in tabs.items():
        tb.to_csv(RES / f"speciation_C0_{int(c0)}.csv", index=False)
    # database discrepancy
    disc = {}
    for c0, tb in tabs.items():
        a, b = (tb[tb.database == d].set_index("pH") for d in ch.DATABASES)
        cols = [c for c in ("f_HCrO4-", "f_CrO4-2", "f_Cr2O7-2", "f_H2CrO4") if c in a and c in b]
        d = (a[cols] - b[cols]).abs()
        disc[int(c0)] = {c: {"max_abs_diff": float(d[c].max()), "at_pH": float(d[c].idxmax())} for c in cols}
    OUT["db_discrepancy"] = disc
    # pKa-like crossover HCrO4-/CrO4-2 at 10 mg/L
    tb = tabs[10.0]
    OUT["crossover_pH"] = {}
    for db in ch.DATABASES:
        s = tb[tb.database == db]
        diff = s["f_HCrO4-"] - s["f_CrO4-2"]
        OUT["crossover_pH"][db] = float(np.interp(0, -diff.values, s["pH"].values))
    # ionic-strength sensitivity of the crossover
    OUT["crossover_vs_I"] = {}
    for I in (0.001, 0.01, 0.1):
        p = np.round(np.arange(5.5, 7.01, 0.05), 2)
        fr = [ch.speciate_cr6(x, 10.0, ionic_strength=I).fractions for x in p]
        d = np.array([f.get("HCrO4-", 0) - f.get("CrO4-2", 0) for f in fr])
        OUT["crossover_vs_I"][I] = float(np.interp(0, -d, p))
    # Cr(III) solubility
    pS = np.round(np.arange(4.5, 10.01, 0.25), 2)
    S = pd.DataFrame({"pH": pS, "S_CrOH3am_mgL": [ch.cr3_solubility_mgL(x) for x in pS],
                      "S_Cr2O3_mgL": [ch.cr3_solubility_mgL(x, phase="Cr2O3") for x in pS]})
    S.to_csv(RES / "cr3_solubility.csv", index=False)
    OUT["cr3_solubility_mgL"] = S.set_index("pH").loc[[5.0, 6.0, 7.0, 8.0]].round(5).to_dict()
    # pe threshold for 50% Cr(VI) at 50 mg/L
    pe_grid = np.arange(-4, 20.01, 0.25)
    thr = []
    for db in ch.DATABASES:
        for p in np.arange(2, 10.01, 0.5):
            pes = pe_grid[pe_grid <= 20.75 - p]  # stay inside the water-stability field
            f6 = np.array([ch.redox_scenario(p, 50.0, pe, database=db)["frac_CrVI"] for pe in pes])
            pe50 = float(np.interp(0.5, f6, pes))
            thr.append({"database": db, "pH": p, "pe50": pe50, "Eh50_V": pe50 * 0.05916})
    thr = pd.DataFrame(thr)
    thr.to_csv(RES / "redox_threshold.csv", index=False)
    OUT["pe50"] = thr.pivot(index="pH", columns="database", values="pe50").loc[[2.0, 4.0, 6.0, 8.0]].round(2).to_dict()

    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5), constrained_layout=True)
    for j, c0 in enumerate((10.0, 500.0)):
        tb = tabs[c0]
        for k, (sp, lab) in enumerate([("f_HCrO4-", "HCrO$_4^-$"), ("f_CrO4-2", "CrO$_4^{2-}$"),
                                       ("f_Cr2O7-2", "Cr$_2$O$_7^{2-}$ (as Cr)")]):
            for db, ls in zip(ch.DATABASES, ("-", "--")):
                s = tb[tb.database == db]
                ax[j].plot(s["pH"], s[sp], ls, color=C[k], lw=1.6,
                           label=f"{lab}" if db == ch.DATABASES[0] else None)
        ax[j].set(title=f"Cr(VI) {int(c0)} mg/L, I = 0.01 M", xlabel="pH", ylabel="Fraction of Cr(VI)")
    ax[0].legend(fontsize=7, loc="center right")
    ax[1].text(0.02, 0.97, "solid: minteq.v4\ndashed: llnl", transform=ax[1].transAxes, ha="left", va="top",
               fontsize=7, color=MUTED)
    for db, ls in zip(ch.DATABASES, ("-", "--")):
        s = thr[thr.database == db]
        ax[2].plot(s["pH"], s["Eh50_V"], ls, color=C[0], label=db.replace(".dat", ""))
    pw = np.array([2, 10])
    ax[2].plot(pw, (20.75 - pw) * 0.05916, color=MUTED, lw=1, ls=":")
    ax[2].text(9.8, (20.75 - 9.8) * 0.05916 + 0.03, "O$_2$/H$_2$O", ha="right", fontsize=7, color=MUTED)
    ax[2].set(title="Cr(VI)/Cr(III) = 1 boundary", xlabel="pH", ylabel="Eh (V)")
    ax[2].legend(fontsize=7)
    fig.savefig(FIG / "fig1_speciation.png")
    plt.close(fig)


def exp_identifiability(corpus, mats):
    """RQ2: do isotherm/kinetic fits on apparent Cr(VI) data reveal the removal mechanism?"""
    rows = []
    for sid, g in corpus.groupby("study_id"):
        m = mats.set_index("study_id").loc[sid]
        iso = ds.endpoints(g[g.series == "isotherm"]).sort_values("Ce_VI")
        fits = fit.fit_isotherms(iso["Ce_VI"], iso["q_app_VI"], T=float(m["T"]))
        sel = fit.selection_table(fits)
        lang = fits.get("Langmuir")
        ci, _ = fit.bootstrap_ci(lang, iso["Ce_VI"].values, T=float(m["T"]), B=300)
        kin = ds.endpoints(g[g.series == "kinetic"]).sort_values("t")
        kf = fit.selection_table(fit.fit_kinetics(kin["t"], kin["q_app_VI"]))
        hi = iso.iloc[-1]
        true_ads_hi = hi["q_adsVI_true"]
        rows.append({
            "study_id": sid, "family": m["family"], "pH": m["pH_opt"],
            "best_iso": sel.loc[0, "model"], "best_iso_w": sel.loc[0, "Akaike_w"],
            "R2_Langmuir": lang.r2, "R2_best": sel["R2"].max(),
            "runs_p_Langmuir": sel.set_index("model").loc["Langmuir", "runs_p"],
            "qmax_L": lang.params[1], "qmax_L_lo": ci[0, 1], "qmax_L_hi": ci[1, 1],
            "qmax_true": m["qmax"], "q_ads_true_at_maxC0": true_ads_hi, "q_app_at_maxC0": hi["q_app_VI"],
            "qmax_true_in_CI": bool(ci[0, 1] <= m["qmax"] <= ci[1, 1]),
            "best_kin": kf.loc[0, "model"] if len(kf) else None,
            "R2_PSO": kf.set_index("model")["R2"].get("PSO", np.nan),
            # mass-weighted over the isotherm series (row-wise ratios explode when removal ~ 0)
            "red_index_aq": float((iso["Ce_total"] - iso["Ce_VI"]).clip(lower=0).sum() / (iso["C0"] - iso["Ce_VI"]).sum()),
            "true_red_frac": float(iso["reduced_mgL"].sum() / (iso["C0"] - iso["Ce_VI_true"]).sum()),
        })
    r = pd.DataFrame(rows)
    r["qmax_ratio"] = r["qmax_L"] / r["qmax_true"]
    r.to_csv(RES / "T3_identifiability.csv", index=False)
    rho, p = stats.spearmanr(r["true_red_frac"], r["R2_best"])
    rho2, p2 = stats.spearmanr(r["true_red_frac"], r["qmax_ratio"])
    rho3, p3 = stats.spearmanr(r["true_red_frac"], r["red_index_aq"])
    OUT["identifiability"] = {
        "n_studies": len(r), "share_R2best_ge_0.95": float((r["R2_best"] >= 0.95).mean()),
        "share_R2Langmuir_ge_0.95": float((r["R2_Langmuir"] >= 0.95).mean()),
        "median_true_red_frac": float(r["true_red_frac"].median()),
        "range_true_red_frac": [float(r["true_red_frac"].min()), float(r["true_red_frac"].max())],
        "spearman_redfrac_vs_R2best": [float(rho), float(p)],
        "spearman_redfrac_vs_qmax_ratio": [float(rho2), float(p2)],
        "median_qmax_ratio_high_red": float(r.loc[r.true_red_frac > 0.3, "qmax_ratio"].median()),
        "median_qmax_ratio_low_red": float(r.loc[r.true_red_frac < 0.1, "qmax_ratio"].median()),
        "n_high_red": int((r.true_red_frac > 0.3).sum()), "n_low_red": int((r.true_red_frac < 0.1).sum()),
        "share_qmax_true_in_CI": float(r["qmax_true_in_CI"].mean()),
        "n_red_gt_0.2": int((r.true_red_frac > 0.2).sum()),
        "n_red_gt_0.2_and_R2_ge_0.95": int(((r.true_red_frac > 0.2) & (r.R2_best >= 0.95)).sum()),
        "max_red_frac_with_R2_ge_0.95": float(r.loc[r.R2_best >= 0.95, "true_red_frac"].max()),
        "share_best_PSO": float((r["best_kin"] == "PSO").mean()),
        "share_runs_flag_p<0.05": float((r["runs_p_Langmuir"] < 0.05).mean()),
        "spearman_redfrac_vs_aq_reduction_index": [float(rho3), float(p3)],
        "median_aq_index_over_true": float((r["red_index_aq"] / r["true_red_frac"]).median()),
        "best_iso_counts": r["best_iso"].value_counts().to_dict(),
    }
    # The population a reader of the literature sees: fits good enough to be published.
    pub = r[r.R2_best >= 0.95]
    rho4, p4 = stats.spearmanr(pub["true_red_frac"], pub["qmax_ratio"])
    OUT["identifiability"]["publishable_fits"] = {
        "n": len(pub), "median_true_red_frac": float(pub.true_red_frac.median()),
        "n_red_gt_0.2": int((pub.true_red_frac > 0.2).sum()),
        "median_qmax_ratio": float(pub.qmax_ratio.median()),
        "qmax_ratio_range": [float(pub.qmax_ratio.min()), float(pub.qmax_ratio.max())],
        "spearman_redfrac_vs_qmax_ratio": [float(rho4), float(p4)],
        "share_qmax_true_in_CI": float(pub.qmax_true_in_CI.mean()),
        "median_R2": float(pub.R2_best.median())}
    low = r[r.R2_best < 0.95]
    OUT["identifiability"]["poor_fits"] = low[["study_id", "family", "R2_best", "true_red_frac",
                                               "q_app_at_maxC0"]].round(3).to_dict("records")
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.5), constrained_layout=True)
    fams = list(sim.FAMILIES)
    for k, f in enumerate(fams):
        s = r[r.family == f]
        lab = f.replace("_", " ")
        ax[0].scatter(s.true_red_frac, s.R2_best, s=22, color=C[k], edgecolor="white", lw=0.8, label=lab)
        s_ok = s[s.R2_best >= 0.95]
        ax[1].scatter(s_ok.true_red_frac, s_ok.qmax_ratio, s=22, color=C[k], edgecolor="white", lw=0.8)
        ax[2].scatter(s.true_red_frac, s.red_index_aq, s=22, color=C[k], edgecolor="white", lw=0.8)
    ax[0].set(xlabel="True share of removal by reduction", ylabel="Best isotherm $R^2$", ylim=(0.5, 1.02))
    ax[1].axhline(1, color=MUTED, lw=1, ls=":")
    ax[1].set(xlabel="True share of removal by reduction", ylabel="Langmuir $q_{max}$ / true $q_{max}$", title="fits with $R^2 \\geq$ 0.95")
    ax[2].plot([0, 1], [0, 1], color=MUTED, lw=1, ls=":")
    ax[2].set(xlabel="True share of removal by reduction", ylabel="Aqueous Cr(III) index")
    for a_ in ax:
        a_.set_xlabel("True reduction share")
    ax[0].legend(fontsize=6.5, loc="lower left")
    fig.savefig(FIG / "fig2_identifiability.png")
    plt.close(fig)

    # Illustrative study (selected, not representative): largest reduction share among
    # studies whose best isotherm fit has R2 >= 0.95.
    good = r[r.R2_best >= 0.95]
    sid = (good if len(good) else r).sort_values("true_red_frac").iloc[-1]["study_id"]
    g = corpus[corpus.study_id == sid]
    iso = ds.endpoints(g[g.series == "isotherm"]).sort_values("Ce_VI")
    fits = fit.fit_isotherms(iso["Ce_VI"], iso["q_app_VI"])
    xx = np.geomspace(iso.Ce_VI.min(), iso.Ce_VI.max(), 100)
    fig, ax = plt.subplots(figsize=(3.5, 2.6), constrained_layout=True)
    ax.scatter(iso.Ce_VI, iso.q_app_VI, s=24, color=C[0], edgecolor="white", label="apparent Cr(VI) removal")
    ax.scatter(iso.Ce_VI, iso.q_adsVI_true, s=24, color=C[1], edgecolor="white", label="true adsorbed Cr(VI)")
    ax.plot(xx, fit.ISOTHERMS["Langmuir"][0](xx, *fits["Langmuir"].params),
            color=C[0], lw=1.5, label=f"Langmuir fit, $R^2$={fits['Langmuir'].r2:.3f}")
    ax.set(xscale="log", xlabel="$C_e$ Cr(VI) (mg/L)", ylabel="q (mg/g)", title=f"Virtual study {sid}")
    ax.legend(fontsize=6.5)
    fig.savefig(FIG / "fig3_example_study.png")
    plt.close(fig)
    OUT["example_study"] = {"study_id": sid, **r.set_index("study_id").loc[sid, ["true_red_frac", "R2_Langmuir",
                                                                                  "qmax_L", "qmax_true"]].to_dict()}
    return r


def add_features(df):
    df = ds.endpoints(df)
    df = ch.speciation_features(df)
    df["pH_minus_pzc"] = df["pH"] - df["pHpzc"]
    return df


def exp_ml(corpus, quick):
    df = add_features(corpus)
    models = ["RF", "GPR"] + ([] if val.XGBRegressor is None else ["XGB"])
    schemes = ["random", "grouped5", "LOSO", "LOFO"] if not quick else ["random", "LOSO"]
    rows, preds = [], {}
    for mname in models:
        for fs_name, fs in (("base", val.BASE), ("chem", val.CHEM)):
            for sc in schemes:
                t0 = time.time()
                res, pred = val.evaluate(df, mname, fs, sc)
                preds[(mname, fs_name, sc)] = pred
                rows.append(res)
                print(f"  {mname:4s} {fs_name:5s} {sc:8s} R2={res['R2']:.3f} RMSE={res['RMSE']:.2f} ({time.time()-t0:.0f}s)")
    T = pd.DataFrame(rows)
    T.to_csv(RES / "T4_validation.csv", index=False)
    OUT["validation"] = T[["model", "features", "scheme", "RMSE", "MAE", "MedAE", "R2"]].round(3).to_dict("records")
    bsv = val.between_study_variance(df, preds[("GPR", "chem", "LOSO")])
    OUT["between_study_GPRchem_LOSO"] = {k: round(v, 4) for k, v in bsv.items()}

    plot_validation(T, models, schemes)
    return df


def plot_validation(T, models, schemes):
    fig, ax = plt.subplots(figsize=(4.8, 3.0), constrained_layout=True)
    w = 0.12
    combos = [(m, f) for m in models for f in ("base", "chem")]
    for i, (m, f) in enumerate(combos):
        s = T[(T.model == m) & (T.features == f)].set_index("scheme").loc[schemes]
        ax.bar(np.arange(len(schemes)) + (i - (len(combos) - 1) / 2) * w, s["R2"], w * 0.9,
               color=C[models.index(m)], alpha=1.0 if f == "chem" else 0.45,
               label=f"{m} {'+speciation' if f == 'chem' else 'base'}")
    ax.set_xticks(range(len(schemes)), [{"random": "Random\n5-fold", "grouped5": "Grouped\n5-fold",
                                  "LOSO": "Leave-one-\nstudy-out", "LOFO": "Leave-one-\nfamily-out"}[s] for s in schemes])
    ax.set(ylabel="$R^2$ (q$_{app}$, mg/g)", ylim=(min(0, T.R2.min() - 0.05), 1))
    ax.legend(fontsize=6.5, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    fig.savefig(FIG / "fig4_validation.png")
    plt.close(fig)


def exp_uq(df):
    rows = []
    for mname, uq in (("RF", "internal"), ("GPR", "internal"), ("GPR", "conformal"), ("RF", "conformal")):
        res, _ = val.evaluate(df, mname, val.CHEM, "LOSO", uq=uq)
        rows.append(res)
        print(f"  UQ {mname} {uq}: PICP90={res['PICP90']:.3f} ACE={res['ACE']:.3f}")
    T = pd.DataFrame(rows)
    T.to_csv(RES / "T5_uncertainty.csv", index=False)
    OUT["uncertainty"] = T.round(3).to_dict("records")


def exp_missing_chem(df):
    """RQ4: speciation-feature spread when I and the pH basis are not reported."""
    rng = np.random.default_rng(7)
    sub = df[df.series == "pH_edge"].copy()
    spreads = []
    for _, r in sub.iterrows():
        f = []
        for _ in range(40):
            I = 10 ** rng.uniform(-3, -1)
            dpH = rng.uniform(-0.5, 0.5)  # initial vs final pH ambiguity
            f.append(ch.speciate_cr6(float(np.clip(r.pH + dpH, 1, 11)), r.C0, ionic_strength=I).fractions.get("HCrO4-", 0))
        spreads.append({"pH": r.pH, "f_HCrO4_nominal": r.f_HCrO4, "f_p05": np.quantile(f, 0.05),
                        "f_p95": np.quantile(f, 0.95)})
    s = pd.DataFrame(spreads).groupby("pH").mean().reset_index()
    s["width90"] = s.f_p95 - s.f_p05
    s.to_csv(RES / "T6_missing_chemistry.csv", index=False)
    # I-only scenario (pH known)
    w_I = []
    for p in (3.0, 5.0, 6.5, 8.0):
        f = [ch.speciate_cr6(p, 100.0, ionic_strength=I).fractions.get("HCrO4-", 0) for I in (0.001, 0.1)]
        w_I.append({"pH": p, "f_I0.001": f[0], "f_I0.1": f[1]})
    OUT["missing_chem"] = {"pH_and_I_unknown": s.round(3).to_dict("records"),
                           "I_only": pd.DataFrame(w_I).round(3).to_dict("records")}


def exp_optimisation(df, quick):
    """Robust vs nominal design: choose pH and dose for a new material to reach Ce_VI <= 0.5 mg/L
    AND dissolved total Cr <= 2 mg/L from C0 = 20 mg/L, minimising dose."""
    from sklearn.ensemble import RandomForestRegressor
    rng = np.random.default_rng(99)
    feats = val.CHEM
    X = df[feats]
    models = {}
    for tgt in ("q_app_VI", "q_total"):
        y = val.to_target(df[tgt].values)
        models[tgt] = RandomForestRegressor(400, min_samples_leaf=2, random_state=0, n_jobs=-1).fit(X, y)
    # group-conformal half widths at 90% from LOSO residuals
    hw = {}
    for tgt in models:
        _, pred = val.evaluate(df, "RF", feats, "grouped5", target=tgt)
        hw[tgt] = float(np.quantile(np.abs(val.to_target(df[tgt].values) - pred), 0.9))
    n_new = 6 if quick else 15
    outcomes = []
    C0, LIM_VI, LIM_TOT = 20.0, 0.5, 2.0
    for j in range(n_new):
        fam = list(sim.FAMILIES)[j % 3]
        m = sim.random_material(rng, fam, 900 + j)
        grid = [(p, D) for p in np.arange(2.0, 7.01, 0.5) for D in (0.5, 1, 2, 3, 4, 6, 8, 10, 15)]
        G = pd.DataFrame([{"C0": C0, "pH": p, "D": D, "t": 1440.0, "T": 298.15, "I": 0.01, "SBET": m.SBET,
                           "Vp": m.Vp, "dp": m.dp, "pHpzc": m.pHpzc} for p, D in grid])
        G = ch.speciation_features(G)
        G["pH_minus_pzc"] = G["pH"] - G["pHpzc"]
        for tgt, col in (("q_app_VI", "Ce_VI"), ("q_total", "Ce_tot")):
            mu = models[tgt].predict(G[feats])
            for lab, yy in (("nom", mu), ("rob", mu - hw[tgt])):  # robust: lower bound on uptake
                G[f"{col}_{lab}"] = np.clip(C0 - val.from_target(yy) * G["D"], 0, None)
        for lab in ("nom", "rob"):
            ok = G[(G[f"Ce_VI_{lab}"] <= LIM_VI) & (G[f"Ce_tot_{lab}"] <= LIM_TOT)]
            if ok.empty:
                outcomes.append({"material": m.material_id, "strategy": lab, "feasible": False})
                continue
            best = ok.sort_values(["D", "pH"]).iloc[0]
            tr = sim.run_batch(m, C0, best.pH, best.D, 1440.0)
            met = tr["Ce_VI_true"] <= LIM_VI and (tr["Ce_VI_true"] + tr["Cr3_aq_true"]) <= LIM_TOT
            outcomes.append({"material": m.material_id, "strategy": lab, "feasible": True, "pH": best.pH,
                             "D": best.D, "true_Ce_VI": tr["Ce_VI_true"],
                             "true_Ce_total": tr["Ce_VI_true"] + tr["Cr3_aq_true"], "target_met": met})
    O = pd.DataFrame(outcomes)
    O.to_csv(RES / "T7_optimisation.csv", index=False)
    summ = {}
    for lab, g in O.groupby("strategy"):
        f = g[g.feasible]
        summ[lab] = {"n_materials": len(g), "n_feasible": len(f),
                     "share_target_met": float(f["target_met"].mean()) if len(f) else None,
                     "median_dose_gL": float(f["D"].median()) if len(f) else None}
    OUT["optimisation"] = summ


def exp_external(df_train, quick):
    ext, _ = sim.virtual_corpus(n_studies=6 if quick else 9, seed=777)
    ext["study_id"] = "EX" + ext["study_id"]
    ext = add_features(ext)
    rows = []
    for mname in ["RF", "GPR"] + ([] if val.XGBRegressor is None else ["XGB"]):
        for fs_name, fs in (("base", val.BASE), ("chem", val.CHEM)):
            mdl = val.make_model(mname).fit(df_train[fs], val.to_target(df_train["q_app_VI"].values))
            p = np.clip(val.from_target(mdl.predict(ext[fs])), 0, None)
            rows.append({"model": mname, "features": fs_name, **val.metrics(ext["q_app_VI"].values, p)})
    T = pd.DataFrame(rows)
    T.to_csv(RES / "T8_external.csv", index=False)
    OUT["external"] = T.round(3).to_dict("records")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--only", choices=["chem_iso"], help="rerun chemistry + identifiability, merge into summary.json")
    a = ap.parse_args()
    t0 = time.time()
    if a.only == "chem_iso":
        OUT.update(json.loads((RES / "summary.json").read_text()))
        exp_chemistry()
        corpus = pd.read_csv(RES / "virtual_corpus.csv")
        mats = pd.read_csv(RES / "virtual_materials.csv")
        exp_identifiability(corpus, mats)
        (RES / "summary.json").write_text(json.dumps(OUT, indent=2, default=float))
        return
    print("[1] chemistry")
    exp_chemistry()
    corpus, mats = sim.virtual_corpus(n_studies=12 if a.quick else 24, seed=2026)
    corpus.to_csv(RES / "virtual_corpus.csv", index=False)
    mats.to_csv(RES / "virtual_materials.csv", index=False)
    OUT["corpus"] = {"n_studies": int(corpus.study_id.nunique()), "n_rows": len(corpus),
                     "families": corpus.groupby("material_family").study_id.nunique().to_dict()}
    print("[2] identifiability"); exp_identifiability(corpus, mats)
    print("[3] ML validation"); df = exp_ml(corpus, a.quick)
    print("[4] uncertainty"); exp_uq(df)
    print("[5] missing chemistry"); exp_missing_chem(df)
    print("[6] external"); exp_external(df, a.quick)
    print("[7] optimisation"); exp_optimisation(df, a.quick)
    OUT["runtime_s"] = round(time.time() - t0, 1)
    (RES / "summary.json").write_text(json.dumps(OUT, indent=2, default=float))
    print(f"done in {OUT['runtime_s']} s")


if __name__ == "__main__":
    main()
