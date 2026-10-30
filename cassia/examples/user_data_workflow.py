"""End-to-end CASSIA workflow for a dataset in the CASSIA schema (data/cassia_schema.csv).

    python examples/user_data_workflow.py path/to/your_data.csv  out_dir/
    python examples/user_data_workflow.py --demo out_dir/        # synthetic demo file, clearly labelled

Steps: schema validation -> observation-level evidence classes -> endpoints (q_app,VI and, where
measured, q_tot,aq) -> PHREEQC speciation features -> per-study isotherm model selection with
bootstrap CIs -> aqueous Cr(III) index (classes A/B) -> study-level CV+ prediction intervals
(when >= 10 studies). Results are written as CSV files plus a short text report.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cassia import chemistry as ch, dataset as ds, fitting as fit, simulate as sim, validation as val  # noqa: E402


def demo_frame():
    """Synthetic records in the schema; every row is marked as synthetic in its provenance fields."""
    df, _ = sim.virtual_corpus(12, seed=123, scenario="E", prefix="DEMO")
    df["doi"] = "synthetic:cassia-demo"
    df["source_location"] = "simulated"
    df["extraction_method"] = "simulated"
    keep = [c for c in ds.SCHEMA if c in df.columns] + ["series"]
    return df[keep]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", nargs="?")
    ap.add_argument("out")
    ap.add_argument("--demo", action="store_true")
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    df = demo_frame() if a.demo else pd.read_csv(a.csv)
    report = [f"CASSIA {__import__('cassia').__version__} workflow; rows = {len(df)}, studies = {df.study_id.nunique()}"]
    if a.demo:
        report.append("INPUT IS SYNTHETIC DEMO DATA - not experimental observations.")

    problems = ds.validate(df)
    report.append("Schema problems: " + ("; ".join(problems) if problems else "none"))
    df = ds.classify(ds.endpoints(df))
    report.append("Observation-level evidence classes: " + df.evidence_class.value_counts().to_dict().__repr__())
    df = ch.speciation_features(df, conc_col="C0")
    df["pH_minus_pzc"] = df["pH"] - df["pHpzc"]
    df["I_imputed"] = df["I_imputed"].astype(float)

    iso_rows = []
    for sid, g in df.groupby("study_id"):
        iso = g[g.get("series", pd.Series("isotherm", index=g.index)) == "isotherm"].sort_values("Ce_VI")
        if len(iso) < 4:
            continue
        fits = fit.fit_isotherms(iso.Ce_VI.values, iso.q_app_VI.values, T=float(iso["T"].iloc[0]))
        tab = fit.selection_table(fits)
        rec = {"study_id": sid, "best_model": tab.loc[0, "model"], "best_R2": tab.R2.max(),
               "endpoint_label": "apparent Cr(VI) removal capacity" if (iso.evidence_class == "C").all()
               else "Cr(VI) disappearance with partition data"}
        if "Langmuir" in fits:
            ci, _ = fit.bootstrap_ci(fits["Langmuir"], iso.Ce_VI.values, T=float(iso["T"].iloc[0]), B=200)
            rec.update(qmax_L=fits["Langmuir"].params[1], qmax_ci=(ci[0, 1], ci[1, 1]))
        if iso.Ce_total.notna().all():
            rec["aq_cr3_index_lower_bound"] = float((iso.Ce_total - iso.Ce_VI).clip(lower=0).sum()
                                                    / (iso.C0 - iso.Ce_VI).sum())
        iso_rows.append(rec)
    pd.DataFrame(iso_rows).to_csv(out / "isotherm_selection.csv", index=False)
    report.append(f"Isotherm selection written for {len(iso_rows)} studies (q_max labelled by evidence class).")

    if df.study_id.nunique() >= 10:
        cp = val.StudyCVPlus("RF", val.CHEM, K=min(10, df.study_id.nunique())).fit(df, "q_app_VI")
        lo, hi, c = cp.interval(df, 0.1, "M1")
        df["q_app_pred"], df["q_app_lo90"], df["q_app_hi90"] = (val.from_target(v) for v in (c, lo, hi))
        report.append("Study-level CV+ (M1) 90% intervals added (in-sample illustration; use LOSO for evaluation).")
    else:
        report.append("Fewer than 10 studies: study-level conformal intervals not computed.")
    df.to_csv(out / "cassia_processed.csv", index=False)
    (out / "report.txt").write_text("\n".join(report) + "\n")
    print("\n".join(report))


if __name__ == "__main__":
    main()
