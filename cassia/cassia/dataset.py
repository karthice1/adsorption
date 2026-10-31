"""Master-dataset schema, endpoint definitions and evidence classification."""
from __future__ import annotations

import numpy as np
import pandas as pd

# column -> (unit, required)
SCHEMA = {
    "study_id": ("-", True), "material_id": ("-", True), "material_family": ("-", True),
    "doi": ("-", True), "source_location": ("table/figure/SI", True), "extraction_method": ("table|digitised|SI", True),
    "C0": ("mg/L", True), "Ce_VI": ("mg/L", False), "Ce_total": ("mg/L", False),
    "pH": ("-", True), "pH_basis": ("initial|final|controlled", False),
    "D": ("g/L", True), "t": ("min", True), "T": ("K", True), "I": ("mol/L", False),
    "electrolyte": ("-", False), "SBET": ("m2/g", False), "Vp": ("cm3/g", False), "dp": ("nm", False),
    "pHpzc": ("-", False), "oxidation_state": ("-", False), "method": ("-", False),
    "cr3_measured": ("bool", False), "solid_phase_cr_speciation": ("XPS|XANES|none", False),
    "xps_cr3_frac": ("-", False), "xanes_cr3_frac": ("-", False), "solid_cr_mgg": ("mg/g", False),
    "cr3_precip_frac": ("-", False), "pH_final": ("-", False), "Ce_VI_below_LOD": ("bool", False),
    "alkalinity": ("mg/L CaCO3", False), "major_anions": ("-", False), "major_cations": ("-", False),
    "Eh_mV": ("mV", False), "reducing_agents": ("-", False),
    "n": ("-", False), "SD": ("same as endpoint", False),
}

ANALYTICAL_CR6 = ("DPC", "diphenylcarbazide", "IC-PCR")
ANALYTICAL_TOTAL = ("AAS", "ICP-OES", "ICP-MS", "FAAS")


def validate(df: pd.DataFrame) -> list[str]:
    problems = [f"missing required column {c}" for c, (_, req) in SCHEMA.items() if req and c not in df]
    for c in ("C0", "D", "t", "T"):
        if c in df and (df[c] <= 0).any():
            problems.append(f"{c} must be > 0")
    if {"Ce_VI", "C0"} <= set(df) and (df["Ce_VI"] > df["C0"] * 1.05).any():
        problems.append("Ce_VI exceeds C0 by >5% for some rows")
    if {"Ce_VI", "Ce_total"} <= set(df):
        bad = df["Ce_total"].notna() & df["Ce_VI"].notna() & (df["Ce_VI"] > df["Ce_total"] * 1.05)
        if bad.any():
            problems.append(f"{int(bad.sum())} rows with Ce_VI > Ce_total (analytical inconsistency)")
    if "T" in df and (df["T"] < 200).any():
        problems.append("T looks like Celsius; schema requires K")
    return problems


def endpoints(df: pd.DataFrame) -> pd.DataFrame:
    """Compute removal endpoints without conflating them.

    q_app_VI  : apparent Cr(VI) disappearance per gram, (C0 - Ce_VI)/D. Always computable.
    q_tot_aq  : total *aqueous* Cr disappearance, (C0 - Ce_total)/D. Only when total Cr was
                measured. It is not adsorption either: it includes Cr(III) retained or
                precipitated on the solid.
    c_reduced_aq : Cr(III) remaining in solution, Ce_total - Ce_VI (>= 0).
    """
    out = df.copy()
    out["q_app_VI"] = (out["C0"] - out["Ce_VI"]) / out["D"]
    if "Ce_total" in out:
        out["q_tot_aq"] = (out["C0"] - out["Ce_total"]) / out["D"]
        out["c_reduced_aq"] = (out["Ce_total"] - out["Ce_VI"]).clip(lower=0)
        out["reduction_index"] = out["c_reduced_aq"] / (out["C0"] - out["Ce_VI"]).replace(0, np.nan)
    return out


def evidence_class(row: pd.Series) -> str:
    """Observation-level evidence class (v0.3).

    C   residual Cr(VI) only (or solid data without a total-aqueous-Cr analysis)
    B   Cr(VI) and total aqueous Cr
    A1  B + solid-phase oxidation-state evidence (XPS or XANES), final/controlled pH and reported
        electrolyte
    A2  A1 + solid Cr inventory (digestion) + phase/coordination evidence distinguishing retained
        from precipitated Cr(III)

    Oxidation-state evidence and chromium mass balance are different analytical tasks: A1 shows
    that Cr(III) is present on the solid; only A2 quantifies how much Cr sits where. Class C rows may
    only be modelled as apparent Cr(VI) disappearance. A publication can contain rows of different
    classes.
    """
    has_total = pd.notna(row.get("Ce_total"))
    has_ox = (str(row.get("solid_phase_cr_speciation", "none")).lower() not in ("none", "nan", "")
              or pd.notna(row.get("xps_cr3_frac")) or pd.notna(row.get("xanes_cr3_frac")))
    has_inv = pd.notna(row.get("solid_cr_mgg"))
    has_phase = pd.notna(row.get("cr3_precip_frac"))
    has_pH_basis = (str(row.get("pH_basis", "")).lower() in ("final", "controlled")
                    or pd.notna(row.get("pH_final")))
    has_I = pd.notna(row.get("I")) or pd.notna(row.get("electrolyte"))
    if not has_total:
        return "C"
    if has_ox and has_pH_basis and has_I:
        return "A2" if (has_inv and has_phase) else "A1"
    return "B"


def classify(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["evidence_class"] = out.apply(evidence_class, axis=1)
    # study-level label = best class available in the study (for reporting only)
    if "study_id" in out:
        order = {"A2": 0, "A1": 1, "B": 2, "C": 3}
        inv = {v: k for k, v in order.items()}
        rank = out["evidence_class"].map(order)
        out["study_best_class"] = rank.groupby(out["study_id"]).transform("min").map(inv)
    out["endpoint_label"] = np.where(out["evidence_class"] == "C", "apparent_CrVI_disappearance",
                                     "CrVI_disappearance_with_partition")
    return out


def template() -> pd.DataFrame:
    return pd.DataFrame({c: pd.Series(dtype="object") for c in SCHEMA})
