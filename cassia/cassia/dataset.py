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
    q_total   : total-Cr uptake, (C0 - Ce_total)/D. Only when total Cr was measured.
    c_reduced_aq : Cr(III) remaining in solution, Ce_total - Ce_VI (>= 0).
    """
    out = df.copy()
    out["q_app_VI"] = (out["C0"] - out["Ce_VI"]) / out["D"]
    if "Ce_total" in out:
        out["q_total"] = (out["C0"] - out["Ce_total"]) / out["D"]
        out["c_reduced_aq"] = (out["Ce_total"] - out["Ce_VI"]).clip(lower=0)
        out["reduction_index"] = out["c_reduced_aq"] / (out["C0"] - out["Ce_VI"]).replace(0, np.nan)
    return out


def evidence_class(row: pd.Series) -> str:
    """A: oxidation-state-resolved mass balance + chemistry; B: partial; C: residual Cr(VI) only.

    Class C rows may only be modelled as apparent Cr(VI) disappearance.
    """
    has_total = pd.notna(row.get("Ce_total"))
    has_solid = str(row.get("solid_phase_cr_speciation", "none")).lower() not in ("none", "nan", "")
    has_pH_basis = str(row.get("pH_basis", "")).lower() in ("final", "controlled")
    has_I = pd.notna(row.get("I")) or pd.notna(row.get("electrolyte"))
    if has_total and has_solid and has_pH_basis and has_I:
        return "A"
    if has_total or has_solid:
        return "B"
    return "C"


def classify(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["evidence_class"] = out.apply(evidence_class, axis=1)
    out["endpoint_label"] = np.where(out["evidence_class"] == "C", "apparent_CrVI_disappearance",
                                     "CrVI_disappearance_with_partition")
    return out


def template() -> pd.DataFrame:
    return pd.DataFrame({c: pd.Series(dtype="object") for c in SCHEMA})
