"""Published pore-characterisation ingestion (pyGAPS adapter).

pyGAPS is used here only for what it was designed for: gas (N2, 77 K) isotherms
digitised or tabulated in the source publications. It is never used to fit aqueous
Cr(VI) data. When a paper reports SBET but no isotherm, the reported value is kept
and flagged ``SBET_source = reported``.
"""
from __future__ import annotations

import numpy as np

N2_CROSS_SECTION_NM2 = 0.162  # IUPAC value for N2 at 77 K


def bet_from_published_isotherm(p_rel, n_mmol_g, p_limits=None):
    """BET area (m2/g) and diagnostics from a published N2 isotherm (relative pressure, mmol/g).

    Uses the Rouquerol consistency criteria implemented in pyGAPS to choose the linear
    range when ``p_limits`` is None (important for microporous carbons, where the
    classical 0.05-0.30 window is often inappropriate).
    """
    from pygaps.characterisation.area_bet import area_BET_raw

    p = np.asarray(p_rel, float)
    n = np.asarray(n_mmol_g, float) / 1000.0  # mol/g
    order = np.argsort(p)
    area, c, n_m, p_m, slope, intercept, pmin, pmax, r = area_BET_raw(
        p[order], n[order], N2_CROSS_SECTION_NM2, p_limits=p_limits)
    return {"SBET": float(area), "C_BET": float(c), "n_m_mmol_g": float(n_m * 1000),
            "p_range": (float(p[order][pmin]), float(p[order][pmax])), "r_BET": float(r),
            "SBET_source": "pyGAPS_from_isotherm"}


def material_descriptors(reported: dict, isotherm: tuple | None = None) -> dict:
    """Merge reported textural properties with recomputed ones, keeping provenance."""
    out = {k: reported.get(k, np.nan) for k in ("SBET", "Vp", "dp", "pHpzc")}
    out["SBET_source"] = "reported" if np.isfinite(out["SBET"]) else "missing"
    if isotherm is not None:
        bet = bet_from_published_isotherm(*isotherm)
        if np.isfinite(out["SBET"]):
            out["SBET_reported"] = out["SBET"]
            out["SBET_rel_diff"] = (bet["SBET"] - out["SBET"]) / out["SBET"]
        out.update(bet)
    return out
