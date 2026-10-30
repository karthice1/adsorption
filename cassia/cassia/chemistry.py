"""Chromium aqueous chemistry adapter built on phreeqpython (PHREEQC).

The adapter keeps the oxidation state explicit. Initial-solution speciation
with ``Cr(6)`` fixes chromium as Cr(VI) (PHREEQC does not re-equilibrate
valence states entered separately), which is the correct description of a
freshly prepared dichromate/chromate test solution. Redox *scenarios* are
run separately by entering total ``Cr`` together with an assumed pe, so that
the Cr(VI)/Cr(III) partition is an explicit, labelled assumption.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from phreeqpython import PhreeqPython

DB_DIR = Path(__file__).resolve().parent.parent / "databases"
DATABASES = ("minteq.v4.dat", "llnl.dat")
CR_MW = 51.996  # g/mol

# Cr(VI) species reported in either database (names are database-specific).
CR6_SPECIES = ("CrO4-2", "HCrO4-", "H2CrO4", "Cr2O7-2", "CrO3Cl-", "NaCrO4-", "KCrO4-")
# Cr(III) solids evaluated for precipitation scenarios, by database.
CR3_SOLIDS = {
    "minteq.v4.dat": ("Cr(OH)3(am)", "Cr(OH)3", "Cr2O3"),
    "llnl.dat": ("Eskolaite",),
}


@lru_cache(maxsize=None)
def engine(database: str) -> PhreeqPython:
    if database not in DATABASES:
        raise ValueError(f"database {database!r} has not been validated for Cr; use one of {DATABASES}")
    return PhreeqPython(database=database, database_directory=DB_DIR)


def _background(electrolyte: str, ionic_strength: float, balance: str = "anion") -> dict:
    """Background electrolyte in mmol/kgw.

    pH adjustment is represented by charge balance: on the anion (acid added, e.g. HCl)
    or, when that is infeasible because the chromate charge exceeds the cation pool,
    on the cation (base added, e.g. NaOH).
    """
    I_mmol = max(ionic_strength, 1e-6) * 1000.0
    cation, anion = {"NaCl": ("Na", "Cl"), "NaNO3": ("Na", "N(5)"), "KCl": ("K", "Cl")}[electrolyte]
    if balance == "anion":
        return {cation: I_mmol, anion: f"{I_mmol} charge"}
    return {cation: f"{I_mmol} charge", anion: I_mmol}


def _add_solution(pp, comp_core: dict, electrolyte: str, ionic_strength: float):
    last = None
    for balance in ("anion", "cation"):
        comp = dict(comp_core)
        comp.update(_background(electrolyte, ionic_strength, balance))
        try:
            return pp.add_solution(comp)
        except Exception as e:  # PHREEQC convergence / negative-concentration failure
            last = e
    raise RuntimeError(f"PHREEQC failed for {comp_core}: {last}")


@dataclass
class Speciation:
    pH: float
    cr_total_mgL: float
    fractions: dict  # species -> fraction of total Cr (as Cr atoms)
    ionic_strength: float
    si: dict

    def as_row(self) -> dict:
        row = {"pH": self.pH, "Cr_mgL": self.cr_total_mgL, "I_calc": self.ionic_strength}
        row.update({f"f_{k}": v for k, v in self.fractions.items()})
        row.update({f"SI_{k}": v for k, v in self.si.items()})
        return row


def _cr_atoms(species: str) -> int:
    return 2 if species.startswith("Cr2") else 3 if species.startswith("Cr3") else 1


def speciate_cr6(pH: float, cr_mgL: float, database: str = "minteq.v4.dat",
                 ionic_strength: float = 0.01, electrolyte: str = "NaCl",
                 temperature_C: float = 25.0) -> Speciation:
    """Speciate a Cr(VI) solution at fixed pH (pH adjusted by charge balance on the anion)."""
    pp = engine(database)
    cr_mmol = cr_mgL / CR_MW
    # Cr(VI) added as K2Cr2O7: one K+ per Cr.
    sol = _add_solution(pp, {"pH": pH, "temp": temperature_C, "Cr(6)": cr_mmol, "K": cr_mmol},
                        electrolyte, ionic_strength)
    try:
        sp = sol.species  # mol/kgw
        tot = sum(_cr_atoms(k) * v for k, v in sp.items() if k.startswith("Cr") or "Cr" in k)
        fr = {k: _cr_atoms(k) * v / tot for k, v in sp.items() if "Cr" in k and _cr_atoms(k) * v / tot > 1e-8}
        return Speciation(pH, cr_mgL, fr, sol.I, {})
    finally:
        sol.forget()


def redox_scenario(pH: float, cr_mgL: float, pe: float, database: str = "minteq.v4.dat",
                   ionic_strength: float = 0.01, electrolyte: str = "NaCl") -> dict:
    """Equilibrium Cr(VI)/Cr(III) partition and Cr(III)-solid saturation at an *assumed* pe.

    This is a scenario calculation, not a prediction: pe is almost never reported
    in batch adsorption papers and is varied over a labelled range.
    """
    pp = engine(database)
    sol = _add_solution(pp, {"pH": pH, "pe": pe, "Cr": cr_mgL / CR_MW, "K": cr_mgL / CR_MW},
                        electrolyte, ionic_strength)
    try:
        sp = sol.species
        cr6 = sum(_cr_atoms(k) * v for k, v in sp.items() if k in CR6_SPECIES)
        crall = sum(_cr_atoms(k) * v for k, v in sp.items() if "Cr" in k)
        phases = sol.phases
        si = {p: phases.get(p, np.nan) for p in CR3_SOLIDS[database]}
        return {"pH": pH, "pe": pe, "database": database, "frac_CrVI": cr6 / crall,
                "frac_CrIII": 1 - cr6 / crall, **{f"SI_{k}": v for k, v in si.items()}}
    finally:
        sol.forget()


def cr3_solubility_mgL(pH: float, database: str = "minteq.v4.dat", phase: str = "Cr(OH)3(am)",
                       ionic_strength: float = 0.01) -> float:
    """Dissolved Cr(III) in equilibrium with a Cr(III) solid at fixed pH (mg/L).

    Uses the PHREEQC SOLUTION phase-boundary option (concentration adjusted to SI = 0).
    Returns ``inf`` where the solid cannot form at any relevant concentration (the
    phase-boundary iteration does not converge, e.g. pH <= 4 for Cr(OH)3(am)), which
    in practice means no precipitation constraint. Only minteq.v4.dat is supported:
    llnl.dat writes Eskolaite in terms of CrO4-2 and O2, so a Cr(III)-only
    solubility is not defined without fixing the redox state.
    """
    if database != "minteq.v4.dat":
        raise NotImplementedError("Cr(III) solubility is implemented for minteq.v4.dat only")
    pp = engine(database)
    comp = {"pH": pH, "Cr(3)": f"0.001 {phase} 0"}
    comp.update(_background("NaCl", ionic_strength))
    try:
        sol = pp.add_solution(comp)
    except Exception:
        return np.inf
    try:
        return sol.total_element("Cr", units="mmol") * CR_MW
    finally:
        sol.forget()


def speciation_table(pH_grid, cr_mgL=50.0, databases=DATABASES, **kw) -> pd.DataFrame:
    rows = []
    for db in databases:
        for pH in pH_grid:
            r = speciate_cr6(pH, cr_mgL, database=db, **kw).as_row()
            r["database"] = db
            rows.append(r)
    return pd.DataFrame(rows).fillna(0.0)


def speciation_features(df: pd.DataFrame, database: str = "minteq.v4.dat",
                        ionic_strength_col: str = "I", default_I: float = 0.01) -> pd.DataFrame:
    """Append PHREEQC-derived Cr(VI) speciation descriptors to a dataset (one call per unique pH/C0/I)."""
    out = df.copy()
    keys = []
    for _, r in out.iterrows():
        I = r.get(ionic_strength_col, np.nan)
        keys.append((round(float(r["pH"]), 2), round(float(r["C0"]), 3), default_I if pd.isna(I) else float(I)))
    cache = {}
    for k in set(keys):
        s = speciate_cr6(k[0], k[1], database=database, ionic_strength=k[2])
        cache[k] = s
    f = lambda k, sp: cache[k].fractions.get(sp, 0.0)
    out["f_HCrO4"] = [f(k, "HCrO4-") for k in keys]
    out["f_CrO4"] = [f(k, "CrO4-2") for k in keys]
    out["f_Cr2O7"] = [f(k, "Cr2O7-2") for k in keys]
    # Mean charge of Cr(VI) per Cr atom: electrostatically relevant descriptor.
    out["z_mean_CrVI"] = -(out["f_HCrO4"] + 2 * out["f_CrO4"] + 1.0 * out["f_Cr2O7"])
    out["I_calc"] = [cache[k].ionic_strength for k in keys]
    return out
