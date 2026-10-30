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
DATABASES = ("minteq.v4.dat", "minteq.dat", "llnl.dat")
CR_MW = 51.996  # g/mol

# Cr(VI) species reported in either database (names are database-specific).
CR6_SPECIES = ("CrO4-2", "HCrO4-", "H2CrO4", "Cr2O7-2", "CrO3Cl-", "NaCrO4-", "KCrO4-")
# Cr(III) solids evaluated for precipitation scenarios, by database.
CR3_SOLIDS = {
    "minteq.v4.dat": ("Cr(OH)3(am)", "Cr(OH)3", "Cr2O3"),
    "minteq.dat": ("Cr(OH)3(A)", "Cr(OH)3(C)", "Cr2O3"),
    "llnl.dat": ("Eskolaite",),
}
# Amorphous Cr(III) hydroxide (the kinetically relevant precipitate) per database.
# llnl.dat has no Cr(III) hydroxide; its only Cr(III) solid is written through CrO4-2 + O2.
AMORPHOUS_CR3 = {"minteq.v4.dat": "Cr(OH)3(am)", "minteq.dat": "Cr(OH)3(A)", "llnl.dat": None}
CR3_INPUT = {"minteq.v4.dat": "Cr(3)", "minteq.dat": "Cr(3)", "llnl.dat": "Cr(+3)"}


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
    if database == "llnl.dat":
        raise NotImplementedError("llnl.dat has no Cr(III) hydroxide; use cr3_solubility_si()")
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


def speciation_features_direct(df: pd.DataFrame, database: str = "minteq.v4.dat",
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


# ---------------------------------------------------------------------------
# v0.2: direct saturation-index evaluation (replaces reliance on phase-boundary
# convergence), concentration-/I-/T-dependent speciation look-up, provenance.
# ---------------------------------------------------------------------------

def si_cr3(pH: float, cr3_mgL: float, database: str = "minteq.v4.dat", ionic_strength: float = 0.01,
           temperature_C: float = 25.0, phase: str | None = None):
    """Saturation index of a Cr(III) solid for a solution of given dissolved Cr(III) at fixed pH.

    Returns (SI, converged). No phase-boundary iteration is involved, so a failure here is a
    speciation failure, reported as converged = False rather than interpreted chemically.
    """
    phase = phase or AMORPHOUS_CR3[database]
    if phase is None:
        return np.nan, False
    pp = engine(database)
    try:
        sol = _add_solution(pp, {"pH": pH, "temp": temperature_C, CR3_INPUT[database]: cr3_mgL / CR_MW},
                            "NaCl", ionic_strength)
    except RuntimeError:
        return np.nan, False
    try:
        si = sol.si(phase)
        return (float(si), si > -998)
    finally:
        sol.forget()


def cr3_solubility_si(pH: float, database: str = "minteq.v4.dat", ionic_strength: float = 0.01,
                      temperature_C: float = 25.0, c_min: float = 1e-6, c_max: float = 5000.0,
                      tol_log: float = 0.01, phase: str | None = None) -> dict:
    """Dissolved Cr(III) at saturation (SI = 0) by bisection on log10 C, with diagnostics.

    status: 'bracketed'                  solubility found inside [c_min, c_max]
            'undersaturated_to_cmax'     SI < 0 even at c_max: no constraint up to c_max (a
                                         thermodynamic statement, evaluated directly)
            'supersaturated_at_cmin'     SI > 0 at c_min
            'speciation_failed'          PHREEQC failed at a bracketing point
    """
    si_lo, ok_lo = si_cr3(pH, c_min, database, ionic_strength, temperature_C, phase)
    # Very high Cr(III) at near-neutral pH can defeat the speciation solver itself; step the
    # upper bracket down (recorded in c_max_used) rather than interpreting the failure.
    for c_try in (c_max, c_max / 5, c_max / 25, c_max / 125):
        si_hi, ok_hi = si_cr3(pH, c_try, database, ionic_strength, temperature_C, phase)
        if ok_hi:
            c_max = c_try
            break
    lo, hi = np.log10(c_min), np.log10(c_max)
    out = {"pH": pH, "database": database, "I": ionic_strength, "T_C": temperature_C,
           "SI_at_cmin": si_lo, "SI_at_cmax": si_hi, "c_max_used": c_max, "n_iter": 0}
    if not (ok_lo and ok_hi):
        return out | {"S_mgL": np.nan, "status": "speciation_failed"}
    if si_hi < 0:
        return out | {"S_mgL": np.inf, "status": "undersaturated_to_cmax"}
    if si_lo > 0:
        return out | {"S_mgL": c_min, "status": "supersaturated_at_cmin"}
    n = 0
    while hi - lo > tol_log and n < 60:
        mid = 0.5 * (lo + hi)
        si, ok = si_cr3(pH, 10 ** mid, database, ionic_strength, temperature_C, phase)
        if not ok:
            return out | {"S_mgL": np.nan, "status": "speciation_failed", "n_iter": n}
        lo, hi = (mid, hi) if si < 0 else (lo, mid)
        n += 1
    S = 10 ** (0.5 * (lo + hi))
    si_S, _ = si_cr3(pH, S, database, ionic_strength, temperature_C, phase)
    return out | {"S_mgL": S, "status": "bracketed", "n_iter": n, "SI_at_S": si_S}


def database_provenance() -> pd.DataFrame:
    """SHA-256, size and the $Id revision line of each bundled database."""
    import hashlib
    rows = []
    for db in DATABASES:
        f = DB_DIR / db
        txt = f.read_text(errors="ignore").splitlines()[:40]
        rev = next((l.strip("# ").strip() for l in txt if "$Id" in l), "")
        rows.append({"database": db, "sha256": hashlib.sha256(f.read_bytes()).hexdigest(),
                     "bytes": f.stat().st_size, "revision": rev,
                     "source": "https://github.com/usgs-coupled/phreeqc3/tree/master/database"})
    return pd.DataFrame(rows)


def logk_table() -> pd.DataFrame:
    """log K (25 C) of the Cr(VI) protonation/dimerisation reactions in each database."""
    import re
    rx = {"HCrO4-": [r"CrO4-2 \+ H\+ = HCrO4-", r"H\+ \+ CrO4-2 = HCrO4-"],
          "H2CrO4": [r"CrO4-2 \+ 2 H\+ = H2CrO4", r"2 H\+ \+ CrO4-2 = H2CrO4"],
          "Cr2O7-2": [r"2 CrO4-2 \+ 2 H\+ = Cr2O7-2 \+ H2O", r"2 H\+ \+ 2 CrO4-2 = Cr2O7-2 \+ H2O"]}
    rows = []
    for db in DATABASES:
        lines = (DB_DIR / db).read_text(errors="ignore").splitlines()
        for sp, pats in rx.items():
            for i, l in enumerate(lines):
                if any(re.fullmatch(r"\s*" + p + r"\s*", l) for p in pats):
                    for l2 in lines[i + 1:i + 4]:
                        m = re.search(r"log_?k\s+(-?[\d.]+)", l2)
                        if m:
                            rows.append({"database": db, "species": sp, "logK": float(m.group(1))})
                            break
                    break
    return pd.DataFrame(rows).pivot(index="species", columns="database", values="logK")


class SpeciationGrid:
    """Cr(VI) speciation on a regular (pH, log10 C, log10 I, T) grid with multilinear interpolation.

    C is the *current* dissolved Cr(VI) concentration (mg/L), so dimerisation follows the
    solution as Cr(VI) is removed. Species: f_HCrO4, f_CrO4, f_Cr2O7 (per Cr atom), f_H2CrO4.
    """
    SPECIES = ("HCrO4-", "CrO4-2", "Cr2O7-2", "H2CrO4")

    def __init__(self, database="minteq.v4.dat",
                 pH=np.round(np.arange(1.5, 10.01, 0.25), 2),
                 C=(0.01, 0.1, 1.0, 5.0, 10.0, 25.0, 50.0, 100.0, 200.0, 400.0, 700.0),
                 I=(0.001, 0.003, 0.01, 0.03, 0.1), T=(288.15, 298.15, 308.15, 318.15), cache=True):
        self.database = database
        self.axes = [np.asarray(pH, float), np.log10(np.asarray(C, float)),
                     np.log10(np.asarray(I, float)), np.asarray(T, float)]
        f = DB_DIR.parent / "results" / f"speciation_grid_{database.replace('.', '_')}.npz"
        if cache and f.exists():
            z = np.load(f)
            if all(np.array_equal(z[f"a{k}"], self.axes[k]) for k in range(4)):
                self.values = z["values"]
                return
        shape = tuple(len(a) for a in self.axes)
        vals = np.zeros(shape + (len(self.SPECIES),))
        for i, p in enumerate(self.axes[0]):
            for j, lc in enumerate(self.axes[1]):
                for k, li in enumerate(self.axes[2]):
                    for l, t in enumerate(self.axes[3]):
                        fr = speciate_cr6(float(p), float(10 ** lc), database=database,
                                          ionic_strength=float(10 ** li), temperature_C=float(t - 273.15)).fractions
                        vals[i, j, k, l] = [fr.get(s, 0.0) for s in self.SPECIES]
        self.values = vals
        if cache:
            f.parent.mkdir(exist_ok=True)
            np.savez_compressed(f, values=vals, **{f"a{k}": self.axes[k] for k in range(4)})

    def __call__(self, pH, C, I=0.01, T=298.15):
        """Multilinear interpolation (inputs clipped to the grid); returns array (..., 4)."""
        x = [np.atleast_1d(np.asarray(v, float)) for v in (pH, np.log10(np.maximum(C, 1e-6)), np.log10(I), T)]
        x = np.broadcast_arrays(*x)
        idx, w = [], []
        for d, a in enumerate(self.axes):
            v = np.clip(x[d], a[0], a[-1])
            i = np.clip(np.searchsorted(a, v) - 1, 0, len(a) - 2)
            idx.append(i)
            w.append((v - a[i]) / (a[i + 1] - a[i]))
        out = 0.0
        for corner in range(16):
            bits = [(corner >> d) & 1 for d in range(4)]
            wt = np.ones_like(w[0])
            for d in range(4):
                wt = wt * (w[d] if bits[d] else 1 - w[d])
            out = out + wt[..., None] * self.values[idx[0] + bits[0], idx[1] + bits[1], idx[2] + bits[2], idx[3] + bits[3]]
        return out

    def validate(self, n=200, seed=0) -> dict:
        """Max/mean absolute interpolation error against direct PHREEQC at random interior points."""
        rng = np.random.default_rng(seed)
        errs = []
        for _ in range(n):
            p = rng.uniform(self.axes[0][0], self.axes[0][-1])
            c = 10 ** rng.uniform(self.axes[1][0], self.axes[1][-1])
            I = 10 ** rng.uniform(self.axes[2][0], self.axes[2][-1])
            T = rng.uniform(self.axes[3][0], self.axes[3][-1])
            fr = speciate_cr6(p, c, database=self.database, ionic_strength=I, temperature_C=T - 273.15).fractions
            direct = np.array([fr.get(s, 0.0) for s in self.SPECIES])
            errs.append(np.abs(self(p, c, I, T)[0] - direct))
        e = np.array(errs)
        return {"n": n, "max_abs_err": float(e.max()), "mean_abs_err": float(e.mean()),
                "p99_abs_err": float(np.quantile(e.max(1), 0.99))}


_GRIDS = {}


def grid(database="minteq.v4.dat") -> SpeciationGrid:
    if database not in _GRIDS:
        _GRIDS[database] = SpeciationGrid(database)
    return _GRIDS[database]


def speciation_features(df: pd.DataFrame, database: str = "minteq.v4.dat", default_I: float = 0.01,
                        conc_col: str = "C0", I_col: str = "I", T_col: str = "T") -> pd.DataFrame:
    """Append speciation descriptors (pH, concentration, I and T dependent) from the PHREEQC grid.

    Missing I is set to ``default_I`` and flagged in ``I_imputed``; it is never written back to
    the source data.
    """
    out = df.copy()
    I = out[I_col] if I_col in out else pd.Series(np.nan, index=out.index)
    out["I_imputed"] = I.isna()
    Iv = I.fillna(default_I).values.astype(float)
    Tv = out[T_col].values.astype(float) if T_col in out else np.full(len(out), 298.15)
    f = grid(database)(out["pH"].values.astype(float), out[conc_col].values.astype(float), Iv, Tv)
    out["f_HCrO4"], out["f_CrO4"], out["f_Cr2O7"], out["f_H2CrO4"] = f[:, 0], f[:, 1], f[:, 2], f[:, 3]
    out["z_mean_CrVI"] = -(out["f_HCrO4"] + 2 * out["f_CrO4"] + out["f_Cr2O7"])
    out["log_I"] = np.log10(Iv)
    return out
