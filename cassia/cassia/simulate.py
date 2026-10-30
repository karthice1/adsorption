"""Known-truth virtual batch experiments for Cr(VI) on carbon.

The generator couples (i) PHREEQC Cr(VI) speciation, (ii) electrostatically
modulated Langmuir uptake of Cr(VI) with linear-driving-force kinetics,
(iii) acid-promoted reduction of Cr(VI) by a finite pool of surface
reductant, and (iv) partitioning of the Cr(III) produced between solution,
surface binding and Cr(OH)3(am) precipitation (PHREEQC solubility). Every
virtual observation therefore carries the *true* split between adsorbed
Cr(VI), reduced chromium and dissolved Cr(III), which no published study
provides for every data point. It is used only to test whether methods can
recover that split; it is not a substitute for published data.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from . import chemistry as ch

_PH_GRID = np.round(np.arange(1.5, 10.01, 0.25), 2)


class _Tables:
    """PHREEQC look-up tables (built once) used inside the ODE right-hand side."""
    _cache = {}

    @classmethod
    def get(cls, database="minteq.v4.dat"):
        if database not in cls._cache:
            fr = [ch.speciate_cr6(p, 50.0, database=database).fractions for p in _PH_GRID]
            fH = np.array([f.get("HCrO4-", 0) for f in fr])
            fD = np.array([f.get("Cr2O7-2", 0) for f in fr])
            S = np.array([ch.cr3_solubility_mgL(p) for p in _PH_GRID])
            S = np.where(np.isfinite(S), S, 1e6)
            cls._cache[database] = (fH, fD, np.log10(S))
        return cls._cache[database]


def f_hcro4(pH, database="minteq.v4.dat"):
    fH, fD, _ = _Tables.get(database)
    return np.interp(pH, _PH_GRID, fH + fD)  # monovalent/dimer anions: the low-pH, high-affinity forms


def cr3_solubility(pH):
    _, _, logS = _Tables.get()
    return 10 ** np.interp(pH, _PH_GRID, logS)


@dataclass
class Material:
    material_id: str
    family: str
    SBET: float      # m2/g
    Vp: float        # cm3/g
    dp: float        # nm
    pHpzc: float
    qmax: float      # mg Cr(VI)/g true adsorption capacity
    K0: float        # L/mg affinity scale
    R0: float        # mg Cr(VI)/g reducible by the surface
    kred: float      # L g-1 min-1 reduction rate constant at full acid promotion
    pH50_red: float  # pH of half-maximal acid promotion of reduction
    pH50_cr3: float  # pH of half-maximal Cr(III) surface binding
    ka: float = 0.02  # 1/min LDF adsorption rate


FAMILIES = {
    # family: (SBET range, reductant scale, qmax per 1000 m2/g)
    "activated_carbon": ((600, 1800), 0.4, 45.0),
    "biochar": ((50, 500), 1.6, 25.0),
    "biomass_porous_carbon": ((800, 2500), 0.8, 40.0),
}


def random_material(rng, family, idx) -> Material:
    (lo, hi), red_scale, q_per = FAMILIES[family]
    S = rng.uniform(lo, hi)
    pzc = rng.uniform(2.5, 8.5) if family != "biochar" else rng.uniform(6.0, 9.5)
    qmax = q_per * S / 1000 * np.exp(rng.normal(0, 0.35))
    R0 = red_scale * qmax * np.exp(rng.normal(0, 0.5))
    return Material(
        material_id=f"{family[:2].upper()}{idx:03d}", family=family, SBET=S,
        Vp=S / 1000 * rng.uniform(0.35, 0.75), dp=rng.uniform(1.6, 6.0), pHpzc=pzc,
        qmax=qmax, K0=10 ** rng.uniform(-1.3, 0.0), R0=R0,
        kred=10 ** rng.uniform(-4.0, -2.3) * red_scale, pH50_red=rng.uniform(2.0, 3.5),
        pH50_cr3=rng.uniform(3.5, 5.0))


def _affinity(m: Material, pH):
    elec = 1.0 / (1.0 + np.exp((pH - m.pHpzc) / 0.8))      # net positive surface below pHpzc
    return m.K0 * (0.05 + f_hcro4(pH)) * (0.1 + elec)


def run_batch(m: Material, C0, pH, D, t_end, T=298.15, n_t=None):
    """Integrate one batch; returns dict of true end-state (or time series if n_t)."""
    K = _affinity(m, pH)
    h = 1.0 / (1.0 + 10 ** (pH - m.pH50_red))               # acid promotion of reduction
    arr = np.exp(-35000 / 8.314 * (1 / T - 1 / 298.15))       # Ea = 35 kJ/mol for reduction
    arr_a = np.exp(-15000 / 8.314 * (1 / T - 1 / 298.15))
    phi = 1.0 / (1.0 + np.exp(-(pH - m.pH50_cr3) / 0.5))    # Cr(III) bound fraction

    def rhs(_, y):
        Cv, q6, Rs, red = y
        Cv = max(Cv, 0.0)
        qeq = m.qmax * K * Cv / (1 + K * Cv)
        r_ads = m.ka * arr_a * (qeq - q6)                        # mg/g/min
        r_red = m.kred * arr * h * max(Rs, 0) / m.R0 * Cv         # mg/g/min (per g carbon)
        return [-D * (r_ads + r_red), r_ads, -r_red, r_red]

    t_eval = np.linspace(0, t_end, n_t) if n_t else [t_end]
    sol = solve_ivp(rhs, (0, t_end), [C0, 0.0, m.R0, 0.0], t_eval=t_eval, method="LSODA",
                    rtol=1e-7, atol=1e-9)
    Cv, q6, _, red = sol.y
    red_mgL = red * D                                          # Cr reduced, mg/L of solution
    c3_unbound = (1 - phi) * red_mgL
    c3_aq = np.minimum(c3_unbound, cr3_solubility(pH))
    precip = c3_unbound - c3_aq
    out = {"t": sol.t, "Ce_VI_true": Cv, "q_adsVI_true": q6, "reduced_mgL": red_mgL,
           "Cr3_aq_true": c3_aq, "Cr3_bound_mgL": phi * red_mgL, "Cr3_precip_mgL": precip}
    if not n_t:
        out = {k: float(v[-1]) if np.ndim(v) else v for k, v in out.items()}
    return out


def measure(rng, true, sd_rel=0.03, lod=0.02, bias=1.0):
    """Analytical model: DPC for Cr(VI), AAS/ICP for total Cr, multiplicative lab bias."""
    cv = true["Ce_VI_true"] * bias * (1 + rng.normal(0, sd_rel))
    ct = (true["Ce_VI_true"] + true["Cr3_aq_true"]) * bias * (1 + rng.normal(0, sd_rel))
    return max(cv, lod / 2), max(ct, lod / 2)


def virtual_study(rng, study_idx, family=None, n_iso=8, n_edge=6, n_kin=8):
    family = family or rng.choice(list(FAMILIES))
    m = random_material(rng, family, study_idx)
    sid = f"VS{study_idx:02d}"
    bias = np.exp(rng.normal(0, 0.05))
    I = float(rng.choice([0.001, 0.01, 0.05, 0.1]))
    T = float(rng.choice([298.15, 303.15, 308.15]))
    D = float(rng.choice([1.0, 2.0, 4.0]))
    pH_opt = float(rng.uniform(1.8, 3.5))
    t_eq = float(rng.choice([360, 720, 1440]))
    rows = []

    def add(kind, C0, pH, t, **extra):
        tr = run_batch(m, C0, pH, D, t, T)
        cv, ct = measure(rng, tr, bias=bias)
        rows.append({"study_id": sid, "material_id": m.material_id, "material_family": family,
                     "series": kind, "C0": C0, "pH": pH, "D": D, "t": t, "T": T, "I": I,
                     "SBET": m.SBET, "Vp": m.Vp, "dp": m.dp, "pHpzc": m.pHpzc,
                     "Ce_VI": cv, "Ce_total": ct, **tr, **extra})

    for C0 in np.geomspace(10, 400, n_iso):
        add("isotherm", float(C0), pH_opt, t_eq)
    for pH in np.linspace(2.0, 8.0, n_edge):
        add("pH_edge", 100.0, float(pH), t_eq)
    for t in np.geomspace(5, t_eq, n_kin):
        add("kinetic", 100.0, pH_opt, float(t))
    df = pd.DataFrame(rows)
    df["frac_removal_by_reduction"] = df["reduced_mgL"] / (df["C0"] - df["Ce_VI_true"]).clip(lower=1e-9)
    return df, asdict(m) | {"study_id": sid, "I": I, "T": T, "D": D, "pH_opt": pH_opt, "bias": bias}


def virtual_corpus(n_studies=24, seed=2026):
    rng = np.random.default_rng(seed)
    fams = list(FAMILIES)
    dfs, mats = [], []
    for i in range(n_studies):
        df, meta = virtual_study(rng, i, family=fams[i % len(fams)])
        dfs.append(df)
        mats.append(meta)
    return pd.concat(dfs, ignore_index=True), pd.DataFrame(mats)
