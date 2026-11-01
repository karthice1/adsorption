"""Known-truth, multi-mechanism virtual batch experiments for Cr(VI) on carbon (CASSIA v0.2).

One general kinetic model contains five removal pathways; the benchmark scenarios
switch pathways on or off:

  A  adsorption only (negative control)
  B  direct surface-mediated reduction of dissolved Cr(VI) (electron transfer at the
     surface without a prior adsorbed state), product partly released
  C  adsorption-coupled reduction: adsorbed Cr(VI) is reduced in place, Cr(III) is
     released kinetically
  D  dissolved-phase reduction by reductants leached from the carbon, followed by
     Cr(III) partitioning
  E  all pathways + Cr(III) re-adsorption + Cr(OH)3(am) precipitation/dissolution +
     proton-balance pH drift

States (per litre of suspension unless stated): Cv dissolved Cr(VI) [mg/L]; q6 adsorbed
Cr(VI) [mg/g]; q3s Cr(III) formed on the surface [mg/g]; q3r Cr(III) re-adsorbed from
solution [mg/g]; C3 dissolved Cr(III) [mg/L]; P3 precipitated Cr(III) [mg/L]; Rs surface
reductant [mmol e-/g]; L dissolved reductant [mmol e-/L]; Lrel reductant leached so far
[mmol e-/L]; pH.

Conservation (unit-tested): Cv + D(q6 + q3s + q3r) + C3 + P3 = C0, and
3 * (Cr reduced, mmol/L) = D (R0 - Rs) + (Lrel - L).

Speciation (pH, current Cr(VI), I, T) and Cr(III) solubility (pH, I) come from PHREEQC
(minteq.v4.dat) through look-up grids. Every parameter and its distribution is listed in
``PARAMETERS`` (exported to the supplementary table). Family-level means of the
reduction-related parameters are drawn per corpus from a common hyper-prior, so which
carbon family is most reductive is itself random across replicates.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict, field

import numpy as np
import pandas as pd

from . import chemistry as ch
from . import measurement as ms

MW = ch.CR_MW
SCENARIOS = ("A", "B", "C", "D", "E")

FAMILIES = {
    # family: (SBET range m2/g, pHpzc range) - physical plausibility only
    "AC_chemical": ((800, 1800), (2.5, 6.0)),
    "AC_physical": ((400, 1200), (6.0, 9.5)),
    "biochar_wood": ((50, 500), (6.5, 9.5)),
    "biochar_manure": ((20, 200), (7.0, 10.0)),
    "hydrochar": ((5, 60), (4.0, 7.0)),
    "porous_carbon": ((800, 2500), (3.0, 8.0)),
}

# name: (unit, distribution text, sampler) - sampler(rng, fam_hyper) -> value
PARAMETERS = {
    "qmax": ("mg/g", "q_f * SBET/1000 * LN(0, 0.35); q_f ~ LN(ln 35, 0.3) per family and corpus"),
    "K0": ("L/mg", "10^U(-1.3, 0)"),
    "ka": ("1/min", "10^U(-2.3, -1.3)"),
    "R0": ("mmol e-/g", "10^(mu_R,f + N(0, 0.3)); mu_R,f ~ U(-1.3, 0.3) per family and corpus"),
    "kB": ("L/(g min)", "10^U(-4.0, -2.3) (scenarios B, E)"),
    "kC": ("1/min", "10^U(-3.5, -1.8) (scenarios C, E)"),
    "L0": ("mmol e-/g", "10^(mu_L,f + N(0, 0.3)); mu_L,f ~ U(-2.5, -0.8) per family and corpus (D, E)"),
    "kl": ("1/min", "10^U(-3.0, -1.5)"),
    "kD": ("L/(mmol e- min)", "10^U(-2.0, 0.0) (scenarios D, E)"),
    "pH50_red": ("-", "U(2.0, 4.0); acid promotion h(pH) = 1/(1 + 10^(pH - pH50_red))"),
    "Ea_red": ("kJ/mol", "U(25, 50)"),
    "pH50_cr3": ("-", "U(3.5, 5.0); Cr(III) surface affinity theta3(pH) = 1/(1 + exp(-(pH - pH50_cr3)/0.5))"),
    "krel": ("1/min", "10^U(-3.0, -1.5); release of surface Cr(III) at rate krel (1 - theta3) q3s"),
    "q3max": ("mg/g", "0.5 qmax LN(0, 0.3)"),
    "K3": ("L/mg", "10^U(-1.0, 0.5)"),
    "k3": ("1/min", "10^U(-2.5, -1.0)"),
    "kp": ("1/min", "0.05 (fixed); precipitation/dissolution toward PHREEQC Cr(OH)3(am) solubility"),
    "nH": ("mol H+/mol Cr", "U(3, 5) net protons consumed per Cr reduced (E only)"),
    "beta": ("mmol/(L pH)", "10^U(0, 1) buffer capacity (E only)"),
    "kbuf": ("1/min", "10^U(-4.5, -3.0); proton exchange with the carbon, flux kbuf*beta*(pHpzc - pH), buffered by beta + 2.303[H+] (E only)"),
    "Ea_ads": ("kJ/mol", "15 (fixed)"),
}

def _measurement_spec() -> dict:
    """Measurement-model rows of Table S1, generated from the constants used by the code."""
    x, n = ms.TECHNIQUES["XPS"], ms.TECHNIQUES["XANES"]
    return {
        "Ce_VI": "DPC: x(1 + N(0, 0.03)) * study bias; LOD 0.01 mg/L, values < LOD reported as LOD/2 and flagged",
        "Ce_total": "AAS/ICP: (Cv + C3)(1 + N(0, 0.03)) * study bias; LOD 0.02 mg/L",
        "xps_cr3_frac": (f"XPS (surface): f_obs = f_s + b (1 - f_s) + N(0, {x['noise']}), clipped to [0, 1]; "
                         f"f_s = e f/(e f + 1 - f); per-study enrichment e ~ LN(0, {x['enrich_logsd']}), "
                         f"beam/vacuum reduction b ~ U(0, {x['beam_max']})"),
        "xanes_cr3_frac": (f"XANES (bulk, LCF): f_obs = f + b (1 - f) + N(0, {n['noise']}), clipped; "
                           f"b ~ U(0, {n['beam_max']}); no enrichment"),
        "solid_cr_mgg": f"solid Cr inventory (digestion + total-Cr analysis): true inventory (1 + N(0, {ms.INVENTORY_REL_SD}))",
        "cr3_precip_frac": f"phase/coordination evidence: precipitated share of solid Cr(III) + N(0, {ms.PHASE_ABS_SD}), clipped",
        "pH_final": "true final pH + N(0, 0.05)",
        "study bias": "exp(N(0, 0.05)), shared by Cr(VI) and total-Cr analyses within a study",
    }


MEASUREMENT = _measurement_spec()



# ---------------------------------------------------------------------------
# PHREEQC look-ups
# ---------------------------------------------------------------------------
class _Sol:
    pH = np.round(np.arange(1.5, 11.01, 0.25), 2)
    I = np.array([0.001, 0.01, 0.1])
    _logS = None

    @classmethod
    def logS(cls):
        if cls._logS is None:
            f = ch.DB_DIR.parent / "results" / "cr3_solubility_grid.npz"
            if f.exists():
                cls._logS = np.load(f)["logS"]
            else:
                S = np.array([[ch.cr3_solubility_si(p, ionic_strength=i)["S_mgL"] for i in cls.I] for p in cls.pH])
                S = np.where(np.isfinite(S), S, 1e5)   # 'undersaturated_to_cmax' -> no constraint
                S = np.where(np.isnan(S), 1e5, S)
                cls._logS = np.log10(S)
                f.parent.mkdir(exist_ok=True)
                np.savez_compressed(f, logS=cls._logS, pH=cls.pH, I=cls.I)
        return cls._logS


def cr3_solubility(pH, I=0.01):
    logS = _Sol.logS()
    li = np.interp(np.log10(I), np.log10(_Sol.I), np.arange(len(_Sol.I)))
    j = int(min(np.floor(li), len(_Sol.I) - 2)); w = li - j
    col = (1 - w) * logS[:, j] + w * logS[:, j + 1]
    return 10 ** np.interp(pH, _Sol.pH, col)


def _spec(pH, C, I, T):
    f = ch.grid()(pH, C, I, T)[0]
    return f[0] + f[2]  # HCrO4- + Cr2O7-2: the high-affinity, low-pH anions (per Cr atom)


# ---------------------------------------------------------------------------
# Materials, hyper-priors and scenarios
# ---------------------------------------------------------------------------
@dataclass
class Material:
    material_id: str
    family: str
    SBET: float
    Vp: float
    dp: float
    pHpzc: float
    qmax: float
    K0: float
    ka: float
    R0: float
    kB: float
    kC: float
    L0: float
    kl: float
    kD: float
    pH50_red: float
    Ea_red: float
    pH50_cr3: float
    krel: float
    q3max: float
    K3: float
    k3: float
    nH: float
    beta: float
    kbuf: float
    kp: float = 0.05
    Ea_ads: float = 15.0
    w_elec: float = 0.8   # electrostatic switch width (pH units)
    f_base: float = 0.05  # affinity floor of CrO4-2-dominated solutions
    e_base: float = 0.1   # affinity floor above pHpzc
    scenario: str = "E"


def corpus_hyperpriors(rng) -> dict:
    """Family-level means drawn per corpus: no family is reductive by construction."""
    return {f: {"q_f": np.exp(rng.normal(np.log(35), 0.3)), "mu_R": rng.uniform(-1.3, 0.3),
                "mu_L": rng.uniform(-2.5, -0.8)} for f in FAMILIES}


def random_material(rng, family, idx, hyper, scenario="E") -> Material:
    (s_lo, s_hi), (z_lo, z_hi) = FAMILIES[family]
    h = hyper[family]
    S = rng.uniform(s_lo, s_hi)
    qmax = h["q_f"] * S / 1000 * np.exp(rng.normal(0, 0.35))
    m = Material(
        material_id=f"{family[:3].upper()}{idx:04d}", family=family, SBET=S,
        Vp=S / 1000 * rng.uniform(0.35, 0.75), dp=rng.uniform(1.6, 6.0), pHpzc=rng.uniform(z_lo, z_hi),
        qmax=qmax, K0=10 ** rng.uniform(-1.3, 0), ka=10 ** rng.uniform(-2.3, -1.3),
        R0=10 ** (h["mu_R"] + rng.normal(0, 0.3)), kB=10 ** rng.uniform(-4.0, -2.3),
        kC=10 ** rng.uniform(-3.5, -1.8), L0=10 ** (h["mu_L"] + rng.normal(0, 0.3)),
        kl=10 ** rng.uniform(-3.0, -1.5), kD=10 ** rng.uniform(-2.0, 0.0),
        pH50_red=rng.uniform(2.0, 4.0), Ea_red=rng.uniform(25, 50), pH50_cr3=rng.uniform(3.5, 5.0),
        krel=10 ** rng.uniform(-3.0, -1.5), q3max=0.5 * qmax * np.exp(rng.normal(0, 0.3)),
        K3=10 ** rng.uniform(-1.0, 0.5), k3=10 ** rng.uniform(-2.5, -1.0),
        nH=rng.uniform(3, 5), beta=10 ** rng.uniform(0, 1), kbuf=10 ** rng.uniform(-4.5, -3.0),
        scenario=scenario)
    return apply_scenario(m, scenario)


def apply_scenario(m: Material, scenario: str) -> Material:
    d = asdict(m)
    d["scenario"] = scenario
    if scenario in ("A", "C", "D"):
        d["kB"] = 0.0
    if scenario in ("A", "B", "D"):
        d["kC"] = 0.0
    if scenario in ("A", "B", "C"):
        d["kD"] = 0.0
    return Material(**d)


# ---------------------------------------------------------------------------
# Kinetic model
# ---------------------------------------------------------------------------
def _tables(pH_fixed, I, T, drift):
    """Per-batch speciation/solubility tables so the ODE right-hand side avoids grid look-ups."""
    g = ch.grid()
    lc = g.axes[1]
    C = 10 ** lc
    pHax = _Sol.pH
    logS_pH = np.log10(cr3_solubility(pHax, I))
    if drift:
        pg = g.axes[0]
        F = g(pg[:, None], C[None, :], I, T)
        F = F[..., 0] + F[..., 2]
        return lc, pg, F, pHax, logS_pH
    f = g(np.full_like(C, pH_fixed), C, I, T)
    return lc, None, f[:, 0] + f[:, 2], None, float(np.interp(pH_fixed, pHax, logS_pH))


def _rhs_factory(m: Material, D, T, I, pH_fixed):
    import math
    from bisect import bisect_right
    R = 8.314e-3
    arr_r = math.exp(-m.Ea_red / R * (1 / T - 1 / 298.15))
    arr_a = math.exp(-m.Ea_ads / R * (1 / T - 1 / 298.15))
    drift = m.scenario == "E"
    k_mmol = 1.0 / MW
    lc, pg, F, pHax, logS = _tables(pH_fixed, I, T, drift)
    lc_l = lc.tolist(); nC = len(lc_l)
    if drift:
        pg_l = pg.tolist(); npg = len(pg_l); F_l = F.tolist(); pHax_l = pHax.tolist(); logS_l = logS.tolist()
    else:
        F_l = F.tolist(); S_fixed = 10 ** logS
    qmax, K0, ka, kB, kC, kD = m.qmax, m.K0, m.ka, m.kB, m.kC, m.kD
    R0, pzc, pH50r, pH50c = m.R0, m.pHpzc, m.pH50_red, m.pH50_cr3
    krel, q3max, K3, k3, kp, kl, L0 = m.krel, m.q3max, m.K3, m.k3, m.kp, m.kl, m.L0
    nH, beta, kbuf = m.nH, m.beta, m.kbuf
    w_el, f_b, e_b = m.w_elec, m.f_base, m.e_base

    def interp1(xs, ys, x, n):
        if x <= xs[0]:
            return ys[0]
        if x >= xs[-1]:
            return ys[-1]
        i = bisect_right(xs, x) - 1
        w = (x - xs[i]) / (xs[i + 1] - xs[i])
        return ys[i] * (1 - w) + ys[i + 1] * w

    def spec(pH, Cv):
        x = math.log10(max(Cv, 1e-3))
        if not drift:
            return interp1(lc_l, F_l, x, nC)
        p = min(max(pH, pg_l[0]), pg_l[-1])
        i = min(bisect_right(pg_l, p) - 1, npg - 2)
        w = (p - pg_l[i]) / (pg_l[i + 1] - pg_l[i])
        return (1 - w) * interp1(lc_l, F_l[i], x, nC) + w * interp1(lc_l, F_l[i + 1], x, nC)

    def rhs(y, _):
        Cv, q6, q3s, q3r, C3, P3, Rs, L, Lrel, pH = y
        Cv = Cv if Cv > 0 else 0.0; q6 = q6 if q6 > 0 else 0.0; q3s = q3s if q3s > 0 else 0.0
        C3 = C3 if C3 > 0 else 0.0; P3 = P3 if P3 > 0 else 0.0; Rs = Rs if Rs > 0 else 0.0
        L = L if L > 0 else 0.0
        if not drift:
            pH = pH_fixed
        f = spec(pH, Cv)
        ze = (pH - pzc) / w_el
        elec = 1.0 / (1.0 + math.exp(ze)) if ze < 700 else 0.0
        K = K0 * (f_b + f) * (e_b + elec)
        qeq = qmax * K * Cv / (1 + K * Cv)
        r_ads = ka * arr_a * (qeq - q6)
        h = arr_r / (1.0 + 10 ** (pH - pH50r))
        sfrac = Rs / R0 if R0 > 0 else 0.0
        r_B = kB * h * sfrac * Cv
        r_C = kC * h * sfrac * q6
        r_D = kD * h * L * Cv
        z = -(pH - pH50c) / 0.5
        th3 = 1.0 / (1.0 + math.exp(z)) if z < 700 else 0.0
        r_rel = krel * (1 - th3) * q3s
        q3eq = q3max * th3 * K3 * C3 / (1 + K3 * C3)
        r_3 = k3 * (q3eq - q3r)
        S = 10 ** interp1(pHax_l, logS_l, pH, 0) if drift else S_fixed
        if C3 > S:
            r_p = kp * (C3 - S)
        elif P3 > 0:
            r_p = max(kp * (C3 - S), -kp * P3)
        else:
            r_p = 0.0
        r_leach = kl * (D * L0 - Lrel)
        dCv = -D * (r_ads + r_B) - r_D
        dq6 = r_ads - r_C
        dq3s = th3 * r_B + r_C - r_rel
        dC3 = D * (1 - th3) * r_B + D * r_rel + r_D - D * r_3 - r_p
        dRs = -3 * k_mmol * (r_B + r_C)
        dL = r_leach - 3 * k_mmol * r_D
        if drift:
            red_mmol = k_mmol * (D * (r_B + r_C) + r_D)
            beta_tot = beta + 2303.0 * 10 ** (-pH)          # carbon + free-proton buffering, mmol/(L pH)
            dpH = (nH * red_mmol - 3 * k_mmol * r_p + kbuf * beta * (pzc - pH)) / beta_tot
        else:
            dpH = 0.0
        return [dCv, dq6, dq3s, r_3, dC3, r_p, dRs, dL, r_leach, dpH]

    return rhs


def run_batch(m: Material, C0, pH0, D, t_points, T=298.15, I=0.01):
    """Integrate one batch; returns a dict of arrays at ``t_points`` (min)."""
    from scipy.integrate import odeint
    t_points = np.atleast_1d(np.asarray(t_points, float))
    y0 = [C0, 0, 0, 0, 0, 0, m.R0, 0, 0, pH0]
    tt = np.concatenate([[0.0], t_points])
    Y, info = odeint(_rhs_factory(m, D, T, I, pH0), y0, tt, rtol=1e-6, atol=1e-8, full_output=True, mxstep=20000)
    if info["message"] != "Integration successful.":
        raise RuntimeError(info["message"])
    Cv, q6, q3s, q3r, C3, P3, Rs, L, Lrel, pH = Y[1:].T
    Cv, C3, P3 = np.maximum(Cv, 0), np.maximum(C3, 0), np.maximum(P3, 0)
    removed = np.maximum(C0 - Cv, 1e-12)
    reduced = D * (q3s + q3r) + C3 + P3
    solid = D * (q6 + q3s + q3r) + P3
    return {"t": t_points, "Ce_VI_true": Cv, "q_adsVI_true": q6, "Cr3_surface_true": D * (q3s + q3r),
            "Cr3_aq_true": C3, "Cr3_precip_true": P3, "reduced_mgL": reduced,
            "Rs": Rs, "L": L, "Lrel": Lrel, "pH_true": pH if m.scenario == "E" else np.full_like(Cv, pH0),
            "share_ads": D * q6 / removed, "share_red": reduced / removed,
            "share_precip": P3 / removed, "share_cr3_aq": C3 / removed,
            "xps_cr3_true": np.where(solid > 0, (solid - D * q6) / np.maximum(solid, 1e-12), 0.0),
            "solid_cr_mgg_true": solid / D,
            "precip_frac_of_cr3_true": np.where(solid - D * q6 > 1e-12, P3 / np.maximum(solid - D * q6, 1e-12), 0.0)}


# ---------------------------------------------------------------------------
# Virtual studies and corpora
# ---------------------------------------------------------------------------
EVIDENCE_MIX = {"C": 0.40, "B": 0.30, "A1": 0.20, "A2": 0.10}
DOSE_SERIES = (0.5, 1.0, 2.0, 4.0, 6.0, 8.0)  # g/L at C0 = 20 mg/L, the study's optimum pH


def virtual_study(rng, study_idx, family, hyper, scenario="E", n_iso=8, n_edge=6, n_kin=8, prefix="VS"):
    """One virtual study. Evidence (observation level):
    C  residual Cr(VI) only;
    B  + total aqueous Cr on isotherm, pH-edge and dose series;
    A1 B + solid oxidation state (XPS or XANES, half each) on 3 solids + final pH + total Cr on kinetics;
    A2 A1 (XANES) + solid Cr inventory and phase/coordination evidence on the same 3 solids.
    """
    m = random_material(rng, family, study_idx, hyper, scenario)
    sid = f"{prefix}{study_idx:03d}"
    ev = rng.choice(list(EVIDENCE_MIX), p=list(EVIDENCE_MIX.values()))
    tech = None if ev in ("C", "B") else ("XANES" if ev == "A2" else str(rng.choice(["XPS", "XANES"])))
    enrich, beam = ms.draw_technique_bias(rng, tech) if tech else (1.0, 0.0)
    bias = np.exp(rng.normal(0, 0.05))
    I = float(rng.choice([0.001, 0.005, 0.01, 0.05, 0.1]))
    T = float(rng.choice([293.15, 298.15, 303.15, 308.15, 313.15]))
    D = float(rng.choice([0.5, 1.0, 2.0, 4.0]))
    pH_opt = float(rng.uniform(1.8, 3.5))
    t_eq = float(rng.choice([360, 720, 1440]))
    full = ev in ("A1", "A2")
    reports_I = True if full else bool(rng.random() < {"B": 0.5, "C": 0.3}[ev])
    reports_final_pH = full or (ev == "B" and rng.random() < 0.5)
    rows = []

    def add(series, C0, pH0, t_arr, total_meas, solid_meas, Dose=None):
        Dz = D if Dose is None else Dose
        tr = run_batch(m, C0, pH0, Dz, t_arr, T, I)
        for k in range(len(tr["t"])):
            one = {key: (v[k] if np.ndim(v) else v) for key, v in tr.items()}
            cv = one["Ce_VI_true"] * bias * (1 + rng.normal(0, 0.03))
            below = cv < 0.01
            row = {"study_id": sid, "material_id": m.material_id, "material_family": family,
                   "scenario": scenario, "study_evidence": ev, "solid_technique": tech or "none",
                   "series": series, "C0": C0, "pH": pH0,
                   "D": Dz, "t": float(one["t"]), "T": T, "I": I if reports_I else np.nan, "I_true": I,
                   "SBET": m.SBET, "Vp": m.Vp, "dp": m.dp, "pHpzc": m.pHpzc,
                   "Ce_VI": 0.005 if below else cv, "Ce_VI_below_LOD": below}
            if total_meas:
                ct = (one["Ce_VI_true"] + one["Cr3_aq_true"]) * bias * (1 + rng.normal(0, 0.03))
                row["Ce_total"] = max(ct, 0.01)
            if solid_meas:
                f = float(ms.measure_oxidation_state(rng, one["xps_cr3_true"], tech, enrich, beam))
                row["xps_cr3_frac" if tech == "XPS" else "xanes_cr3_frac"] = f
                if ev == "A2":
                    row["solid_cr_mgg"] = float(ms.measure_inventory(rng, one["solid_cr_mgg_true"]))
                    row["cr3_precip_frac"] = float(ms.measure_phase(rng, one["precip_frac_of_cr3_true"]))
            row["pH_final"] = one["pH_true"] + rng.normal(0, 0.05) if reports_final_pH else np.nan
            row.update({k2: one[k2] for k2 in ("Ce_VI_true", "q_adsVI_true", "Cr3_surface_true", "Cr3_aq_true",
                                               "Cr3_precip_true", "reduced_mgL", "share_ads", "share_red",
                                               "share_precip", "share_cr3_aq", "pH_true", "xps_cr3_true",
                                               "solid_cr_mgg_true", "precip_frac_of_cr3_true")})
            rows.append(row)

    iso_C0 = np.geomspace(5, 300, n_iso)
    for j, C0 in enumerate(iso_C0):
        add("isotherm", float(C0), pH_opt, [t_eq], ev != "C", full and j in (0, n_iso - 1))
    for j, pH in enumerate(np.linspace(2.0, 8.0, n_edge)):
        add("pH_edge", 50.0, float(pH), [t_eq], ev != "C", full and j == 0)
    add("kinetic", 50.0, pH_opt, np.geomspace(5, t_eq, n_kin), full, False)
    for Dose in DOSE_SERIES:
        add("dose", 20.0, pH_opt, [t_eq], ev != "C", False, Dose=Dose)
    df = pd.DataFrame(rows)
    meta = asdict(m) | {"study_id": sid, "I": I, "T": T, "D": D, "pH_opt": pH_opt, "bias": bias,
                        "t_eq": t_eq, "study_evidence": ev, "solid_technique": tech or "none",
                        "xps_enrich": enrich, "beam_reduction": beam}
    return df, meta


def virtual_corpus(n_studies=30, seed=2026, scenario="E", prefix="VS"):
    rng = np.random.default_rng(seed)
    hyper = corpus_hyperpriors(rng)
    fams = list(FAMILIES)
    dfs, mats = [], []
    for i in range(n_studies):
        df, meta = virtual_study(rng, i, fams[i % len(fams)], hyper, scenario, prefix=prefix)
        dfs.append(df)
        mats.append(meta)
    mats = pd.DataFrame(mats)
    for f, h in hyper.items():
        mats.loc[mats.family == f, "mu_R_family"] = h["mu_R"]
    return pd.concat(dfs, ignore_index=True), mats


STUDY_DESIGN = {
    "SBET, pHpzc": ("m2/g, -", "uniform within family ranges: " + "; ".join(
        f"{f}: SBET {a}-{b}, pHpzc {c}-{d}" for f, ((a, b), (c, d)) in FAMILIES.items())),
    "Vp, dp": ("cm3/g, nm", "Vp = SBET/1000 U(0.35, 0.75); dp ~ U(1.6, 6.0)"),
    "affinity": ("-", "K = K0 (f_base + f[HCrO4- + Cr2O7-2]) (e_base + 1/(1 + exp((pH - pHpzc)/w_elec))); "
                      "f_base 0.05, e_base 0.1, w_elec 0.8 (fixed)"),
    "evidence class of a study": ("-", ", ".join(f"{k} {v:.2f}" for k, v in EVIDENCE_MIX.items())
                                      + "; B: total aqueous Cr on isotherm, pH-edge and dose series; A1: B + solid oxidation "
                                        "state on 3 solids (XPS or XANES, probability 0.5 each) + final pH + total Cr on kinetics; "
                                        "A2: A1 with XANES + solid Cr inventory + phase evidence on the same 3 solids"),
    "I": ("mol/L", "one of 0.001, 0.005, 0.01, 0.05, 0.1; reported always (A1, A2), 50% (B), 30% (C)"),
    "final pH reported": ("-", "always (A1, A2), 50% (B), never (C)"),
    "corpus composition": ("-", "30 studies per corpus, families assigned in rotation (5 studies per family); one material per study"),
    "family hyper-priors": ("-", "redrawn per corpus and family: q_f ~ LN(ln 35, 0.3), mu_R,f ~ U(-1.3, 0.3), mu_L,f ~ U(-2.5, -0.8)"),
    "scenario switches": ("-", "A: kB = kC = kD = 0; B: kC = kD = 0; C: kB = kD = 0; D: kB = kC = 0; E: all pathways"),
    "T": ("K", "one of 293.15-313.15 in 5 K steps"),
    "D": ("g/L", "one of 0.5, 1, 2, 4 (dose series: 0.5, 1, 2, 4, 6, 8 at C0 = 20 mg/L)"),
    "pH_opt": ("-", "U(1.8, 3.5): pH of isotherm, kinetic and dose series"),
    "t_eq": ("min", "one of 360, 720, 1440"),
    "series": ("-", "isotherm 8 x C0 5-300 mg/L (log-spaced); pH edge 6 x pH 2-8 at 50 mg/L; kinetics 8 x t 5-t_eq; dose 6"),
}


def parameter_table() -> pd.DataFrame:
    rows = [{"parameter": k, "unit": u, "distribution": d} for k, (u, d) in PARAMETERS.items()]
    rows += [{"parameter": k, "unit": u, "distribution": d} for k, (u, d) in STUDY_DESIGN.items()]
    rows += [{"parameter": k, "unit": "-", "distribution": v} for k, v in MEASUREMENT.items()]
    rows += [{"parameter": "numerical", "unit": "-",
              "distribution": "LSODA (scipy odeint), rtol 1e-6, atol 1e-8, mxstep 20000; y0 = [C0, 0, 0, 0, 0, 0, R0, 0, 0, pH0]"}]
    return pd.DataFrame(rows)
