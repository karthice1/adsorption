"""Mechanism-identifiability engine (CASSIA Layer 4).

Question: given a measurement design, how well can the share of Cr(VI) disappearance
caused by reduction be inferred, and which extra measurements make it identifiable?

Designs (cumulative):
  D1  residual Cr(VI) at equilibrium (isotherm + pH edge)            -> evidence class C
  D2  D1 + Cr(VI) kinetics                                           -> class C
  D3  D2 + total dissolved Cr at every point                         -> class B
  D4  D3 + XPS Cr(III) fraction on 3 solids + final pH at every point -> class A

Tools:
  * local (Fisher-information) identifiability with a weakly informative prior; the
    Laplace posterior SD of the reduction share is compared with its prior SD;
  * detection tests that need no mechanistic model fit (lack of isotherm fit, dissolved
    Cr(III) mass balance, XPS Cr(III));
  * profile likelihood of the reductant capacity with re-optimisation of all others;
  * first-order bias of the inferred share from misreported pH or ionic strength.
"""
from __future__ import annotations

from dataclasses import asdict

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import least_squares

from . import fitting as fit
from . import measurement as ms
from . import simulate as sim

LOG_PARAMS = ["qmax", "K0", "ka", "R0", "kB", "kC", "L0", "kD", "krel", "q3max", "K3", "k3"]
LIN_PARAMS = ["pH50_red", "pH50_cr3"]
PARAMS = LOG_PARAMS + LIN_PARAMS
# prior SD: log10 units for LOG_PARAMS, pH units for LIN_PARAMS (wide, weakly informative)
PRIOR_SD = {**{p: 0.5 for p in LOG_PARAMS}, "pH50_red": 0.6, "pH50_cr3": 0.45}
DESIGNS = ("D1", "D2", "D3", "D4")                      # v0.2 designs (D4 = unbiased solid oxidation state)
DESIGNS_V3 = ("D1", "D2", "D3", "D4x", "D4n", "D5")    # v0.3 designs with technique-specific measurement models
TECH_OF = {"D4x": "XPS", "D4n": "XANES", "D5": "XANES"}
TECH_NOMINAL = {"XPS": (1.0, 0.075), "XANES": (1.0, 0.015)}   # (enrichment, beam reduction) at prior means
DESIGN_LABEL = {"D1": "Cr(VI) at equilibrium (class C)", "D2": "+ Cr(VI) kinetics (class C)",
                "D4x": "+ XPS oxidation state, final pH (class A1)", "D4n": "+ XANES oxidation state, final pH (class A1)",
                "D5": "+ XANES, solid Cr inventory, phase evidence (class A2)",
                "D3": "+ total dissolved Cr (class B)", "D4": "+ XPS and final pH (class A)"}


def _theta(m: sim.Material):
    d = asdict(m)
    return np.array([np.log10(max(d[p], 1e-12)) for p in LOG_PARAMS] + [d[p] for p in LIN_PARAMS])


def _with(m: sim.Material, th):
    d = asdict(m)
    for i, p in enumerate(LOG_PARAMS):
        if d[p] > 0:  # pathways switched off by the scenario stay off
            d[p] = 10 ** th[i]
    for j, p in enumerate(LIN_PARAMS):
        d[p] = th[len(LOG_PARAMS) + j]
    return sim.Material(**d)


def active_mask(m: sim.Material):
    d = asdict(m)
    return np.array([d[p] > 0 for p in LOG_PARAMS] + [True] * len(LIN_PARAMS))


def design_points(meta: dict, n_iso=8, n_edge=6, n_kin=8):
    iso = [(float(c), meta["pH_opt"]) for c in np.geomspace(5, 300, n_iso)]
    edge = [(50.0, float(p)) for p in np.linspace(2.0, 8.0, n_edge)]
    kin_t = np.geomspace(5, meta["t_eq"], n_kin)
    return iso, edge, kin_t


def observe(m: sim.Material, meta: dict, pH_shift=0.0, I_factor=1.0):
    """Noise-free observables for every design element + true reduction share (isotherm, mass-weighted)."""
    D, T, I, teq = meta["D"], meta["T"], meta["I"] * I_factor, meta["t_eq"]
    iso, edge, kin_t = design_points(meta)
    cv, ct, ph, xps, inv, pfr, red, rem, prec = [], [], [], [], [], [], 0.0, 0.0, 0.0
    for k, (C0, pH0) in enumerate(iso + edge):
        tr = sim.run_batch(m, C0, pH0 + pH_shift, D, [teq], T, I)
        cv.append(tr["Ce_VI_true"][0]); ct.append(tr["Ce_VI_true"][0] + tr["Cr3_aq_true"][0])
        ph.append(tr["pH_true"][0])
        if k in (0, len(iso) - 1, len(iso)):
            xps.append(tr["xps_cr3_true"][0])
            inv.append(tr["solid_cr_mgg_true"][0]); pfr.append(tr["precip_frac_of_cr3_true"][0])
        if k < len(iso):
            red += tr["reduced_mgL"][0]; rem += C0 - tr["Ce_VI_true"][0]; prec += tr["Cr3_precip_true"][0]
    tk = sim.run_batch(m, 50.0, meta["pH_opt"] + pH_shift, D, kin_t, T, I)
    out = {"cv_eq": np.array(cv), "cv_kin": tk["Ce_VI_true"], "ct_eq": np.array(ct),
           "ct_kin": tk["Ce_VI_true"] + tk["Cr3_aq_true"], "xps": np.array(xps), "pH": np.array(ph),
           "inv": np.array(inv), "pfrac": np.array(pfr)}
    out["share_red"] = red / max(rem, 1e-12)
    out["share_precip"] = prec / max(rem, 1e-12)
    return out


def stack(obs: dict, design: str, log_bias=0.0, tech=None):
    """Observation vector and SDs for one design (concentrations carry the study bias).

    ``tech`` = (enrichment, beam reduction) of the solid oxidation-state technique; defaults to
    the technique's prior mean. Legacy design D4 is the v0.2 unbiased solid oxidation state.
    """
    b = 10 ** log_bias
    parts, sds = [b * obs["cv_eq"]], [0.03 * b * obs["cv_eq"] + 0.01]
    if design != "D1":
        parts.append(b * obs["cv_kin"]); sds.append(0.03 * b * obs["cv_kin"] + 0.01)
    if design not in ("D1", "D2"):
        parts += [b * obs["ct_eq"], b * obs["ct_kin"]]
        sds += [0.03 * b * obs["ct_eq"] + 0.02, 0.03 * b * obs["ct_kin"] + 0.02]
    if design == "D4":
        parts += [obs["xps"], obs["pH"]]
        sds += [np.full(len(obs["xps"]), 0.05), np.full(len(obs["pH"]), 0.05)]
    elif design in TECH_OF:
        t = TECH_OF[design]
        e, bm = tech if tech is not None else TECH_NOMINAL[t]
        parts += [ms.oxidation_state_bias(obs["xps"], e, bm), obs["pH"]]
        sds += [np.full(len(obs["xps"]), ms.TECHNIQUES[t]["noise"]), np.full(len(obs["pH"]), 0.05)]
        if design == "D5":
            parts += [obs["inv"], obs["pfrac"]]
            sds += [ms.INVENTORY_REL_SD * obs["inv"] + 0.01, np.full(len(obs["pfrac"]), ms.PHASE_ABS_SD)]
    return np.concatenate(parts), np.concatenate(sds)


def fisher_study(m: sim.Material, meta: dict, h=0.02):
    """Jacobians of all observables and of the reduction share w.r.t. (theta, log study bias)."""
    th0 = _theta(m)
    act = active_mask(m)
    base = observe(m, meta)
    jac, gshare = {}, []
    cols = []
    for i in range(len(th0)):
        if not act[i]:
            continue
        th = th0.copy(); th[i] += h
        o = observe(_with(m, th), meta)
        cols.append((i, o))
        gshare.append((o["share_red"] - base["share_red"]) / h)
    out = {"base": base, "active": [PARAMS[i] for i, _ in cols], "g_share": np.array(gshare + [0.0])}
    for d in DESIGNS:
        y0, sd = stack(base, d)
        J = np.column_stack([(stack(o, d)[0] - y0) / h for _, o in cols] +
                            [y0 * np.log(10)])  # d y / d log10 bias
        out[d] = (J, sd)
    return out


def posterior_share_sd(fs: dict, design: str):
    J, sd = fs[design]
    W = 1.0 / sd ** 2
    prior = np.array([PRIOR_SD[p] for p in fs["active"]] + [0.05 / np.log(10)])
    S0inv = np.diag(1.0 / prior ** 2)
    Spost = np.linalg.inv(J.T @ (W[:, None] * J) + S0inv)
    g = fs["g_share"]
    post = float(np.sqrt(max(g @ Spost @ g, 0)))
    pri = float(np.sqrt(g @ np.diag(prior ** 2) @ g))
    ev = np.linalg.eigvalsh(J.T @ (W[:, None] * J))
    return {"design": design, "post_sd_share": post, "prior_sd_share": pri, "sd_ratio": post / pri if pri > 0 else np.nan,
            "fim_cond": float(ev.max() / max(ev.min(), 1e-300)), "fim_rank_1e-6": int((ev > ev.max() * 1e-6).sum()),
            "n_params": J.shape[1]}


def metadata_bias(m: sim.Material, meta: dict, fs: dict, design="D3", pH_shift=0.5, I_factor=10.0):
    """First-order bias in the inferred reduction share when pH or I is misreported.

    The data are generated at the true conditions; the analyst models them at pH + shift
    (or I * factor). Linearised estimator bias: d_theta = S_post J' W (y_true - y_assumed).
    """
    J, sd = fs[design]
    W = 1.0 / sd ** 2
    prior = np.array([PRIOR_SD[p] for p in fs["active"]] + [0.05 / np.log(10)])
    Spost = np.linalg.inv(J.T @ (W[:, None] * J) + np.diag(1.0 / prior ** 2))
    y_true, _ = stack(fs["base"], design)
    res = {}
    for name, kw in (("pH", {"pH_shift": pH_shift}), ("I", {"I_factor": I_factor})):
        y_ass, _ = stack(observe(m, meta, **kw), design)
        dth = Spost @ J.T @ (W * (y_true - y_ass))
        res[f"bias_share_{name}"] = float(fs["g_share"] @ dth)
    return res


def detection_tests(rng, fs: dict, meta: dict, n_rep=20):
    """Power of three model-free diagnostics, using simulated noisy observations."""
    o = fs["base"]
    iso_n = 8
    C0_iso = np.geomspace(5, 300, iso_n)
    hits = {"lack_of_fit": 0, "dissolved_cr3": 0, "xps_cr3": 0}
    for _ in range(n_rep):
        b = np.exp(rng.normal(0, 0.05))
        cv = np.maximum(o["cv_eq"] * b * (1 + rng.normal(0, 0.03, o["cv_eq"].size)), 0.005)
        ct = np.maximum(o["ct_eq"] * b * (1 + rng.normal(0, 0.03, o["ct_eq"].size)), 0.01)
        q = (C0_iso - cv[:iso_n]) / meta["D"]
        fits = fit.fit_isotherms(cv[:iso_n], q, T=meta["T"])
        if fits:
            tab = fit.selection_table(fits)
            r2 = tab["R2"].max()
            runs = tab.loc[0, "runs_p"]
            hits["lack_of_fit"] += int(r2 < 0.95 or (runs == runs and runs < 0.05))
        diff = ct - cv
        sdd = np.sqrt((0.03 * ct) ** 2 + (0.03 * cv) ** 2 + 0.02 ** 2)
        z = diff.sum() / np.sqrt((sdd ** 2).sum())
        hits["dissolved_cr3"] += int(z > stats.norm.ppf(0.99))
        x = np.clip(o["xps"] + rng.normal(0, 0.05, o["xps"].size), 0, 1)
        hits["xps_cr3"] += int(x.max() > 0.15)
    return {k: v / n_rep for k, v in hits.items()}


def profile_likelihood(m: sim.Material, meta: dict, design: str, grid_log_R0, seed=0, max_nfev=60):
    """Profile of log10 R0 (reductant capacity) with all other active parameters re-optimised.

    Data: one noisy realisation of ``design``. Returns -2 log L profile and the implied
    reduction share along it.
    """
    rng = np.random.default_rng(seed)
    base = observe(m, meta)
    y, sd = stack(base, design)
    y = y * (1 + rng.normal(0, 0.03, y.size))
    th_true = _theta(m)
    act = active_mask(m)
    iR = LOG_PARAMS.index("R0")
    free = [i for i in range(len(th_true)) if act[i] and i != iR]
    prior = np.array([PRIOR_SD[PARAMS[i]] for i in free])
    rows = []
    x_prev = th_true[free] + rng.normal(0, 0.1, len(free))
    for lr in grid_log_R0:
        def resid(x):
            th = th_true.copy(); th[free] = x; th[iR] = lr
            try:
                yy, _ = stack(observe(_with(m, th), meta), design)
            except RuntimeError:
                return np.full(y.size + len(free), 1e3)
            return np.concatenate([(y - yy) / sd, (x - th_true[free]) / prior])
        # two starts (previous profile point, true parameters); keep the better optimum
        best = None
        for x0 in (x_prev, th_true[free]):
            r = least_squares(resid, x0, max_nfev=max_nfev, x_scale=0.3)
            if best is None or r.cost < best.cost:
                best = r
        th = th_true.copy(); th[free] = best.x; th[iR] = lr
        share = observe(_with(m, th), meta)["share_red"]
        rows.append({"design": design, "log10_R0": lr, "m2logL": float(2 * best.cost), "share_red": share,
                     "n_obs": int(y.size), "nfev": int(best.nfev), "converged": bool(best.status > 0)})
        x_prev = best.x
    df = pd.DataFrame(rows)
    df["delta"] = df["m2logL"] - df["m2logL"].min()
    return df


def _empirical_predict(x, meta, base, design):
    """Adsorption-only empirical model: pH-dependent Sips equilibrium + PFO kinetics.

    x = [log10 qmax, n, a0, a1, p0, log10 w, log10 k1, log10 bias, delta_pH, log10 qe_kin, nA]
    log10 K(pH) = a0 + a1 / (1 + exp((pH - p0)/w)); q = qmax (K C)^n / (1 + (K C)^n).
    Kinetics are fitted as analysts usually do, separately: Avrami-type
    q(t) = qe_kin (1 - exp(-(k1 t)^nA)) with its own qe_kin. Dissolved and solid Cr(III) are
    zero by construction; final pH relaxes by a free fraction delta toward pHpzc.
    """
    from scipy.optimize import brentq
    qmax, n, a0, a1, p0, w, k1, b, dlt = 10 ** x[0], x[1], x[2], x[3], x[4], 10 ** x[5], 10 ** x[6], 10 ** x[7], x[8]
    qek, nA = 10 ** x[9], x[10]
    D = meta["D"]
    iso, edge, kin_t = design_points(meta)

    def qeq(C, pH):
        K = 10 ** (a0 + a1 / (1 + np.exp(np.clip((pH - p0) / w, -50, 50))))
        u = (K * max(C, 0)) ** n
        return qmax * u / (1 + u)

    def ce(C0, pH):
        f = lambda C: C0 - C - D * qeq(C, pH)
        return C0 if f(C0) >= 0 else brentq(f, 0.0, C0, xtol=1e-10)

    cv = np.array([ce(C0, pH) for C0, pH in iso + edge])
    C0k = 50.0
    cvk = np.maximum(C0k - D * qek * (1 - np.exp(-(k1 * kin_t) ** nA)), 0.0)
    parts = [b * cv]
    if design in ("D2", "D3", "D4"):
        parts.append(b * cvk)
    if design in ("D3", "D4"):
        parts += [b * cv, b * cvk]
    if design == "D4":
        pH0 = np.array([p for _, p in iso + edge])
        parts += [np.zeros(len(base["xps"])), pH0 + dlt * (meta["pHpzc"] - pH0)]
    return np.concatenate(parts)


def equivalence_test(m: sim.Material, meta: dict, design: str, sd_rel=0.03, seed=0, n_starts=4):
    """Can an adsorption-only model reproduce data generated with reduction?

    One noisy realisation of ``design`` (relative analytical SD ``sd_rel``) is generated
    from the true multi-mechanism model. The empirical adsorption-only model above
    (9 free parameters, multi-start) is fitted by weighted least squares; p_A > 0.05 from
    the chi-square goodness-of-fit test means the data cannot distinguish 'adsorption only'
    from the true mixture of mechanisms.
    """
    rng = np.random.default_rng(seed)
    base = observe(m, meta)
    y, _ = stack(base, design)
    nx = len(base["xps"]) + len(base["pH"]) if design == "D4" else 0
    nc_eq = len(base["cv_eq"]) + len(base["cv_kin"]) - 1  # last kinetic point
    nc = y.size - nx
    floor = np.concatenate([np.full(nc, 0.01), np.full(nx, 0.05)])
    sd = np.concatenate([sd_rel * y[:nc], np.zeros(nx)]) + floor
    bias = np.exp(rng.normal(0, 0.05))
    yobs = y.copy()
    yobs[:nc] = y[:nc] * bias * (1 + rng.normal(0, sd_rel, nc))
    yobs[nc:] = y[nc:] + rng.normal(0, 0.05, nx)
    lo = np.array([-1.0, 0.2, -4.0, -2.0, 1.0, np.log10(0.1), -4.0, -0.15, 0.0, -2.0, 0.2])
    hi = np.array([3.0, 3.0, 2.0, 6.0, 10.0, np.log10(5.0), 0.0, 0.15, 1.0, 3.0, 2.5])

    def resid(x):
        return (yobs - _empirical_predict(x, meta, base, design)) / sd

    best = None
    for k in range(n_starts):
        x0 = lo + (hi - lo) * rng.uniform(0.2, 0.8, lo.size)
        x0[0] = np.log10(max(m.qmax, 0.1)); x0[7] = 0.0
        x0[9] = np.log10(max((50.0 - yobs[nc_eq]) / meta["D"], 0.05)) if design != "D1" else 0.0
        r = least_squares(resid, np.clip(x0, lo + 1e-6, hi - 1e-6), bounds=(lo, hi), max_nfev=400, x_scale=0.5)
        if best is None or r.cost < best.cost:
            best = r
    chi2 = float(2 * best.cost)
    n_free = 8 + (3 if design != "D1" else 0) - (0 if design == "D4" else 1)  # unused params do not count
    dof = max(y.size - n_free, 1)
    p = float(stats.chi2.sf(chi2, dof))
    return {"design": design, "sd_rel": sd_rel, "n_obs": int(y.size), "chi2_adsorption_only": chi2, "dof": dof,
            "p_adsorption_only": p, "indistinguishable": bool(p > 0.05), "share_red": base["share_red"]}


# ---------------------------------------------------------------------------
# v0.3: technique-aware Fisher analysis and a formal equivalence test
# ---------------------------------------------------------------------------
TECH_PRIOR_SD = {"XPS": (0.4 / np.log(10), 0.15 / np.sqrt(12)), "XANES": (None, 0.03 / np.sqrt(12))}


def fisher_study_v3(m: sim.Material, meta: dict, h=0.02, designs=DESIGNS_V3):
    """Jacobians for the v0.3 designs w.r.t. (theta, log10 study bias, technique nuisances).

    Derived quantities: reduction share and precipitation share (isotherm series).
    """
    th0 = _theta(m)
    act = active_mask(m)
    base = observe(m, meta)
    cols, g_red, g_pre = [], [], []
    for i in range(len(th0)):
        if not act[i]:
            continue
        th = th0.copy(); th[i] += h
        o = observe(_with(m, th), meta)
        cols.append((i, o))
        g_red.append((o["share_red"] - base["share_red"]) / h)
        g_pre.append((o["share_precip"] - base["share_precip"]) / h)
    out = {"base": base, "active": [PARAMS[i] for i, _ in cols]}
    for d in designs:
        y0, sd = stack(base, d)
        Jc = [(stack(o, d)[0] - y0) / h for _, o in cols] + [y0 * 0 + np.r_[y0[:_n_conc(base, d)] * np.log(10), np.zeros(len(y0) - _n_conc(base, d))]]
        prior = [PRIOR_SD[p] for p in out["active"]] + [0.05 / np.log(10)]
        extra = 0
        if d in TECH_OF:
            t = TECH_OF[d]
            e0, b0 = TECH_NOMINAL[t]
            sd_le, sd_b = TECH_PRIOR_SD[t]
            if sd_le is not None:
                Jc.append((stack(base, d, tech=(e0 * 10 ** h, b0))[0] - y0) / h); prior.append(sd_le); extra += 1
            Jc.append((stack(base, d, tech=(e0, b0 + h * 0.1))[0] - y0) / (h * 0.1)); prior.append(sd_b); extra += 1
        J = np.column_stack(Jc)
        out[d] = (J, sd, np.array(prior), np.array(g_red + [0.0] + [0.0] * extra), np.array(g_pre + [0.0] + [0.0] * extra))
    return out


def _n_conc(base, design):
    n = len(base["cv_eq"])
    if design != "D1":
        n += len(base["cv_kin"])
    if design not in ("D1", "D2"):
        n += len(base["ct_eq"]) + len(base["ct_kin"])
    return n


def posterior_sd_v3(fs: dict, design: str):
    J, sd, prior, g_red, g_pre = fs[design]
    W = 1.0 / sd ** 2
    Spost = np.linalg.inv(J.T @ (W[:, None] * J) + np.diag(1.0 / prior ** 2))
    S0 = np.diag(prior ** 2)
    f = lambda g, S: float(np.sqrt(max(g @ S @ g, 0)))
    return {"design": design, "post_sd_red": f(g_red, Spost), "prior_sd_red": f(g_red, S0),
            "post_sd_precip": f(g_pre, Spost), "prior_sd_precip": f(g_pre, S0)}


def detection_tests_v3(rng, fs: dict, meta: dict, n_rep=20):
    """Model-free diagnostics with technique-specific measurement bias (XPS can false-alarm)."""
    o = fs["base"]
    hits = {"dissolved_cr3": 0, "xps_cr3": 0, "xanes_cr3": 0}
    for _ in range(n_rep):
        b = np.exp(rng.normal(0, 0.05))
        cv = np.maximum(o["cv_eq"] * b * (1 + rng.normal(0, 0.03, o["cv_eq"].size)), 0.005)
        ct = np.maximum(o["ct_eq"] * b * (1 + rng.normal(0, 0.03, o["ct_eq"].size)), 0.01)
        sdd = np.sqrt((0.03 * ct) ** 2 + (0.03 * cv) ** 2 + 0.02 ** 2)
        hits["dissolved_cr3"] += int((ct - cv).sum() / np.sqrt((sdd ** 2).sum()) > stats.norm.ppf(0.99))
        hits["xps_cr3"] += int(ms.measure_oxidation_state(rng, o["xps"], "XPS").max() > 0.20)
        hits["xanes_cr3"] += int(ms.measure_oxidation_state(rng, o["xps"], "XANES").max() > 0.10)
    return {k: v / n_rep for k, v in hits.items()}


def _empirical_predict_v3(x, meta, base, design):
    """Adsorption-only empirical model for the v0.3 designs (see _empirical_predict).

    Extra parameter x[11]: apparent solid Cr(III) fraction produced by the measurement itself
    (beam/vacuum reduction), bounded by the technique's assumed maximum. Under adsorption only, the
    solid inventory equals (C0 - Ce,VI)/D; phase evidence is undefined and is not used.
    """
    core = _empirical_predict(x[:11], meta, base, "D2" if design != "D1" else "D1")
    if design in ("D1", "D2"):
        return core
    b = 10 ** x[7]
    iso, edge, kin_t = design_points(meta)
    n_eq = len(iso) + len(edge)
    cv = core[:n_eq] / b
    cvk = core[n_eq:] / b
    parts = [core, b * cv, b * cvk]
    if design in TECH_OF:
        pH0 = np.array([p for _, p in iso + edge])
        parts += [np.full(len(base["xps"]), x[11]), pH0 + x[8] * (meta["pHpzc"] - pH0)]
        if design == "D5":
            C0s = np.array([c for c, _ in iso + edge])
            idx = [0, len(iso) - 1, len(iso)]
            parts.append((C0s[idx] - cv[idx]) / meta["D"])
    return np.concatenate(parts)


def _stack_equiv(base, design):
    """Observation vector used by the equivalence test (phase evidence excluded; see above)."""
    y, sd = stack(base, design)
    if design == "D5":
        n = len(base["pfrac"])
        y, sd = y[:-n], sd[:-n]
    return y, sd


def equivalence_test_v3(m: sim.Material, meta: dict, design: str, sd_rel=0.03, seed=0, n_starts=3,
                        margins=(0.05, 0.10, 0.20), margin_ox=0.05, alpha=0.05, x_init=None):
    """Goodness of fit and a formal equivalence test of an adsorption-only model.

    Weighted least squares with independent, heteroscedastic errors: SD = sd_rel * y + floor for
    concentrations (floor 0.01 mg/L Cr(VI), 0.02 total), the technique noise for oxidation state,
    0.05 for pH; the shared calibration bias is a free parameter. Observations below the LOD enter
    as LOD/2 with SD = LOD. dof = n_obs - n_free.

    Returns chi2_obs (noisy data), lambda_true (chi2 of the best adsorption-only fit to the
    noise-free data = non-centrality of the misfit, the operating-characteristic truth) and, for
    each margin eps, lambda0(eps) = sum((eps*y_i or margin_ox)/sd_i)^2 and the TOST-type decision
    'equivalent within eps' if chi2_obs <= ncx2.ppf(alpha, dof, lambda0).
    """
    rng = np.random.default_rng(seed)
    base = observe(m, meta)
    y, _ = _stack_equiv(base, design)
    nc = _n_conc(base, design)
    nrest = y.size - nc
    # true technique bias for this study (drawn), measurement noise
    tech = TECH_OF.get(design)
    e, bm = ms.draw_technique_bias(rng, tech) if tech else (1.0, 0.0)
    y_true_obs, _ = stack(base, design, tech=(e, bm)) if tech else (y, None)
    if design == "D5":
        y_true_obs = y_true_obs[:-len(base["pfrac"])]
    floors = np.r_[np.full(nc, 0.01), np.zeros(nrest)]
    sd = np.r_[sd_rel * y[:nc] + floors[:nc], _stack_equiv(base, design)[1][nc:]]
    bias = np.exp(rng.normal(0, 0.05))
    yobs = y_true_obs.copy()
    yobs[:nc] = y_true_obs[:nc] * bias * (1 + rng.normal(0, sd_rel, nc))
    yobs[nc:] = y_true_obs[nc:] + rng.normal(0, 1, nrest) * sd[nc:]
    lod = yobs[:nc] < 0.01
    yobs[:nc][lod] = 0.005
    sd_fit = sd.copy(); sd_fit[:nc][lod] = 0.01
    beam_max = ms.TECHNIQUES[tech]["beam_max"] if tech else 0.0
    lo = np.array([-1.0, 0.2, -4.0, -2.0, 1.0, np.log10(0.1), -4.0, -0.15, 0.0, -2.0, 0.2, 0.0])
    hi = np.array([3.0, 3.0, 2.0, 6.0, 10.0, np.log10(5.0), 0.0, 0.15, 1.0, 3.0, 2.5, max(beam_max, 1e-6)])
    n_free = 8 + (0 if design == "D1" else 3) - (1 if design not in TECH_OF else 0) + (1 if tech else 0)

    def fit_to(target, sdv, inits=()):
        """Bounded WLS from random starts plus warm starts; returns (chi2, x)."""
        f = lambda x: (target - _empirical_predict_v3(x, meta, base, design)) / sdv
        starts = []
        for k in range(n_starts):
            x0 = lo + (hi - lo) * rng.uniform(0.2, 0.8, lo.size)
            x0[0] = np.log10(max(m.qmax, 0.1)); x0[7] = 0.0
            starts.append(x0)
        starts += [np.asarray(x, float) for x in inits if x is not None]
        best = None
        for x0 in starts:
            r = least_squares(f, np.clip(x0, lo + 1e-9, hi - 1e-9), bounds=(lo, hi), max_nfev=600, x_scale=0.5)
            if best is None or r.cost < best.cost:
                best = r
        return float(2 * best.cost), best.x

    chi2, x_noisy = fit_to(yobs, sd_fit, (x_init,))
    lam_true, x_true = fit_to(y_true_obs, sd, (x_noisy, x_init))      # noise-free misfit of the best fit
    chi2_b, x_b = fit_to(yobs, sd_fit, (x_true,))                      # cross-seeded refit of the noisy data
    if chi2_b < chi2:
        chi2, x_noisy = chi2_b, x_b
    dof = max(y.size - n_free, 1)
    out = {"design": design, "sd_rel": sd_rel, "n_obs": int(y.size), "n_free": int(n_free), "dof": int(dof),
           "chi2": chi2, "lambda_true": lam_true, "p_gof": float(stats.chi2.sf(chi2, dof)),
           "share_red": base["share_red"], "share_precip": base["share_precip"],
           "tech_enrich": e, "tech_beam": bm, "_x": x_noisy}
    for eps in margins:
        mvec = np.r_[eps * y[:nc], np.full(nrest, margin_ox)]
        lam0 = float(np.sum((mvec / sd) ** 2))
        crit = float(stats.ncx2.ppf(alpha, dof, lam0))
        out[f"lambda0_{eps}"] = lam0
        out[f"equiv_{eps}"] = bool(chi2 <= crit)
        out[f"truly_within_{eps}"] = bool(lam_true < lam0)
    return out
