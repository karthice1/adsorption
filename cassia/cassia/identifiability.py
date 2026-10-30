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
from . import simulate as sim

LOG_PARAMS = ["qmax", "K0", "ka", "R0", "kB", "kC", "L0", "kD", "krel", "q3max", "K3", "k3"]
LIN_PARAMS = ["pH50_red", "pH50_cr3"]
PARAMS = LOG_PARAMS + LIN_PARAMS
# prior SD: log10 units for LOG_PARAMS, pH units for LIN_PARAMS (wide, weakly informative)
PRIOR_SD = {**{p: 0.5 for p in LOG_PARAMS}, "pH50_red": 0.6, "pH50_cr3": 0.45}
DESIGNS = ("D1", "D2", "D3", "D4")
DESIGN_LABEL = {"D1": "Cr(VI) at equilibrium (class C)", "D2": "+ Cr(VI) kinetics (class C)",
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
    cv, ct, ph, xps, red, rem, prec = [], [], [], [], 0.0, 0.0, 0.0
    for k, (C0, pH0) in enumerate(iso + edge):
        tr = sim.run_batch(m, C0, pH0 + pH_shift, D, [teq], T, I)
        cv.append(tr["Ce_VI_true"][0]); ct.append(tr["Ce_VI_true"][0] + tr["Cr3_aq_true"][0])
        ph.append(tr["pH_true"][0])
        if k in (0, len(iso) - 1, len(iso)):
            xps.append(tr["xps_cr3_true"][0])
        if k < len(iso):
            red += tr["reduced_mgL"][0]; rem += C0 - tr["Ce_VI_true"][0]; prec += tr["Cr3_precip_true"][0]
    tk = sim.run_batch(m, 50.0, meta["pH_opt"] + pH_shift, D, kin_t, T, I)
    out = {"cv_eq": np.array(cv), "cv_kin": tk["Ce_VI_true"], "ct_eq": np.array(ct),
           "ct_kin": tk["Ce_VI_true"] + tk["Cr3_aq_true"], "xps": np.array(xps), "pH": np.array(ph)}
    out["share_red"] = red / max(rem, 1e-12)
    out["share_precip"] = prec / max(rem, 1e-12)
    return out


def stack(obs: dict, design: str, log_bias=0.0):
    """Observation vector and SDs for one design (concentrations carry the study bias)."""
    b = 10 ** log_bias
    parts, sds = [b * obs["cv_eq"]], [0.03 * b * obs["cv_eq"] + 0.01]
    if design in ("D2", "D3", "D4"):
        parts.append(b * obs["cv_kin"]); sds.append(0.03 * b * obs["cv_kin"] + 0.01)
    if design in ("D3", "D4"):
        parts += [b * obs["ct_eq"], b * obs["ct_kin"]]
        sds += [0.03 * b * obs["ct_eq"] + 0.02, 0.03 * b * obs["ct_kin"] + 0.02]
    if design == "D4":
        parts += [obs["xps"], obs["pH"]]
        sds += [np.full(len(obs["xps"]), 0.05), np.full(len(obs["pH"]), 0.05)]
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
