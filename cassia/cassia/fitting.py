"""Nonlinear isotherm/kinetic fitting and model selection for aqueous Cr data.

Model equations are taken from the ``adsorption`` package (``adsorption.models``)
so that CASSIA and the upstream package share one definition. Two upstream
behaviours are deliberately *not* reused for aqueous data: the gas-phase
concentration transform in ``Isotherms.set_inlet`` (P*Mr/RT), and the Temkin
form with T fixed at 298.15 K. CASSIA fits Ce in mg/L directly and passes T.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit

from adsorption.models import Isotherms, Kinetics

_ISO = Isotherms()
_KIN = Kinetics()
R = 8.314


def temkin_T(T):
    """Temkin isotherm at the experimental temperature (upstream fixes T = 298.15 K)."""
    return lambda Ce, a, b: (R * T / b) * np.log(a * np.asarray(Ce))


ISOTHERMS = {
    # name: (function or factory(T), p0 builder, lower, upper)
    "Langmuir": (_ISO.langmuir, lambda x, y: [1.0 / np.median(x), 1.2 * y.max()], [1e-8, 1e-6], [1e4, 1e5]),
    "Freundlich": (_ISO.freundlich, lambda x, y: [y.max() / x.max() ** 0.5, 2.0], [1e-8, 0.1], [1e5, 50]),
    "Sips": (_ISO.sips, lambda x, y: [1.0 / np.median(x), 1.0, 1.2 * y.max()], [1e-8, 0.1, 1e-6], [1e4, 20, 1e5]),
    "Temkin": (temkin_T, lambda x, y: [10.0 / np.median(x), 100.0], [1e-8, 1e-3], [1e6, 1e7]),
}
KINETICS = {
    "PFO": (_KIN.pfo, lambda x, y: [0.05, y.max()], [1e-8, 1e-6], [10, 1e5]),
    "PSO": (_KIN.pso, lambda x, y: [1e-3, 1.1 * y.max()], [1e-10, 1e-6], [10, 1e5]),
}


@dataclass
class Fit:
    model: str
    params: np.ndarray
    cov: np.ndarray
    n: int
    k: int
    rss: float
    y: np.ndarray
    yhat: np.ndarray
    ci: dict = field(default_factory=dict)

    @property
    def r2(self):
        return 1 - self.rss / np.sum((self.y - self.y.mean()) ** 2)

    @property
    def aicc(self):
        n, k = self.n, self.k + 1  # +1 for the error variance
        aic = n * np.log(self.rss / n) + 2 * k
        return aic + 2 * k * (k + 1) / max(n - k - 1, 1)

    @property
    def bic(self):
        return self.n * np.log(self.rss / self.n) + (self.k + 1) * np.log(self.n)

    def diagnostics(self):
        res = self.y - self.yhat
        sw = stats.shapiro(res).pvalue if self.n >= 3 else np.nan
        signs = np.sign(res[res != 0])
        runs = 1 + np.sum(signs[1:] != signs[:-1]) if len(signs) > 1 else np.nan
        npos, nneg = np.sum(signs > 0), np.sum(signs < 0)
        if npos and nneg:
            mu = 2 * npos * nneg / (npos + nneg) + 1
            var = (mu - 1) * (mu - 2) / (npos + nneg - 1)
            runs_p = 2 * stats.norm.sf(abs((runs - mu) / np.sqrt(var))) if var > 0 else np.nan
        else:
            runs_p = np.nan
        rel_se = np.sqrt(np.clip(np.diag(self.cov), 0, None)) / np.abs(self.params)
        return {"shapiro_p": sw, "runs_p": runs_p, "max_rel_SE": float(np.nanmax(rel_se))}


def _fit(name, func, p0, lo, hi, x, y, sigma=None):
    p0 = np.clip(p0, np.asarray(lo) * 1.01, np.asarray(hi) * 0.99)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p, cov = curve_fit(func, x, y, p0=p0, bounds=(lo, hi), sigma=sigma, maxfev=20000)
    yhat = func(x, *p)
    return Fit(name, p, cov, len(y), len(p), float(np.sum((y - yhat) ** 2)), y, yhat)


def fit_isotherms(Ce, qe, T=298.15, models=None, sigma=None):
    Ce, qe = np.asarray(Ce, float), np.asarray(qe, float)
    out = {}
    for name in models or ISOTHERMS:
        f, p0, lo, hi = ISOTHERMS[name]
        func = f(T) if name == "Temkin" else f
        try:
            out[name] = _fit(name, func, p0(Ce, qe), lo, hi, Ce, qe, sigma)
        except (RuntimeError, ValueError):
            continue
    return out


def fit_kinetics(t, qt, models=None):
    t, qt = np.asarray(t, float), np.asarray(qt, float)
    out = {}
    for name in models or KINETICS:
        f, p0, lo, hi = KINETICS[name]
        try:
            out[name] = _fit(name, f, p0(t, qt), lo, hi, t, qt)
        except (RuntimeError, ValueError):
            continue
    return out


def selection_table(fits: dict) -> pd.DataFrame:
    if not fits:
        return pd.DataFrame(columns=["model", "k", "R2", "RMSE", "AICc", "BIC", "dAICc", "Akaike_w"])
    rows = []
    for name, f in fits.items():
        rows.append({"model": name, "k": f.k, "R2": f.r2, "RMSE": np.sqrt(f.rss / f.n),
                     "AICc": f.aicc, "BIC": f.bic, **f.diagnostics()})
    df = pd.DataFrame(rows).sort_values("AICc")
    df["dAICc"] = df["AICc"] - df["AICc"].min()
    w = np.exp(-0.5 * df["dAICc"])
    df["Akaike_w"] = w / w.sum()
    return df.reset_index(drop=True)


def bootstrap_ci(fit: Fit, x, kind="isotherm", T=298.15, B=500, seed=0, level=0.95):
    """Residual bootstrap percentile CIs for the parameters of one fitted model."""
    rng = np.random.default_rng(seed)
    table = ISOTHERMS if kind == "isotherm" else KINETICS
    f, _, lo, hi = table[fit.model]
    func = f(T) if fit.model == "Temkin" else f
    res = fit.y - fit.yhat
    draws = []
    for _ in range(B):
        yb = fit.yhat + rng.choice(res, size=len(res), replace=True)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                p, _ = curve_fit(func, x, yb, p0=fit.params, bounds=(lo, hi), maxfev=5000)
            draws.append(p)
        except (RuntimeError, ValueError):
            pass
    draws = np.array(draws)
    a = (1 - level) / 2
    return np.quantile(draws, [a, 1 - a], axis=0), draws
