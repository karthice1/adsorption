"""Study-aware validation, ablation, cluster-level conformal prediction (CASSIA v0.2).

Conventions
-----------
* Target: y = log10(max(q, 0) + 0.1 mg/g). Negative apparent uptake (an analytical artefact
  when C_e exceeds C0 within noise) is censored at q = 0 and counted; values below the LOD
  enter as LOD/2 upstream. Back-transformation uses Duan's smearing factor estimated on the
  training residuals (``retransform='smear'``) or the naive inverse.
* Missing ionic strength: imputed at 0.01 M with an ``I_imputed`` indicator feature; the
  source data are never modified.
* Hyperparameters are fixed a priori (``HYPER``) unless ``tune=True`` (nested grouped CV).
* Seeds: every model uses ``random_state=seed`` (default 0).

Conformal methods (all calibrated on out-of-study residuals; unit = study)
-------------------------------------------------------------------------
M0  pooled-residual: half-width = ceil((n+1)(1-a))/n quantile of pooled absolute residuals
    from inner grouped K-fold CV, centred on the full-model prediction. Heuristic: residuals of
    one study are dependent, so no formal guarantee for a new study.
M1  study-level CV+ with study-weighted observation scores: each calibration study carries
    total weight 1 (weight 1/n_j per observation); interval from the CV+ construction of
    Barber et al. (2021) at the level of studies. Targets marginal coverage of a random
    observation of a new study.
M2  study-level CV+ with study scores s_j = beta-quantile of |r_ij| within study j. Targets
    'at least a fraction beta of a new study's observations covered' with probability
    >= 1 - 2a (CV+ bound) under exchangeability of studies.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.model_selection import GroupKFold, KFold, LeaveOneGroupOut
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

try:
    from xgboost import XGBRegressor
except ImportError:  # optional
    XGBRegressor = None

BASE = ["C0", "pH", "D", "t", "T", "log_I", "I_imputed", "SBET", "Vp", "dp", "pHpzc"]
CHEM = BASE + ["f_HCrO4", "f_CrO4", "f_Cr2O7", "z_mean_CrVI", "pH_minus_pzc"]
OFFSET = 0.1
LEVELS = np.array([0.5, 0.6, 0.7, 0.8, 0.9, 0.95])
HYPER = {
    "RF": {"n_estimators": 300, "min_samples_leaf": 2, "max_features": 1.0},
    "XGB": {"n_estimators": 500, "max_depth": 4, "learning_rate": 0.05, "subsample": 0.8, "colsample_bytree": 0.8},
    "GPR": {"kernel": "Constant * Matern(nu=2.5, ARD) + White", "n_restarts_optimizer": 0, "normalize_y": True},
}
TUNE_GRID = {
    "RF": [{"min_samples_leaf": a, "max_features": b} for a in (1, 3, 5) for b in (0.33, 0.66, 1.0)],
    "XGB": [{"max_depth": a, "learning_rate": b} for a in (3, 5) for b in (0.03, 0.1)],
}


def to_target(q):
    return np.log10(np.clip(q, 0, None) + OFFSET)


def from_target(y, smear=1.0):
    return np.clip(10 ** y * smear - OFFSET, 0, None)


def duan_smear(resid_log10):
    """Duan (1983) smearing factor for a log10-scale model."""
    return float(np.mean(10 ** np.asarray(resid_log10)))


class _ARDGPR(GaussianProcessRegressor):
    def __init__(self, kernel=None):
        super().__init__(kernel=kernel, normalize_y=True, n_restarts_optimizer=0, random_state=0)

    def fit(self, X, y):
        self.kernel = ConstantKernel(1.0) * Matern(length_scale=np.ones(X.shape[1]), nu=2.5,
                                                   length_scale_bounds=(1e-2, 1e3)) + WhiteKernel(0.01, (1e-5, 1))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return super().fit(X, y)


def make_model(name, seed=0, **over):
    if name == "RF":
        return RandomForestRegressor(**{**HYPER["RF"], **over}, random_state=seed, n_jobs=1)
    if name == "XGB":
        if XGBRegressor is None:
            raise ImportError("xgboost not installed")
        return XGBRegressor(**{**HYPER["XGB"], **over}, random_state=seed, n_jobs=1, verbosity=0)
    if name == "GPR":
        return make_pipeline(StandardScaler(), _ARDGPR())
    raise ValueError(name)


def prepare(df: pd.DataFrame, target="q_app_VI"):
    """Rows usable for a target, with censoring bookkeeping."""
    d = df[df[target].notna()].copy()
    d["target_censored"] = d[target] < 0
    return d


def splitter(scheme, df, seed=0):
    g = df["study_id"].values
    if scheme == "random":
        return list(KFold(5, shuffle=True, random_state=seed).split(df))
    if scheme == "grouped5":
        return list(GroupKFold(5).split(df, groups=g))
    if scheme == "LOSO":
        return list(LeaveOneGroupOut().split(df, groups=g))
    if scheme == "LOFO":
        return list(LeaveOneGroupOut().split(df, groups=df["material_family"].values))
    raise ValueError(scheme)


def metrics(q_true, q_pred):
    e = q_pred - q_true
    ss = np.sum((q_true - q_true.mean()) ** 2)
    return {"RMSE": float(np.sqrt(np.mean(e ** 2))), "MAE": float(np.mean(np.abs(e))),
            "MedAE": float(np.median(np.abs(e))), "bias": float(np.mean(e)),
            "R2": float(1 - np.sum(e ** 2) / ss) if ss > 0 else np.nan}


def cluster_bootstrap(df, q_true, q_pred, B=1000, seed=0):
    """Study-level bootstrap percentile CIs for pooled R2 and RMSE."""
    rng = np.random.default_rng(seed)
    sid = df["study_id"].values
    groups = {s: np.where(sid == s)[0] for s in np.unique(sid)}
    keys = list(groups)
    r2, rmse = [], []
    for _ in range(B):
        idx = np.concatenate([groups[k] for k in rng.choice(keys, len(keys), replace=True)])
        m = metrics(q_true[idx], q_pred[idx])
        r2.append(m["R2"]); rmse.append(m["RMSE"])
    return {"R2_lo": float(np.nanquantile(r2, 0.025)), "R2_hi": float(np.nanquantile(r2, 0.975)),
            "RMSE_lo": float(np.quantile(rmse, 0.025)), "RMSE_hi": float(np.quantile(rmse, 0.975))}


def per_study(df, q_true, q_pred):
    out = []
    for s in np.unique(df["study_id"]):
        i = (df["study_id"] == s).values
        m = metrics(q_true[i], q_pred[i])
        out.append({"study_id": s, "family": df["material_family"].values[i][0], **m})
    return pd.DataFrame(out)


def _tuned_params(name, X, y, groups, seed):
    best, best_s = {}, np.inf
    for p in TUNE_GRID[name]:
        errs = []
        for tr, te in GroupKFold(4).split(X, y, groups):
            mdl = make_model(name, seed, **p).fit(X.iloc[tr], y[tr])
            errs.append(np.mean((mdl.predict(X.iloc[te]) - y[te]) ** 2))
        if np.mean(errs) < best_s:
            best, best_s = p, np.mean(errs)
    return best


def evaluate(df, model_name, features, scheme, target="q_app_VI", seed=0, tune=False, retransform="smear"):
    """Out-of-fold predictions and metrics (mg/g) with study-level bootstrap CIs."""
    d = prepare(df, target)
    X, y = d[features], to_target(d[target].values)
    groups = d["study_id"].values
    pred_log = np.full(len(y), np.nan)
    smear = np.ones(len(y))
    chosen = []
    for tr, te in splitter(scheme, d, seed):
        over = _tuned_params(model_name, X.iloc[tr], y[tr], groups[tr], seed) if tune else {}
        chosen.append(over)
        mdl = make_model(model_name, seed, **over).fit(X.iloc[tr], y[tr])
        pred_log[te] = mdl.predict(X.iloc[te])
        if retransform == "smear":
            smear[te] = duan_smear(y[tr] - mdl.predict(X.iloc[tr]))
    q_true = np.clip(d[target].values, 0, None)
    q_pred = from_target(pred_log, smear)
    res = metrics(q_true, q_pred) | cluster_bootstrap(d, q_true, q_pred)
    res.update({"model": model_name, "features": "chem" if "f_HCrO4" in features else "base", "scheme": scheme,
                "target": target, "n_obs": len(d), "n_studies": int(d["study_id"].nunique()),
                "n_censored": int(d["target_censored"].sum()), "tuned": tune})
    return res, pred_log, d, q_pred


# ---------------------------------------------------------------------------
# Cluster-level conformal prediction
# ---------------------------------------------------------------------------
def _wquantile(v, w, q):
    o = np.argsort(v)
    v, w = v[o], w[o]
    c = np.cumsum(w) / w.sum()
    return v[min(np.searchsorted(c, q), len(v) - 1)]


class StudyCVPlus:
    """Study-level CV+ for one endpoint (fit on training studies, predict intervals at new X)."""

    def __init__(self, model_name="RF", features=None, K=10, seed=0):
        self.model_name, self.features, self.K, self.seed = model_name, features or CHEM, K, seed

    def fit(self, d: pd.DataFrame, target: str):
        d = prepare(d, target)
        X, y = d[self.features], to_target(d[target].values)
        g = d["study_id"].values
        K = min(self.K, len(np.unique(g)))
        self.fold_models, self.fold_of_study, self.res = [], {}, []
        for k, (tr, te) in enumerate(GroupKFold(K).split(X, y, g)):
            mdl = make_model(self.model_name, self.seed).fit(X.iloc[tr], y[tr])
            self.fold_models.append(mdl)
            r = y[te] - mdl.predict(X.iloc[te])
            for s in np.unique(g[te]):
                i = g[te] == s
                self.res.append({"study": s, "fold": k, "r": r[i]})
        self.full = make_model(self.model_name, self.seed).fit(X, y)
        self.n_studies = len(self.res)
        return self

    def _fold_preds(self, X):
        return np.stack([m.predict(X[self.features]) for m in self.fold_models])  # (K, n)

    def interval(self, X, alpha=0.1, method="M1", beta=0.9):
        """Return (lo, hi, centre) on the log10 target scale."""
        centre = self.full.predict(X[self.features])
        if method == "M0":
            r = np.abs(np.concatenate([e["r"] for e in self.res]))
            n = len(r)
            qa = np.quantile(r, min(1.0, np.ceil((n + 1) * (1 - alpha)) / n))
            return centre - qa, centre + qa, centre
        P = self._fold_preds(X)
        if method == "M1":
            vals_lo, vals_hi, w = [], [], []
            for e in self.res:
                mu = P[e["fold"]]
                a = np.abs(e["r"])
                vals_lo.append(mu[None, :] - a[:, None]); vals_hi.append(mu[None, :] + a[:, None])
                w.append(np.full(len(a), 1.0 / len(a)))
            L, H, W = np.vstack(vals_lo), np.vstack(vals_hi), np.concatenate(w)
            W = np.concatenate([W, [1.0]])  # +1 study mass for the test study (conservative)
            lo = np.array([_wquantile(np.append(L[:, j], -np.inf), W, alpha) for j in range(L.shape[1])])
            hi = np.array([_wquantile(np.append(H[:, j], np.inf), W, 1 - alpha) for j in range(H.shape[1])])
            return lo, hi, centre
        if method == "M2":
            s = np.array([np.quantile(np.abs(e["r"]), beta) for e in self.res])
            mu = np.stack([P[e["fold"]] for e in self.res])            # (J, n)
            L, H = mu - s[:, None], mu + s[:, None]
            J = len(s)
            k_lo = int(np.floor(alpha * (J + 1)))
            k_hi = int(np.ceil((1 - alpha) * (J + 1)))
            Ls, Hs = np.sort(L, 0), np.sort(H, 0)
            lo = Ls[k_lo - 1] if k_lo >= 1 else np.full(L.shape[1], -np.inf)
            hi = Hs[k_hi - 1] if k_hi <= J else np.full(H.shape[1], np.inf)
            return lo, hi, centre
        raise ValueError(method)


def conformal_loso(df, target="q_app_VI", model_name="RF", features=None, K=10, alphas=(0.1,), beta=0.9, seed=0):
    """Outer leave-one-study-out evaluation of M0/M1/M2; per-study coverage table."""
    features = features or CHEM
    d = prepare(df, target)
    y = to_target(d[target].values)
    g = d["study_id"].values
    rows = []
    for test in np.unique(g):
        tr, te = g != test, g == test
        cp = StudyCVPlus(model_name, features, K, seed).fit(d[tr], target)
        for a in alphas:
            for meth in ("M0", "M1", "M2"):
                lo, hi, _ = cp.interval(d[te], a, meth, beta)
                cov = (y[te] >= lo) & (y[te] <= hi)
                rows.append({"study_id": test, "family": d["material_family"].values[te][0], "alpha": a,
                             "method": meth, "coverage": float(cov.mean()), "n": int(te.sum()),
                             "width_log10": float(np.mean(np.minimum(hi - lo, 10))),
                             "n_cal_studies": cp.n_studies})
    return pd.DataFrame(rows)


def between_study_variance(df, pred_log, target="q_app_VI"):
    d = prepare(df, target)
    r = to_target(d[target].values) - pred_log
    gr = pd.Series(r).groupby(d["study_id"].values)
    n_bar = gr.size().mean()
    msb = n_bar * gr.mean().var(ddof=1)
    msw = gr.var(ddof=1).mean()
    tau2 = max((msb - msw) / n_bar, 0.0)
    return {"tau2": tau2, "sigma2": msw, "ICC": tau2 / (tau2 + msw)}
