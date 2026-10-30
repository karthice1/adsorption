"""Study-aware validation, chemistry-feature ablation and uncertainty calibration."""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from scipy import stats
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

BASE = ["C0", "pH", "D", "t", "T", "I", "SBET", "Vp", "dp", "pHpzc"]
CHEM = BASE + ["f_HCrO4", "f_CrO4", "f_Cr2O7", "z_mean_CrVI", "I_calc", "pH_minus_pzc"]
OFFSET = 0.1  # mg/g, for the log target
LEVELS = np.array([0.5, 0.6, 0.7, 0.8, 0.9, 0.95])


def to_target(q):
    return np.log10(np.clip(q, 0, None) + OFFSET)


def from_target(y):
    return 10 ** y - OFFSET


def make_model(name, seed=0):
    if name == "RF":
        return RandomForestRegressor(n_estimators=400, min_samples_leaf=2, random_state=seed, n_jobs=-1)
    if name == "XGB":
        if XGBRegressor is None:
            raise ImportError("xgboost not installed")
        return XGBRegressor(n_estimators=500, max_depth=4, learning_rate=0.05, subsample=0.8,
                            colsample_bytree=0.8, random_state=seed, n_jobs=4)
    if name == "GPR":
        k = ConstantKernel(1.0) * Matern(length_scale=np.ones(1), nu=2.5) + WhiteKernel(0.01)
        return make_pipeline(StandardScaler(), _ARDGPR(k))
    raise ValueError(name)


class _ARDGPR(GaussianProcessRegressor):
    """GPR whose Matern length scale is expanded to one per feature at fit time (ARD)."""

    def __init__(self, kernel=None):
        super().__init__(kernel=kernel, normalize_y=True, n_restarts_optimizer=1, random_state=0)

    def fit(self, X, y):
        k = ConstantKernel(1.0) * Matern(length_scale=np.ones(X.shape[1]), nu=2.5,
                                         length_scale_bounds=(1e-2, 1e3)) + WhiteKernel(0.01, (1e-5, 1))
        self.kernel = k
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return super().fit(X, y)


def _predict_with_sd(model, X):
    """Point prediction plus a model-internal predictive SD (target scale)."""
    if hasattr(model, "steps"):  # GPR pipeline
        Xs = model[:-1].transform(X)
        mu, sd = model[-1].predict(Xs, return_std=True)
        return mu, sd
    if isinstance(model, RandomForestRegressor):
        per_tree = np.stack([t.predict(np.asarray(X)) for t in model.estimators_])
        return per_tree.mean(0), per_tree.std(0)
    return model.predict(X), None


def splitter(scheme, df, seed=0):
    g = df["study_id"].values
    if scheme == "random":
        return KFold(5, shuffle=True, random_state=seed).split(df)
    if scheme == "grouped5":
        return GroupKFold(5).split(df, groups=g)
    if scheme == "LOSO":
        return LeaveOneGroupOut().split(df, groups=g)
    if scheme == "LOFO":
        return LeaveOneGroupOut().split(df, groups=df["material_family"].values)
    raise ValueError(scheme)


def metrics(q_true, q_pred):
    e = q_pred - q_true
    return {"RMSE": float(np.sqrt(np.mean(e ** 2))), "MAE": float(np.mean(np.abs(e))),
            "MedAE": float(np.median(np.abs(e))),
            "R2": float(1 - np.sum(e ** 2) / np.sum((q_true - q_true.mean()) ** 2))}


def _group_conformal_scores(model_name, X, y, groups, seed=0):
    """Absolute out-of-study residuals (target scale) from inner grouped CV on training studies."""
    scores, study_means = np.empty(len(y)), []
    for tr, te in GroupKFold(min(5, len(np.unique(groups)))).split(X, y, groups):
        m = make_model(model_name, seed).fit(X.iloc[tr], y[tr])
        r = y[te] - m.predict(X.iloc[te])
        scores[te] = np.abs(r)
        for gname in np.unique(groups[te]):
            study_means.append(r[groups[te] == gname].mean())
    return scores, float(np.std(study_means, ddof=1))


def evaluate(df, model_name, features, scheme, target="q_app_VI", uq=None, seed=0):
    """Out-of-fold predictions; ``uq='conformal'`` adds study-aware intervals."""
    X, y = df[features], to_target(df[target].values)
    groups = df["study_id"].values
    pred = np.full(len(y), np.nan)
    lo = np.full((len(LEVELS), len(y)), np.nan)
    hi = lo.copy()
    taus = []
    for tr, te in splitter(scheme, df, seed):
        m = make_model(model_name, seed).fit(X.iloc[tr], y[tr])
        mu, sd = _predict_with_sd(m, X.iloc[te])
        pred[te] = mu
        if uq == "internal" and sd is not None:
            z = stats.norm.ppf(0.5 + LEVELS / 2)[:, None]
            lo[:, te], hi[:, te] = mu - z * sd, mu + z * sd
        elif uq == "conformal":
            s, tau = _group_conformal_scores(model_name, X.iloc[tr], y[tr], groups[tr], seed)
            taus.append(tau)
            n = len(s)
            for i, a in enumerate(LEVELS):
                qa = np.quantile(s, min(1.0, np.ceil((n + 1) * a) / n))
                lo[i, te], hi[i, te] = mu - qa, mu + qa
    q_true = df[target].values
    res = metrics(q_true, np.clip(from_target(pred), 0, None))
    res.update({"model": model_name, "features": "chem" if "f_HCrO4" in features else "base",
                "scheme": scheme, "uq": uq or "none"})
    if uq:
        cov = ((y >= lo) & (y <= hi)).mean(1)
        res.update({"PICP90": float(cov[LEVELS == 0.9][0]),
                    "ACE": float(np.mean(np.abs(cov - LEVELS))),
                    "MPIW90_log": float(np.mean((hi - lo)[LEVELS == 0.9]))})
        if taus:
            res["tau_between_study_log"] = float(np.mean(taus))
    return res, pred


def between_study_variance(df, pred, target="q_app_VI"):
    """Random-intercept decomposition of out-of-study residuals (log scale), method of moments."""
    r = to_target(df[target].values) - pred
    g = pd.Series(r).groupby(df["study_id"].values)
    k, n_bar = g.ngroups, g.size().mean()
    msb = n_bar * g.mean().var(ddof=1)
    msw = g.var(ddof=1).mean()
    tau2 = max((msb - msw) / n_bar, 0.0)
    return {"tau2": tau2, "sigma2": msw, "ICC": tau2 / (tau2 + msw)}
