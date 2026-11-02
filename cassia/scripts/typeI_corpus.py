"""Corpus-level analysis of the type-I error of the adsorption-only rejection test (Table S14b).

    PYTHONPATH=. python scripts/typeI_corpus.py   -> results/v3/T_typeI_corpus.csv, typeI_corpus.json

Uses only existing outputs: the 360 negative controls of 18 corpora (results/v3/equiv3.csv.gz, arm control_A,
6 corpora x 20; results/v3/typeI_ctrl_part*.csv.gz, 12 corpora x 20). No model is refitted (~1 min).
All control corpora are adsorption-only corpora drawn i.i.d. from the same generator (scenario switches are
set to A after the material draw), so corpora are exchangeable; studies within a corpus share hyper-priors.

Four analyses per design x analytical noise:
  1 study-pooled threshold, leave-one-corpus-out (LOCO): the procedure of scripts/typeI_validation.py,
    reported as EMPIRICAL average cross-corpus type-I error; per-corpus rates, the number of corpora at
    >= 15 %, and a binomial reference (what per-corpus maxima would look like if every control had an
    independent 5 % rejection probability, given 20 controls per corpus).
  2 corpus-level split conformal (cassia.validation.cluster_conformal_threshold): one control per corpus.
    With K = 18 corpora the 95 % threshold is unbounded (needs K >= 19); the smallest attainable level is
    1/(K+1) = 5.26 %, whose threshold is the maximum of the 18 draws. Reported: the threshold (mean over
    B random draws), its LOCO type-I (K = 17, level 1/18) and power on the study arm of stage equiv3.
  3 between-corpus heterogeneity of the null statistic: ICC of log chi2 (one-way ANOVA), F test.
  4 calibration-set sensitivity, to explain the D1 rates on the independent cohort: repeated random splits
    into 6 calibration corpora (120 controls) and 12 validation corpora; study-pooled conformal threshold.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from cassia.validation import cluster_conformal_threshold

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results/v3"
B = 2000          # random one-per-corpus draws (analysis 2) and random calibration splits (analysis 4)
ALPHA = 0.05


def conformal_threshold(x, alpha=ALPHA):
    x = np.sort(np.asarray(x)); k = int(np.ceil((1 - alpha) * (len(x) + 1)))
    return float(x[k - 1]) if k <= len(x) else np.inf


def load():
    Q = pd.read_csv(RES / "equiv3.csv.gz")
    V = pd.concat([pd.read_csv(f) for f in sorted(RES.glob("typeI_ctrl_part*.csv.gz"))], ignore_index=True)
    Qc = Q[Q.arm == "control_A"].assign(corpus=lambda x: x.scenario + "-" + x.rep.astype(str), cohort="calibration")
    Vc = V.assign(corpus="V-" + V.rep.astype(str), cohort="validation")
    cols = ["corpus", "cohort", "sd_rel", "design", "chi2"]
    return pd.concat([Qc[cols], Vc[cols]], ignore_index=True), Q[Q.arm == "true"]


def icc_log(g):
    y = np.log(g.chi2.clip(lower=1e-9)); grp = g.corpus
    m = y.groupby(grp).mean(); n = y.groupby(grp).size(); K, N = len(m), len(y)
    msb = (n * (m - y.mean()) ** 2).sum() / (K - 1); msw = ((y - grp.map(m)) ** 2).sum() / (N - K)
    n0 = (N - (n ** 2).sum() / N) / (K - 1)
    return float((msb - msw) / (msb + (n0 - 1) * msw)), float(stats.f.sf(msb / msw, K - 1, N - K))


def main():
    A, T = load()
    rng = np.random.default_rng(20261102)
    corpora = np.array(sorted(A.corpus.unique()))
    cal_corpora = sorted(A[A.cohort == "calibration"].corpus.unique())
    K = len(corpora)
    rows, cells_split = [], {}
    for (sd, d), g in A.groupby(["sd_rel", "design"]):
        rec = {"sd_rel": sd, "design": d, "n_corpora": K, "n_controls": len(g)}
        by = {c: h.chi2.values for c, h in g.groupby("corpus")}
        # ---- 1 study-pooled LOCO (empirical)
        kk, nn = [], []
        for c in corpora:
            thr = conformal_threshold(g[g.corpus != c].chi2.values)
            kk.append(int((by[c] > thr).sum())); nn.append(len(by[c]))
        kk, nn = np.array(kk), np.array(nn); r = kk / nn
        boot = [kk[i].sum() / nn[i].sum() for i in (rng.integers(0, K, K) for _ in range(4000))]
        p_bar = kk.sum() / nn.sum(); n_c = int(np.median(nn))
        rec |= {"pooled_loco_typeI": p_bar, "pooled_loco_lo": float(np.quantile(boot, 0.025)),
                "pooled_loco_hi": float(np.quantile(boot, 0.975)), "pooled_max_corpus": float(r.max()),
                "pooled_min_corpus": float(r.min()), "pooled_n_corpora_ge_0.15": int((r >= 0.15 - 1e-12).sum()),
                "pooled_n_corpora_le_0.05": int((r <= 0.05 + 1e-12).sum()),
                # binomial reference: P(max over K corpora >= observed max) and expected number of corpora at >= 15 %
                # if each of the n_c controls in a corpus were rejected independently with probability 0.05
                "binom_P_max_ge_observed": float(1 - stats.binom.cdf(round(r.max() * n_c) - 1, n_c, ALPHA) ** K),
                "binom_expected_n_ge_0.15": float(K * stats.binom.sf(np.ceil(0.15 * n_c) - 1, n_c, ALPHA)),
                # overdispersion of per-corpus rejection counts relative to binomial(n_c, p_bar)
                "dispersion_chi2_p": float(stats.chi2.sf(((kk - nn * p_bar) ** 2 / (nn * p_bar * (1 - p_bar))).sum(), K - 1))}
        # ---- 3 heterogeneity of the null statistic
        rec["icc_log_chi2"], rec["icc_F_p"] = icc_log(g)
        # ---- 2 corpus-level split conformal
        sc, cl = g.chi2.values, g.corpus.values
        rec["corpus_thr_alpha0.05"] = cluster_conformal_threshold(sc, cl, ALPHA, rng)   # inf with K = 18
        rec["corpus_alpha_min"] = 1 / (K + 1)
        thr_c = np.array([cluster_conformal_threshold(sc, cl, 1 / (K + 1), rng) for _ in range(B)])
        rec["corpus_thr_mean"] = float(thr_c.mean()); rec["corpus_thr_p05"] = float(np.quantile(thr_c, 0.05))
        rec["corpus_thr_p95"] = float(np.quantile(thr_c, 0.95))
        loco_c = []
        for c in corpora:
            m = cl != c
            t = np.array([cluster_conformal_threshold(sc[m], cl[m], 1 / K, rng) for _ in range(B // 10)])
            loco_c.append(float((by[c][None, :] > t[:, None]).mean()))
        loco_c = np.array(loco_c)
        rec |= {"corpus_loco_level": 1 / K, "corpus_loco_typeI": float(loco_c.mean()), "corpus_loco_max_corpus": float(loco_c.max())}
        tru = T[(T.sd_rel == sd) & (T.design == d)]
        thr_pooled = conformal_threshold(sc)
        for lo_, hi_, lab in ((0.05, 0.3, "5-30%"), (0.3, 1.01, ">30%")):
            x = tru[(tru.share_red > lo_) & (tru.share_red <= hi_)].chi2.values
            rec[f"n_true|{lab}"] = len(x)
            rec[f"power_pooled|{lab}"] = float((x > thr_pooled).mean())
            rec[f"power_corpus|{lab}"] = float((x[None, :] > thr_c[:, None]).mean())
        # ---- 4 calibration-set sensitivity (6 calibration corpora, 12 validation corpora)
        ti = []
        for _ in range(B):
            cs = rng.choice(corpora, 6, replace=False); m = np.isin(cl, cs)
            ti.append(float((sc[~m] > conformal_threshold(sc[m])).mean()))
        ti = np.array(ti)
        m = np.isin(cl, cal_corpora); obs = float((sc[~m] > conformal_threshold(sc[m])).mean())
        cells_split[(sd, d)] = ti
        # same split sizes but studies (not corpora) allocated at random
        ts = []
        for _ in range(B):
            idx = rng.permutation(len(sc)); ts.append(float((sc[idx[120:]] > conformal_threshold(sc[idx[:120]])).mean()))
        rec |= {"split6_observed_typeI": obs, "split6_median": float(np.median(ti)), "split6_p025": float(np.quantile(ti, 0.025)),
                "split6_p975": float(np.quantile(ti, 0.975)), "split6_P_ge_observed": float((ti >= obs - 1e-12).mean()),
                "split6_P_ge_0.10": float((ti >= 0.10 - 1e-12).mean()),
                "splitstudy_p975": float(np.quantile(ts, 0.975)), "splitstudy_P_ge_observed": float((np.array(ts) >= obs - 1e-12).mean())}
        rows.append(rec)
    R = pd.DataFrame(rows)
    R.to_csv(RES / "T_typeI_corpus.csv", index=False)
    # probability that a random 6-corpus calibration gives >= 14.2 % (the observed D1 rate) in at least one of
    # the 12 design x noise cells (same random splits are not shared between cells, so this is approximate)
    M = np.column_stack([cells_split[k] for k in sorted(cells_split)])
    summ = {"K_corpora": K, "alpha": ALPHA, "corpus_level_finite_needs_K_ge": int(np.ceil((1 - ALPHA) / ALPHA)),
            "alpha_min_K18": 1 / (K + 1), "P_any_cell_ge_0.142_approx": float((M >= 0.142 - 1e-12).any(axis=1).mean()),
            "calibration_corpora": cal_corpora, "B": B}
    (RES / "typeI_corpus.json").write_text(json.dumps(summ, indent=1))
    pd.set_option("display.width", 250)
    print(R.round(3).T.to_string()); print(summ)


if __name__ == "__main__":
    main()
