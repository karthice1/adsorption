## Appendix A. Coverage guarantee of the study-level split-conformal interval (M2s)

**Setting.** Let the studies used in one leave-one-study-out fold be divided at random into a *proper-training* set and $J$ *calibration* studies, and let $Z_{J+1}$ be the held-out test study. Each study is a set of observations,

$$Z_j=\{(x_{ji},y_{ji}):\ i=1,\dots,n_j\}.$$

The study size $n_j$ is part of the study: it is random, can vary between studies, and the observations within a study may be arbitrarily dependent. Let $\hat\mu$ be a model fitted to the proper-training studies only. For a fixed $\beta\in(0,1]$ define:

- the residuals $R_{ji}=|y_{ji}-\hat\mu(x_{ji})|$;
- the study score $s_j$, the $k_j$-th smallest of the residuals $R_{ji}$ of study $j$, with $k_j=\lceil\beta n_j\rceil$.

For a miscoverage level $\alpha$, let $k=\lceil(1-\alpha)(J+1)\rceil$. Set $\hat q=s_{(k)}$, the $k$-th smallest of the $J$ calibration scores, if $k\le J$, and $\hat q=\infty$ otherwise. The M2s interval at any $x$ is $\hat\mu(x)\pm\hat q$.

**Assumption.** Conditional on the proper-training studies, $Z_1,\dots,Z_J,Z_{J+1}$ are exchangeable. Independent draws of studies from a common population of studies, as in the simulated corpora, satisfy this assumption.

**Theorem.** Under the assumption,

$$P\Big(\#\{i:\ |y_{J+1,i}-\hat\mu(x_{J+1,i})|\le\hat q\}\ \ge\ \lceil\beta n_{J+1}\rceil\Big)\ \ge\ 1-\alpha .$$

In words: with probability at least $1-\alpha$, at least $\lceil\beta n_{J+1}\rceil$ of the $n_{J+1}$ observations of a new study lie inside its intervals, so its within-study coverage is at least $\lceil\beta n\rceil/n\ge\beta$. If the scores are almost surely distinct, the probability is also at most $1-\alpha+1/(J+1)$.

**Proof.**
1. $\hat\mu$ does not depend on $Z_1,\dots,Z_{J+1}$, and every score is the same measurable function of its study. Hence $s_1,\dots,s_{J+1}$ are exchangeable.
2. By the standard quantile lemma of split conformal prediction (Vovk et al., 2005; Lei et al., 2018), the rank of $s_{J+1}$ among the $J+1$ scores, with ties broken at random, is uniform on $\{1,\dots,J+1\}$. Therefore $P(s_{J+1}\le s_{(k)})\ge k/(J+1)\ge 1-\alpha$.
3. On the event $s_{J+1}\le\hat q$, by the definition of $s_{J+1}$ as an order statistic, at least $\lceil\beta n_{J+1}\rceil$ residuals of the test study are $\le s_{J+1}\le\hat q$.
4. The upper bound follows from the uniform rank when the scores are distinct. ∎

**Remarks.**
1. *Marginal, not conditional.* The guarantee is marginal over the population of studies, including its mixture of material families and study sizes. It does not hold conditionally within a family. The family-level coverage reported in §3.5 and Table S5 is therefore descriptive.
2. *Effect of small studies.* $\lceil\beta n\rceil/n$ is conservative for small studies. For $\beta=0.9$, a study of 28 observations needs 26 inside (92.9%), and a study of 8 observations needs all 8.
3. *Minimum number of calibration studies.* The interval is finite only if $J\ge(1-\alpha)/\alpha$: at least 9 calibration studies for $\alpha=0.10$ and 19 for $\alpha=0.05$. With half of the 29 training studies used for calibration (about 14), M2s is bounded at 90% and unbounded at 95%. For total-Cr endpoints only about 8 calibration studies are available, so M2s is unbounded even at 90%. This finite-sample limit is reported as a result, not hidden.
4. *Relation to Dunn et al. (2023).* Their methods target coverage of a *single random observation* from a new group, using pooling, subsampling or repeated subsampling. M1 is our analogue of their weighted pooling. M2s targets a different, group-level $(\beta,1-\alpha)$ tolerance statement: a fraction $\beta$ of a new study's observations, with confidence $1-\alpha$.
5. *M2 (CV+) has no proven guarantee.* M2 uses the same study scores, but with K-fold models instead of a single split, as in the CV+ construction (Barber et al., 2021). The CV+ guarantee is proven for scalar residual scores of exchangeable data points. We have not shown that it extends to within-study quantile scores, so no $1-2\alpha$ guarantee is claimed for M2. M2 is reported as an empirical method only.

## Appendix B. Goodness-of-fit and equivalence tests of the adsorption-only model

**Observations and weights.**
- *Vector.* The observation vector of a design (§2.6) contains residual Cr(VI) and, as the design requires, total aqueous Cr, the solid oxidation-state fraction, final pH and solid Cr inventory. Phase evidence is excluded, because it is undefined under adsorption only.
- *Error model.* Each observation $y_i$ carries an independent error with SD $\sigma_i$:
  - $\sigma_{\rm rel}y_i+0.01$ mg L⁻¹ for Cr(VI), with $\sigma_{\rm rel}$ = 0.03 or 0.10;
  - $\sigma_{\rm rel}y_i+0.02$ mg L⁻¹ for total Cr;
  - the technique noise (0.05 XPS, 0.03 XANES) for oxidation state;
  - 0.05 for pH;
  - 5% + 0.01 mg g⁻¹ for inventory.
- *Shared calibration error.* A calibration bias shared by all concentrations of a study is a free parameter, so residual correlation from calibration is absorbed rather than ignored.
- *Below the detection limit.* Observations below the limit of detection are set to LOD/2 with $\sigma_i$ = LOD.

**The adsorption-only model and its fit.**
- *Model.* The empirical adsorption-only model (§2.6) has:
  - a pH-dependent Sips equilibrium (5 parameters plus $q_{\max}$);
  - a free kinetic curve (3 parameters; D2 and later);
  - the study bias;
  - a free pH-relaxation fraction (A1/A2 designs);
  - a free apparent solid Cr(III) fraction produced by the measurement itself, bounded by the technique's assumed maximum beam reduction (A1/A2).
- *Fit.* The model is fitted by bounded weighted least squares with two to three random starts.
- *Statistic.* The misfit statistic is $\chi^2=\sum_i r_i^2/\sigma_i^2$ at the optimum, with $\nu=n_{\rm obs}-n_{\rm free}$ degrees of freedom.

**Calibrated rejection.**
- *Why the nominal test fails.* The empirical model is not exactly the data-generating adsorption model, so $\chi^2$ under the null is inflated relative to $\chi^2_\nu$. The nominal test therefore over-rejects (6.7–15% at nominal 5%; §3.3).
- *Calibration.* The rejection threshold is instead the 95th percentile of $\chi^2$ in the paired negative controls (same materials, reduction switched off), separately for each design and noise level.
- *Cross-validated type-I error.* The type-I error is estimated by split-half cross-validation of the controls: the threshold is estimated on one half and the error is evaluated on the other.

**Equivalence within a margin.**
- *Hypotheses.* Non-rejection is not evidence of equivalence. We therefore test $H_0{:}\ \lambda\ge\lambda_0$ against $H_1{:}\ \lambda<\lambda_0$, where $\lambda$ is the non-centrality of the misfit, i.e. the $\chi^2$ of the best adsorption-only approximation to the noise-free data.
- *Margin.* The margin is $\lambda_0(\varepsilon)=\sum_i(m_i/\sigma_i)^2$, with $m_i=\varepsilon y_i$ for concentrations and inventory and $m_i$ = 0.05 for oxidation-state fractions. In words, the margin is a systematic discrepancy of relative size $\varepsilon$ in every concentration, with $\varepsilon\in\{0.05,0.10,0.20\}$.
- *Decision.* 'Equivalent within $\varepsilon$' is declared when $\chi^2\le F^{-1}_{\nu,\lambda_0}(\alpha)$, the $\alpha$-quantile of the noncentral $\chi^2$ with $\nu$ degrees of freedom and non-centrality $\lambda_0$, at $\alpha$ = 0.05. This is the goodness-of-fit form of an equivalence test (Wellek, 2010). As in two one-sided tests (Schuirmann, 1987), equivalence must be demonstrated rather than presumed.

**Operating characteristics.**
- *Truth.* Because the simulator gives the noise-free observations, the true non-centrality $\lambda_{\rm true}$ of each study is known.
- *Reported.* We report the false-equivalence rate, $P(\text{declared equivalent}\mid\lambda_{\rm true}\ge\lambda_0)$, and the equivalence power, $P(\text{declared equivalent}\mid\lambda_{\rm true}<\lambda_0)$.
- *Three outcomes per study.* Each study is reported as one of: rejected by the calibrated test; equivalent within $\varepsilon$; or inconclusive, i.e. neither. Only 'rejected' and 'equivalent' are conclusions. "Not rejected" is never interpreted as "structurally equivalent".
