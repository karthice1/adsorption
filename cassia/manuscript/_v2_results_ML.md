### 3.4 Transferability of ML models across studies (synthetic)

Fig. 4 and Table 5 summarise 20 independent scenario-E corpora (30 studies, six families, 840 observations each). The Gaussian process was run on 3 corpora. The target is apparent Cr(VI) uptake. About 40 observations per corpus (4.8%) had negative apparent uptake and were censored at zero.

**Random splits overstate transferability.**

| Validation scheme | RF (base features) R² mean ± SD | 95% range over corpora | RMSE (mg g⁻¹) |
|---|---|---|---|
| Random 5-fold | 0.892 ± 0.024 | – | 5.0 |
| Grouped 5-fold | 0.638 | – | 9.1 |
| Leave-one-study-out | 0.631 ± 0.137 | 0.41–0.82 | 9.2 |
| Leave-one-family-out | 0.559 | 0.14–0.79 | 10.0 |
| External corpora (different hyper-priors) | 0.562 | 0.28–0.77 | 9.7 |

- The study-level cluster bootstrap gave mean 95% intervals of 0.44–0.76 for the LOSO R² of a single corpus. One corpus of 30 studies therefore cannot pin transferability down to better than about ±0.15 in R².
- The same ordering held for XGBoost (random 0.873, LOSO 0.596, LOFO 0.515).
- The Gaussian process collapsed outside random splits (LOSO 0.26 base, 0.33 chem; 3 corpora).

**Apparent model superiority is not significant.**
- RF exceeded XGBoost under LOSO by a mean of 0.038 R² (SD 0.090). It was better in 14 of 20 corpora (sign test p = 0.12).
- Nested grouped tuning changed LOSO R² by +0.001 (RF) and −0.040 (XGBoost) on three corpora. Tuning on 24–29 training studies overfits the inner folds.

**Speciation features do not improve prediction.** With the concentration-, I- and T-dependent features, the change in LOSO R² relative to the base set was:
- RF: −0.009 ± 0.031 (positive in 7 of 20 corpora);
- XGBoost: −0.013 ± 0.060 (positive in 8 of 20 corpora);
- Gaussian process: +0.08 ± 0.49 (3 corpora).

The hypothesis that chemistry-informed features improve predictive reliability (RQ1 of v0.1) is not supported. The speciation fractions are smooth functions of pH, C₀, I and T, which are already inputs, so tree models gain nothing from them (§4.3).

**Most error on a new study is between-study.**
- Forty-five per cent (ICC 0.45 ± 0.14) of the variance of out-of-study residuals was a study-level offset.
- Per-study RMSE under LOSO ranged from 1.9 to 14.8 mg g⁻¹ (10th–90th percentile) and scaled with material capacity: median 9.0 mg g⁻¹ for chemically activated carbons, 3.5 mg g⁻¹ for hydrochars.

**Retransformation.** Duan smearing reduced the mean back-transformation bias from −1.86 to −1.30 mg g⁻¹, but did not remove it, because the residuals are heteroscedastic.

**More studies help, more rows do not.**
- LOSO R² rose from 0.47 ± 0.14 with 12 studies to 0.59 ± 0.14 (24), 0.56 ± 0.11 (48) and 0.75 ± 0.08 (96).
- Random-split R² stayed between 0.80 and 0.90 throughout.
- The effective sample size for transferability is the number of independent studies. On this evidence, a literature corpus of 20–30 studies yields LOSO R² of about 0.5–0.6.

**Table 5.** Out-of-fold R² for apparent Cr(VI) uptake: mean (2.5–97.5% range across 20 corpora; GPR: 3 corpora).

| Model | Features | Random 5-fold | Grouped 5-fold | LOSO | LOFO | External |
|---|---|---|---|---|---|---|
| RF | base | 0.892 (0.85–0.93) | 0.638 (0.40–0.81) | 0.631 (0.41–0.82) | 0.559 (0.14–0.79) | 0.562 (0.28–0.77) |
| RF | chem | 0.882 (0.81–0.92) | 0.629 (0.40–0.79) | 0.622 (0.40–0.79) | 0.534 (0.08–0.74) | 0.525 (0.19–0.75) |
| XGB | base | 0.873 (0.76–0.93) | 0.624 (0.42–0.80) | 0.596 (0.36–0.78) | 0.515 (0.16–0.73) | 0.561 (0.40–0.74) |
| XGB | chem | 0.865 (0.76–0.93) | 0.607 (0.39–0.76) | 0.584 (0.29–0.78) | 0.509 (0.17–0.74) | 0.497 (0.17–0.69) |
| GPR | base | 0.511 (0.03–0.82) | – | 0.256 (−0.39–0.64) | – | – |
| GPR | chem | 0.736 (0.68–0.83) | – | 0.332 (0.20–0.47) | – | – |

![Fig. 4](figures/fig4_ml_validation.png)

**Fig. 4.** Transferability of ML models (synthetic; 20 corpora of 30 studies). (a) R² for apparent Cr(VI) uptake under five validation schemes (bars: mean; whiskers: 2.5–97.5% range across corpora). (b) Learning curve: R² against the number of independent studies (RF with speciation features; mean ± SD of 4 replicates). (c) Per-study RMSE under leave-one-study-out by family (ACc/ACp: chemically/physically activated carbon; BCw/BCm: wood/manure biochar; HC: hydrochar; PC: biomass-derived porous carbon).
