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

The hypothesis that chemistry-informed features improve predictive reliability is not supported. The speciation fractions are smooth functions of pH, C₀, I and T, which are already inputs, so tree models gain nothing from them (§4.3).

**Most error on a new study is between-study.**
- Forty-five per cent (ICC 0.45 ± 0.14) of the variance of out-of-study residuals was a study-level offset.
- Per-study RMSE under LOSO ranged from 1.9 to 14.8 mg g⁻¹ (10th–90th percentile) and scaled with material capacity: median 9.0 mg g⁻¹ for chemically activated carbons, 3.5 mg g⁻¹ for hydrochars.

**The conclusions hold across mechanisms, but transfer between mechanisms is not free (Fig. 5a; Table S6).** RF with speciation features was evaluated on five corpora for each mechanism scenario:
- Random-split R² was 0.86–0.89 for every scenario. LOSO R² was 0.48–0.65 and LOFO R² 0.40–0.63, so the random-split inflation is not specific to scenario E.
- Speciation features again lowered LOSO R² slightly in every scenario (−0.010 to −0.033).
- Models trained on scenario-E corpora and tested on corpora generated under other mechanisms reached R² = 0.65 (A), 0.53 (B), 0.60 (C) and 0.55 (D). For B and D this is 0.05–0.10 below the within-scenario LOSO R²; for A and C it is equal or higher.
- These mechanism shifts were tested only within one model family; they say nothing about mechanisms outside it.

**Retransformation.** Duan smearing reduced the mean back-transformation bias from −1.86 to −1.30 mg g⁻¹, but did not remove it, because the residuals are heteroscedastic.

**More studies help, more rows do not (Fig. 5b).** The learning curve was computed with 8 replicate corpora per size. An earlier run with 4 replicates (Fig. 4b) was not monotonic (0.47, 0.59, 0.56, 0.75):
- LOSO R² rose from 0.31 ± 0.09 (mean ± SE) with 12 studies to 0.58 ± 0.05 (24), 0.66 ± 0.02 (48) and 0.66 ± 0.02 (96). The rank correlation between corpus size and LOSO R² was ρ = 0.58 (p < 0.001, 32 corpora).
- With 8 replicates the curve is monotonic in its means but levels off beyond about 48 studies. The 48- and 96-study means differ by less than their SE. The dip at 48 studies in the 4-replicate run and its higher value at 96 were therefore replicate noise.
- The between-replicate SD fell from 0.26 (12 studies) to 0.06 (96). Small corpora are unreliable as well as less accurate.
- Random-split R² stayed between 0.85 and 0.90 throughout.
- The effective sample size for transferability is the number of independent studies. On this evidence, a literature corpus of 20–30 studies yields LOSO R² of about 0.5–0.6, and the plateau near 0.66 reflects irreducible between-study heterogeneity in the simulator, not a general property of real data.

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

**Fig. 4.** Transferability of ML models (synthetic; 20 corpora of 30 studies). (a) R² for apparent Cr(VI) uptake under five validation schemes (bars: mean; whiskers: 2.5–97.5% range across corpora). (b) Learning curve with 4 replicates: R² against the number of independent studies (RF with speciation features; mean ± SD of 4 replicates); superseded by the 8-replicate curve in Fig. 5b. (c) Per-study RMSE under leave-one-study-out by family (ACc/ACp: chemically/physically activated carbon; BCw/BCm: wood/manure biochar; HC: hydrochar; PC: biomass-derived porous carbon).

![Fig. 5](figures/fig8_mechanism_transfer.png)

**Fig. 5.** Robustness of the ML and interval results across mechanisms (synthetic). (a) Out-of-fold R² for apparent uptake by mechanism scenario (RF with speciation features; 5 corpora per scenario; bars: mean, whiskers: SD): random 5-fold, leave-one-study-out, leave-one-family-out, and models trained on scenario-E corpora and tested on the given scenario. (b) Learning curve: leave-one-study-out R² against the number of studies (8 replicate corpora per size; mean ± SE; points: replicates) and random-split R². (c) Fraction of unseen studies with ≥ 90% within-study coverage at nominal 90%, by scenario and interval method (dotted: 0.90).
