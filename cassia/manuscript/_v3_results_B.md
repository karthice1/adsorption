### 3.5 Study-level prediction intervals (synthetic)

Table 7 and Fig. 6 report outer leave-one-study-out evaluations on 10 corpora: 300 test studies for Cr(VI) endpoints and 182 for total-Cr endpoints, which are measured only in Class A/B studies.

**Average coverage is the wrong criterion.**
- The pooled-residual interval M0 (the conventional approach) reached its nominal *average* coverage for every endpoint (0.89–0.90 at 90%).
- Only 65–75% of individual test studies were covered at 90%, and the worst study had 11–29% of its observations covered.
- The observations of one study share an offset, so pooled calibration averages over-coverage of 'easy' studies with severe under-coverage of others.

**The study-level CV+ intervals address exactly this.**
- M1 (study-weighted) raised the share of studies covered at ≥ 90% to 76–87%.
- M2 (one score per study, CV+) raised it to 88–93%. M2 has no proven guarantee (Appendix A), so these figures are empirical only.
- M2 is conservative at lower nominal levels (mean coverage 0.77–0.82 at a nominal 0.50).

**The price is width, and it is large.**
- At 90%, M2 intervals spanned 2.1–3.1 log₁₀ units, against 1.3–1.7 for M0: a factor of 100–1,000 between the bounds.
- At 95%, intervals for the total-Cr endpoints were effectively unbounded. Only about 17 studies with total-Cr data were available for calibration, fewer than the 19 that a study-level 95% quantile requires.
- This is a concrete data requirement: certifying total-Cr predictions at 95% needs at least 19 independent studies that report total Cr.

**The exact split-conformal interval M2s behaved as its theorem predicts (Table 7b, Fig. 6).**
- *Calibration.* Training studies were split in half, giving 14 calibration studies for Cr(VI) endpoints and 8.6 on average for total-Cr endpoints. The fraction of unseen studies with at least 90% of their observations covered tracked the nominal level for apparent uptake: 0.93, 0.80, 0.72 and 0.50 at nominal 0.90, 0.80, 0.70 and 0.50. Every value lies within the theoretical bounds of 1 − α and 1 − α + 1/(J + 1).
- *Cost.* The price of exactness is width, 3.3 log₁₀ units at 90% for apparent uptake. Only half of the training studies are used to fit the model.
- *Unbounded intervals.* At 95%, M2s is unbounded for every endpoint, because it needs J ≥ 19. At 90% it is unbounded for 63% of test studies on the total-Cr endpoints, because it needs J ≥ 9.
- *Across mechanisms.* The guarantee held in all five scenarios (≥ 0.92 of studies covered at ≥ 90%; Fig. 5c).
- *Interpretation.* Exact study-level guarantees are available, but with literature-scale numbers of studies they are too wide for low-effluent engineering decisions. This is a substantive negative result, not a defect of the method.

**Table 7b.** M2s, exact study-level split conformal: fraction of unseen studies with ≥ 90% within-study coverage (10 corpora, leave-one-study-out).

| Endpoint | Calibration studies (mean) | 1 − α = 0.50 | 0.70 | 0.80 | 0.90 | 0.95 | Width at 0.90 (log₁₀; bounded intervals) |
|---|---|---|---|---|---|---|---|
| q_app,VI | 14 | 0.50 | 0.72 | 0.80 | 0.93 | unbounded | 3.28 |
| C_e,VI | 14 | 0.57 | 0.76 | 0.81 | 0.95 | unbounded | 3.94 |
| q_tot,aq | 8.6 | 0.52 | 0.78 | 0.88 | 0.97 (63% unbounded) | unbounded | 3.81 |
| C_e,tot | 8.6 | 0.52 | 0.72 | 0.87 | 0.98 (63% unbounded) | unbounded | 2.18 |

**Coverage is marginal over families, not conditional.** For apparent uptake, M0 covered 0.77 of hydrochar observations (0.85 with M1, 0.91 with M2), compared with 0.96–0.99 for the other families. Under M2s, only 84% of hydrochar studies (42 of 50; Wilson 95% interval 72–92%) met the 90% within-study target, against 90–98% for the other families (Table S5b; mean within-study coverage 0.95 against 0.97–0.995). With 50 test studies per family the intervals of the families overlap, so this is a tendency, not a demonstrated family effect; it is consistent with the M0–M1 results (Table S5) and with the guarantee being marginal, not family-conditional. Study-level calibration does not repair heterogeneity between families.

**Table 7.** Study-level coverage at a nominal 90% under leave-one-study-out validation (10 corpora; RF with speciation features). Coverage is the fraction of a test study's observations inside its interval.

| Endpoint | Method | Mean | Median | 10th pct | Worst | Studies ≥ 0.90 | Width (log₁₀) | Calibration studies |
|---|---|---|---|---|---|---|---|---|
| q_app,VI | M0 | 0.898 | 0.96 | 0.71 | 0.11 | 0.70 | 1.40 | 29 |
| q_app,VI | M1 | 0.937 | 1.00 | 0.82 | 0.11 | 0.80 | 1.76 | 29 |
| q_app,VI | M2 | 0.960 | 1.00 | 0.89 | 0.21 | 0.88 | 2.32 | 29 |
| q_tot,aq | M0 | 0.893 | 1.00 | 0.58 | 0.15 | 0.75 | 1.66 | 17 |
| q_tot,aq | M1 | 0.953 | 1.00 | 0.89 | 0.15 | 0.87 | 2.27 | 17 |
| q_tot,aq | M2 | 0.976 | 1.00 | 0.90 | 0.40 | 0.91 | 3.11 | 17 |
| C_e,VI | M0 | 0.900 | 1.00 | 0.64 | 0.29 | 0.65 | 1.61 | 29 |
| C_e,VI | M1 | 0.940 | 1.00 | 0.79 | 0.39 | 0.76 | 1.94 | 29 |
| C_e,VI | M2 | 0.974 | 1.00 | 0.89 | 0.46 | 0.90 | 2.68 | 29 |
| C_e,tot | M0 | 0.901 | 1.00 | 0.70 | 0.11 | 0.73 | 1.28 | 17 |
| C_e,tot | M1 | 0.949 | 1.00 | 0.85 | 0.18 | 0.87 | 1.68 | 17 |
| C_e,tot | M2 | 0.973 | 1.00 | 0.90 | 0.25 | 0.93 | 2.07 | 17 |

![Fig. 6](figures/fig5_conformal.png)

**Fig. 6.** Study-level conformal prediction (synthetic; outer leave-one-study-out on 10 corpora). (a) Nominal against mean per-study coverage for apparent uptake (dotted: 10th-percentile study). (b) Distribution of per-study coverage at 90%. (c) Mean per-study coverage at 90% for the four endpoints. M0: pooled residuals; M1: study-level CV+ with study-weighted scores; M2: study-level CV+ with one score per study.

### 3.6 Treatment design with two endpoints (synthetic)

**Set-up.** Of the 300 new materials per target (20 cohorts × 15), 169 (S1) and 181 (S2) could truly meet both limits somewhere on the design grid (Table 8, Fig. 7).

**Point predictions are unsafe.**
- They recommended a design for 105 materials (S1). Only 41% of these complied in truth.
- Sixty-two materials (20.7%; Wilson 95% interval 16.5–25.6%) received an unsafe recommendation: 38 failed the Cr(VI) limit and 45 the total-Cr limit (some failed both).
- For S2 the unsafe rate was 7.0% (4.6–10.5%).

**Calibrated bounds are safe but useless.**
- Marginal-90% and joint-90% rules recommended nothing for any material, so they never failed and never helped.
- With the study-level uncertainty available from 30 heterogeneous training studies, no effluent at or below 0.5 mg L⁻¹ Cr(VI) can be certified for a new material.
- This is not a defect of the rule. It is the honest consequence of the between-study variability in §3.4.

**Confirmatory batches help only partly.**
- Estimating the new material's offset from four confirmatory batches (pH 3, 0.5–4 g L⁻¹) narrowed the bounds.
- At joint 90% (95% per endpoint) the rule recommended a design for 1 material only, because the total-Cr bound remains unbounded with about 17 calibration studies (§3.5).
- At joint 80% it recommended designs for 14 (S1) and 6 (S2) materials, with 86% and 83% compliance and 2 and 1 unsafe recommendations.

**Model-guided, test-confirmed design worked.**
- In this protocol the model proposes the pH and starting dose. Duplicate confirmatory batches are then run at the proposed design, and the dose is escalated until both measured concentrations fall below 90% of their limits (at most four rounds).
- These designs were obtained with the **simple acceptance rule** (duplicates; mean ≤ 90% of each limit), with one realised verification sequence per material.
- It accepted designs for 105 (S1) and 98 (S2) of 300 materials per target. Every accepted design complied in truth: no unsafe accepted designs were observed in the simulated test cohorts.
- The upper confidence bound depends on the denominator (Table 8c). Over all 300 materials per target, the Wilson 95% upper bound on the unsafe rate is 1.3%. Among accepted designs only (0/105 and 0/98), it is 3.5% (S1) and 3.8% (S2). These are observed frequencies, not guarantees.
- It found compliant designs for 62% (S1) and 54% (S2) of attainable materials, after a mean of 1.4–1.5 rounds.
- The median dose was twice the true minimum dose, which is the price of escalating on a coarse dose grid.

**Verification is a sequential decision, so its false-acceptance risk was evaluated explicitly (Table 8b, Fig. 7b,c).**
- *Simulation.* An independent cohort of 300 new materials per target was drawn from the same hyper-priors as the design cohort (20 training corpora × 15 materials; the material draws differ from those of the design stage). For every material, target and rule, 1,000 verification sequences were simulated, with up to four dose rounds and both endpoints measured each round. Each sequence drew a calibration bias and correlated replicate errors.
- *Per material.* The probability of accepting a non-compliant design was:
  - simple rule (duplicates, mean ≤ 90% of each limit): mean 0.01% across materials, worst material 1.4–1.7%;
  - error-controlled rule (triplicates, one-sided upper confidence bound ≤ limit, Bonferroni over two endpoints and four rounds): worst material 0.1–0.7%.
- *Overall.* No material exceeded 5% under either rule, with shared or independent calibration bias. The controlled rule needed about 50% more batches per material (4.7 against 3.1 for S1; 4.0 against 2.6 for S2).
- *Operating characteristic.* Because these averages depend on how close true concentrations sit to the limits, Table 8b gives the acceptance probability over four rounds for a design whose true concentrations stay at a fixed ratio to both limits:

| True concentration / limit | Simple rule | Error-controlled rule |
|---|---|---|
| 1.00 (at the limit) | 3.7% | 0.8% |
| 1.05 | 0.2% | 0.04% |
| 0.90 (truly compliant) | 59% | 37% |

  Both figures assume shared calibration bias. The controlled rule is about five times safer at the limit and less likely to accept compliant designs close to the limit.
- *Scope.* These are properties under the assumed analytical-error model (3% random, 5% calibration uncertainty), not guarantees for real laboratories.

**Table 8b.** Sequential verification: probability of accepting a non-compliant design (1,000 simulated sequences per material; 300 materials per target).

| Target | Rule | Bias model | Mean P(false acceptance) | Worst material | Materials > 5% | Mean batches |
|---|---|---|---|---|---|---|
| S1 | simple (n = 2, mean ≤ 0.9·limit) | shared | 0.0001 | 0.014 | 0 | 3.1 |
| S1 | controlled (n = 3, UCB ≤ limit) | shared | 0.0000 | 0.001 | 0 | 4.7 |
| S2 | simple | shared | 0.0001 | 0.015 | 0 | 2.6 |
| S2 | controlled | shared | 0.0000 | 0.007 | 0 | 4.0 |
| S1 | simple | independent | 0.0001 | 0.015 | 0 | 3.1 |
| S1 | controlled | independent | 0.0000 | 0.002 | 0 | 4.7 |

**Table 8c.** Accepted designs and unsafe acceptances by rule, with explicit denominators (synthetic; 300 materials per target in each cohort). Design cohort: one realised sequence per material; verification cohort: expected numbers over 1,000 sequences per material (shared calibration bias).

| Cohort | Target | Rule | Accepted designs | Unsafe accepted | Unsafe rate, all materials (Wilson 95% upper) | Unsafe rate among accepted (Wilson 95% upper) |
|---|---|---|---|---|---|---|
| Design | S1 | simple | 105 | 0 | 0/300 (1.3%) | 0/105 (3.5%) |
| Design | S2 | simple | 98 | 0 | 0/300 (1.3%) | 0/98 (3.8%) |
| Verification | S1 | simple | 82.5 | 0.018 | 0.006% | 0.022% |
| Verification | S1 | error-controlled | 81.6 | 0.002 | 0.001% | 0.003% |
| Verification | S2 | simple | 100.2 | 0.036 | 0.012% | 0.036% |
| Verification | S2 | error-controlled | 98.4 | 0.013 | 0.004% | 0.013% |

The error-controlled rule accepted almost as many designs as the simple rule (81.6 against 82.5 for S1; 98.4 against 100.2 for S2), at about 50% more batches, and reduced the expected number of unsafe acceptances three- to ninefold. The numbers of accepted designs differ between the two cohorts because they contain different materials. Wilson bounds are given only for observed counts; the verification-cohort rates are model expectations under the assumed analytical-error model.

**Sensitivity of accepted designs (S1).**
- They tolerated a dose increase of 20% (100% compliance) and a pH increase of 0.25 (95%).
- They were sensitive to lower pH: compliance fell to 81% at −0.25 and 62% at −0.5 pH units, because dissolved Cr(III) release rises as pH falls.
- A 20% dose reduction gave 83% compliance.

**Out-of-study accuracy of the endpoint models (Table S13).** The design layer uses four RF models: Cr(VI) and total Cr, as uptake and as effluent concentration. The total-Cr models are trained only on rows with a total-Cr measurement, from about 18 Class A/B studies per corpus. Under leave-one-study-out (10 corpora), the R² on the log scale was:

| Endpoint | LOSO R² (log₁₀) | Random-split R² (log₁₀) | LOSO RMSE (log₁₀) |
|---|---|---|---|
| C_e,VI | 0.74 ± 0.06 | 0.94 | 0.49 |
| C_e,tot | 0.65 ± 0.09 | 0.91 | 0.37 |
| q_app,VI | 0.39 ± 0.16 | 0.77 | 0.51 |
| q_tot,aq | 0.42 ± 0.17 | 0.75 | 0.52 |

R² on the log scale weights low concentrations, which matter for compliance, more than the linear-scale R² of Table 6 (linear-scale LOSO R² here: 0.87, 0.85, 0.56 and 0.52). The total-Cr models are less accurate than the Cr(VI) models, as expected from fewer training studies, and both lose about 0.2–0.4 in R² between random and study-level validation. Their study-level intervals are reported in Table 7.

**Effluent predictions.** Predicted and true effluent concentrations on a design sub-grid were rank-correlated (Spearman 0.65 for Cr(VI), 0.55 for total Cr), with median absolute errors of 0.40 and 0.22 log₁₀ units. The 95% upper bounds held for 97.9% (Cr(VI)) and 100% (total Cr) of grid points.

**Table 8.** Design outcomes for 300 new virtual materials per target (S1: C₀ = 20 mg L⁻¹, limits Cr(VI) ≤ 0.5 and total Cr ≤ 2 mg L⁻¹; S2: C₀ = 5 mg L⁻¹, limits 0.1 and 0.5 mg L⁻¹).

| Target | Rule | Designs recommended | Compliant (true) | Unsafe (of all materials) | Attainable materials found | Median dose (g L⁻¹) |
|---|---|---|---|---|---|---|
| S1 | point prediction | 105 | 41% | 62 (20.7%) | 25% | 4 |
| S1 | marginal 90% | 0 | – | 0 | 0% | – |
| S1 | joint 90% (Bonferroni) | 0 | – | 0 | 0% | – |
| S1 | joint 90% + 4 confirmatory batches | 1 | 100% | 0 | 0.6% | 6 |
| S1 | joint 80% + 4 confirmatory batches | 14 | 86% | 2 (0.7%) | 7% | 4 |
| S1 | point + verification batches | 105 | 100% | 0 | 62% | 8 |
| S2 | point prediction | 39 | 46% | 21 (7.0%) | 10% | 4 |
| S2 | joint 90% (Bonferroni) | 0 | – | 0 | 0% | – |
| S2 | joint 80% + 4 confirmatory batches | 6 | 83% | 1 (0.3%) | 3% | 4 |
| S2 | point + verification batches | 98 | 100% | 0 | 54% | 8 |

![Fig. 7](figures/fig6_design.png)

**Fig. 7.** Treatment design for new virtual materials (synthetic; 20 cohorts × 15 materials per target). (a) Fraction of all materials receiving an unsafe recommendation, and fraction of attainable materials for which a compliant design was found, by design rule (point prediction; marginal and joint 90% bounds; joint 90% and 80% bounds with up to four confirmatory batches; point prediction with sequential verification); solid: S1, hatched: S2. (b) Distribution over materials of the probability of accepting a non-compliant design under the two verification rules (shared calibration bias; materials whose dose ladder contains a non-compliant dose; dotted: 5%). (c) Mean probability of acceptance within four rounds; numbers in the bars: mean batches per material.

### 3.7 Propagation of unreported chemistry

**Prediction level.** Unknown ionic strength (10⁻³–10⁻¹ M) and an unknown pH basis (±0.5) were propagated with 1,000 draws per observation for 462 pH-edge observations from studies that did not report them.
- The 90% range of f(HCrO₄⁻) was 0.035 at pH ≤ 3.5, 0.045 at pH 3.5–5, 0.24 at pH 5–6.5 and 0.08 at pH 6.5–8.5.
- The resulting spread of the point prediction was 0.08–0.14 log₁₀ units, small compared with the study-level interval width (1.18–1.28 log₁₀ units).
- The union of the M1 intervals over all scenarios was 6–14% wider than the nominal interval.

**Design level.** Twenty-four new materials were each re-designed 300 times with the test-confirmed protocol, under different assumed ionic strengths.
- An unknown ionic strength changed the accepted dose for 3 of the 10 materials with accepted designs (e.g. 4–8 g L⁻¹).
- It never produced a non-compliant accepted design, because the confirmatory batches are run at the true conditions.
- The 95th percentile of true residual total Cr at the accepted designs was 0.08–1.38 mg L⁻¹ across materials, below the 2 mg L⁻¹ limit.

**Capacity level.** A conventional isotherm fit uses neither ionic strength nor the pH basis, so unreported solution chemistry does not change the fitted q_max itself. It changes what q_max means: the inferred reduction share and the speciation-based affinity (mechanistic level below).

**Mechanistic level.** Unreported pH matters more than unreported ionic strength (§3.3): a 0.5-unit pH error biased the inferred reduction share by up to 0.16 (90th percentile).

**Conclusion.** For prediction and design, missing solution chemistry is a second-order uncertainty compared with between-study heterogeneity. For mechanistic interpretation, a missing pH basis is first-order.

![Fig. 8](figures/fig7_missing_chemistry.png)

**Fig. 8.** Propagation of unreported ionic strength and pH basis (1,000 draws per observation; synthetic studies). (a) 90% range of f(HCrO₄⁻). (b) Width of the nominal M1 interval, of its union over scenarios, and spread of the point prediction (log₁₀). (c) Accepted dose under unknown ionic strength for each new material (box over 300 draws; numbers: compliance probability).

### 3.8 Empirical application: status

**Not done.** No published experimental observation was analysed. The seed dataset (Ramirez et al., 2020; Data in Brief 29, 105292) is openly available through PubMed Central (PMC7044665) with its raw isotherm and kinetic files as supplementary material. It reports Cr(VI) removal by a ZnCl₂-activated teakwood-sawdust carbon; its interpretation is in a companion article (J. Environ. Chem. Eng. 8, 103702; doi:10.1016/j.jece.2020.103702). It measures residual Cr(VI) only and is therefore a Class C benchmark: it can test whether CASSIA processes genuine observations, not whether it recovers the true Cr(VI)/Cr(III) distribution. Neither this dataset nor the articles of the systematic search could be retrieved, because the hosts that serve them (PubMed Central, NCBI, ScienceDirect, doi.org, Crossref, Europe PMC, Mendeley Data, Zenodo, Figshare and the institutional repository) were not reachable from the computational environment. No values were transcribed from memory, estimated or digitised from figures.

**Ready.** The empirical pipeline is complete and tested:
- schema validation and QC flags;
- observation-level evidence classes;
- endpoints with Class C labelling;
- isotherm selection with bootstrap intervals;
- the aqueous Cr(III) index;
- study-level CV+ intervals;
- the design protocol.

It runs end to end on files in the CASSIA schema (`examples/user_data_workflow.py`).

**Next steps, in order:**
1. Extract the Ramirez et al. data (isotherm and kinetic series with pH, C₀, dose, volume, temperature and contact time) and run the full workflow as a Class C benchmark.
2. Extract at least 20 further Cr(VI)–carbon studies by the registered protocol (`data/PROTOCOL.md`), reconstructing pH (initial and final), C₀, dose, solution volume, temperature, contact time, background electrolyte and ionic strength.
3. Classify every observation (C, B, A1, A2), identifying the studies that report total Cr and those with XPS or XANES evidence.
4. Apply the equivalence and identifiability tests (§3.3), the LOSO/LOFO validation (§3.4) and the study-level intervals (§3.5) to the real corpus, keeping evidence classes separate.
5. Report empirical findings in a section separate from the synthetic benchmark.
