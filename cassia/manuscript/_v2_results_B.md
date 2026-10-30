### 3.5 Study-level prediction intervals (synthetic)

Table 6 and Fig. 5 report outer leave-one-study-out evaluations on 10 corpora: 300 test studies for Cr(VI) endpoints and 182 for total-Cr endpoints, which are measured only in Class A/B studies.

**Average coverage is the wrong criterion.**
- The pooled-residual interval M0 (the v0.1 method) reached its nominal *average* coverage for every endpoint (0.89–0.90 at 90%).
- Only 65–75% of individual test studies were covered at 90%, and the worst study had 11–29% of its observations covered.
- The observations of one study share an offset, so pooled calibration averages over-coverage of 'easy' studies with severe under-coverage of others.

**The study-level CV+ intervals address exactly this.**
- M1 (study-weighted) raised the share of studies covered at ≥ 90% to 76–87%.
- M2 (one score per study) raised it to 88–93%. That exceeds the 1 − 2α = 0.80 lower bound implied by its construction.
- M2 is conservative at lower nominal levels (mean coverage 0.77–0.82 at a nominal 0.50).

**The price is width, and it is large.**
- At 90%, M2 intervals spanned 2.1–3.1 log₁₀ units, against 1.3–1.7 for M0: a factor of 100–1,000 between the bounds.
- At 95%, intervals for the total-Cr endpoints were effectively unbounded. Only about 17 studies with total-Cr data were available for calibration, fewer than the 19 that a study-level 95% quantile requires.
- This is a concrete data requirement: certifying total-Cr predictions at 95% needs at least 19 independent studies that report total Cr.

**Coverage is marginal over families, not conditional.** For apparent uptake, M0 covered 0.77 of hydrochar observations (0.85 with M1, 0.91 with M2), compared with 0.96–0.99 for the other families. Study-level calibration does not repair heterogeneity between families.

**Table 6.** Study-level coverage at a nominal 90% under leave-one-study-out validation (10 corpora; RF with speciation features). Coverage is the fraction of a test study's observations inside its interval.

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

![Fig. 5](figures/fig5_conformal.png)

**Fig. 5.** Study-level conformal prediction (synthetic; outer leave-one-study-out on 10 corpora). (a) Nominal against mean per-study coverage for apparent uptake (dotted: 10th-percentile study). (b) Distribution of per-study coverage at 90%. (c) Mean per-study coverage at 90% for the four endpoints. M0: pooled residuals; M1: study-level CV+ with study-weighted scores; M2: study-level CV+ with one score per study.

### 3.6 Treatment design with two endpoints (synthetic)

**Set-up.** Of the 300 new materials per target (20 cohorts × 15), 169 (S1) and 181 (S2) could truly meet both limits somewhere on the design grid (Table 7, Fig. 6).

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
- It accepted designs for 105 (S1) and 98 (S2) materials. Every accepted design complied in truth (Wilson 95% lower bounds 0.965 and 0.962), with no unsafe recommendation.
- It found compliant designs for 62% (S1) and 54% (S2) of attainable materials, after a mean of 1.4–1.5 rounds.
- The median dose was twice the true minimum dose, which is the price of escalating on a coarse dose grid.

**Sensitivity of accepted designs (S1).**
- They tolerated a dose increase of 20% (100% compliance) and a pH increase of 0.25 (95%).
- They were sensitive to lower pH: compliance fell to 81% at −0.25 and 62% at −0.5 pH units, because dissolved Cr(III) release rises as pH falls.
- A 20% dose reduction gave 83% compliance.

**Effluent predictions.** Predicted and true effluent concentrations on a design sub-grid were rank-correlated (Spearman 0.65 for Cr(VI), 0.55 for total Cr), with median absolute errors of 0.40 and 0.22 log₁₀ units. The 95% upper bounds held for 97.9% (Cr(VI)) and 100% (total Cr) of grid points.

**Table 7.** Design outcomes for 300 new virtual materials per target (S1: C₀ = 20 mg L⁻¹, limits Cr(VI) ≤ 0.5 and total Cr ≤ 2 mg L⁻¹; S2: C₀ = 5 mg L⁻¹, limits 0.1 and 0.5 mg L⁻¹).

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

![Fig. 6](figures/fig6_design.png)

**Fig. 6.** Treatment design for new virtual materials (synthetic; 20 cohorts × 15 materials). (a) Fraction of recommended designs that comply in truth, with Wilson 95% intervals; numbers in the bars are median doses. (b) Predicted against true total Cr at the selected designs (S1). (c) Compliance of point + verification designs under perturbations of pH and dose (S1).

### 3.7 Propagation of unreported chemistry

**Prediction level.** Unknown ionic strength (10⁻³–10⁻¹ M) and an unknown pH basis (±0.5) were propagated with 1,000 draws per observation for 462 pH-edge observations from studies that did not report them.
- The 90% range of f(HCrO₄⁻) was 0.035 at pH ≤ 3.5, 0.045 at pH 3.5–5, 0.24 at pH 5–6.5 and 0.08 at pH 6.5–8.5.
- The resulting spread of the point prediction was 0.08–0.14 log₁₀ units, small compared with the study-level interval width (1.18–1.28 log₁₀ units).
- The union of the M1 intervals over all scenarios was 6–14% wider than the nominal interval.

**Design level.** Twenty-four new materials were each re-designed 300 times with the test-confirmed protocol, under different assumed ionic strengths.
- An unknown ionic strength changed the accepted dose for 3 of the 10 materials with accepted designs (e.g. 4–8 g L⁻¹).
- It never produced a non-compliant accepted design, because the confirmatory batches are run at the true conditions.

**Mechanistic level.** Unreported pH matters more than unreported ionic strength (§3.3): a 0.5-unit pH error biased the inferred reduction share by up to 0.16 (90th percentile).

**Conclusion.** For prediction and design, missing solution chemistry is a second-order uncertainty compared with between-study heterogeneity. For mechanistic interpretation, a missing pH basis is first-order.

![Fig. 7](figures/fig7_missing_chemistry.png)

**Fig. 7.** Propagation of unreported ionic strength and pH basis (1,000 draws per observation; synthetic studies). (a) 90% range of f(HCrO₄⁻). (b) Width of the nominal M1 interval, of its union over scenarios, and spread of the point prediction (log₁₀). (c) Accepted dose under unknown ionic strength for each new material (box over 300 draws; numbers: compliance probability).

### 3.8 Empirical application: status

**Not done.** No published experimental observation was analysed in this revision. The seed dataset (Ramirez et al., 2020; Data in Brief 29, 105292) and the systematic search were blocked: the hosts required to retrieve articles and supplementary files (PubMed Central, ScienceDirect, Crossref, Europe PMC, Zenodo, Mendeley Data, Figshare) were not reachable from the computational environment. No values were transcribed from memory or estimated.

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
1. Extract the Ramirez et al. data as a Class C benchmark of the workflow.
2. Carry out the registered search and extraction, identifying the subset of studies with total-Cr and XPS/XANES evidence.
3. Apply the identifiability and validation analyses of §3.3–3.6 to the real corpus, keeping the evidence classes separate.
