# CASSIA: A Provenance-Aware Computational Framework for Evidence-Constrained Analysis of Cr(VI) Removal by Carbonaceous Materials

## Highlights

- Residual Cr(VI) data seldom rule out adsorption-only fits when reduction dominates
- Total-Cr and XPS data make reduction detectable, defining evidence classes B and A
- Three PHREEQC databases differ by 63–83 mV in the Cr(VI)/Cr(III) redox boundary
- Study-level conformal intervals cover 88–93% of new studies but span ≥2 log units
- Model-guided, test-confirmed design avoided all unsafe effluent recommendations

## Abstract

Loss of dissolved Cr(VI) on carbonaceous materials is routinely reported as adsorption, although Cr(VI) is also reduced to Cr(III). CASSIA is an open, provenance-aware framework that keeps apparent removal, total-aqueous-Cr removal and adsorption apart, assigns observation-level evidence classes, and asks which mechanisms a given measurement design can identify. Three PHREEQC databases agreed on Cr(VI) speciation (within 0.044 in species fraction) but differed by 63–83 mV in the Cr(VI)/Cr(III) boundary; Cr(III) remained undersaturated up to 5,000 mg L⁻¹ at pH ≤ 3.5–4.0. In known-truth simulations (five mechanism scenarios, 15,000 virtual studies), isotherm fits of residual Cr(VI) with R² ≥ 0.95 occurred in 41–52% of studies with more than 20% reduction under surface, adsorption-coupled and combined reduction, and Langmuir q_max then exceeded the adsorption capacity in proportion to the reduction share. A flexible adsorption-only model could not be rejected for 79% of high-reduction studies observed through residual Cr(VI) alone; adding total Cr and XPS reduced this to 32% and 12.5%. Across 20 replicate corpora, random cross-validation overstated transferability (random-forest R² 0.89 versus 0.63 leave-one-study-out), speciation features did not improve prediction, and 45% of residual variance was between-study. Study-level CV+ intervals covered at least 90% of observations in 88–93% of unseen studies but were too wide to certify low effluent concentrations; a model-guided, test-confirmed design protocol found compliant designs for 54–62% of attainable materials with no unsafe recommendation. No empirical literature data were analysed; the extraction protocol and pipeline are released for that step.

**Keywords:** hexavalent chromium; carbonaceous adsorbents; PHREEQC; mechanism identifiability; leave-one-study-out validation; conformal prediction; data provenance

---

## 1. Introduction

Hexavalent chromium is a priority contaminant in effluents from electroplating, leather tanning, pigment and wood-preservation industries. Its toxicity and mobility differ sharply from those of Cr(III) (Kotaś and Stasicka, 2000). Activated carbons, biochars and biomass-derived porous carbons are the most widely studied low-cost sorbents for Cr(VI) (Mohan and Pittman, 2006; Mohan et al., 2011; Dong et al., 2011), and pooled literature data are increasingly used to train machine-learning (ML) models.

**Redox behaviour on carbon.** Carbon surfaces are redox active towards Cr(VI). Experimental studies on biomass and biochars have shown that, at low pH, Cr(VI) is reduced to Cr(III). The reduction can proceed directly from solution or after adsorption, through oxygen-containing functional groups and the carbon matrix. The resulting Cr(III) is partly released to solution, partly retained on the surface and, at higher pH, precipitated (Park et al., 2004; Dong et al., 2011; Mohan et al., 2011).

**What residual Cr(VI) measures.** The quantity most often reported, residual Cr(VI) by the 1,5-diphenylcarbazide method, records the *disappearance* of Cr(VI), not its adsorption. Isotherm constants fitted to q = (C₀ − C_e,VI)/m therefore mix adsorption, reduction and precipitation in unknown proportions. This adds to well-documented problems of aqueous adsorption modelling: linearised fitting, R² as the only selection criterion, and over-interpreted kinetic models (Tran et al., 2017; Simonin, 2016; Wang and Guo, 2020).

**Two further gaps.**
- *Chemistry.* Cr(VI) speciation depends on concentration through dichromate formation, and Cr(III) solubility depends strongly on pH. Both rest on thermodynamic data whose sources differ (Rai et al., 1987; Ball and Nordstrom, 1998). Yet adsorption studies rarely report the background electrolyte, the final pH or the redox state, and ML studies seldom propagate this missing information.
- *Statistics.* Observations from one publication share an adsorbent, a laboratory and a protocol. Random train/test splits therefore leak study identity (Roberts et al., 2017). Prediction intervals calibrated on pooled observations carry no guarantee for a new study, because exchangeability holds, at best, between studies (Dunn et al., 2023).

**Contribution.** This paper does not propose a new removal mechanism or another chromium-uptake predictor. It contributes an open framework, CASSIA (*Chromium Adsorption, Speciation and Surface-reduction Identifiability Analysis*), built to answer three questions:

- **Q1.** What can actually be established about Cr(VI) removal from a given set of measurements?
- **Q2.** What remains mechanistically unidentifiable without additional evidence, and which additional measurement resolves it?
- **Q3.** How much confidence should an engineer place in predicted removal performance, and in a treatment design, for a new carbonaceous material?

**Evidence used.** Two kinds of evidence are used, and they are kept separate throughout:
1. *Thermodynamic calculations* with three chromium-bearing PHREEQC databases.
2. *Known-truth simulation benchmarks.* A multi-mechanism kinetic model generates virtual batch studies in which the true split between adsorption, reduction and precipitation is known. Five alternative mechanism scenarios, at least 100 replicate corpora and randomised family-level priors are used, so that conclusions do not rest on one assumed mechanism.

**Status of empirical data.** The application to published experimental data is defined (schema, evidence classes, extraction protocol, seed dataset) but has not yet been carried out (Section 3.8). No result in this paper is an empirical finding about a real carbon.

## 2. Methods

### 2.1 Endpoints, evidence classes and terminology

Three endpoints are computed and never merged:
- **Apparent Cr(VI) uptake**, q_app,VI = (C₀ − C_e,VI)/D. It is always computable.
- **Total-aqueous-Cr uptake**, q_tot,aq = (C₀ − C_e,tot)/D. It requires a total-Cr analysis (AAS or ICP). It is not adsorption either, because it includes Cr(III) retained or precipitated on the solid.
- **Aqueous Cr(III) index**, Σ(C_e,tot − C_e,VI)/Σ(C₀ − C_e,VI). It is a lower bound on the share of Cr(VI) disappearance due to reduction.

Evidence classes are assigned **per observation** (`cassia.dataset.classify`), because a single paper may report total Cr or XPS for some experiments only:
- **Class A:** Cr(VI) and total Cr measured, solid-phase Cr speciation (XPS/XANES), final or controlled pH, and a reported electrolyte.
- **Class B:** total Cr or solid-phase speciation.
- **Class C:** residual Cr(VI) only.

A study-level label (the best class present) is kept for reporting only. Class C observations may be labelled only as *apparent* Cr(VI) removal.

Negative apparent uptake (C_e,VI > C₀ within analytical error) is retained in the data, flagged, and censored at q = 0 for model targets. Values below the limit of detection (LOD) enter as LOD/2 and are flagged.

### 2.2 Software architecture and release

CASSIA 0.2.0 is organised in six layers (Table 1).

**Table 1.** CASSIA layers, modules and responsibilities.

| Layer | Module | Responsibility |
|---|---|---|
| 1 Data | `cassia.dataset` | 34-field schema with DOI and source location per row; validation and QC flags; endpoints; observation-level evidence classes |
| 2 Chemistry | `cassia.chemistry` | PHREEQC adapter (phreeqpython 1.6.2): Cr(VI) speciation, a 4-D speciation grid, redox scenarios, Cr(III) saturation, database provenance |
| 3 Simulator | `cassia.simulate` | Known-truth multi-mechanism batch model (scenarios A–E), mass and electron balances, measurement model |
| 4 Identifiability | `cassia.identifiability` | Fisher information by measurement design, detection tests, equivalence tests, profile likelihood, metadata-bias propagation |
| 5 Models and uncertainty | `cassia.fitting`, `cassia.validation` | Nonlinear isotherm/kinetic selection; RF, XGBoost and GPR under four validation schemes; study-level conformal prediction |
| 6 Decision support | `scripts/run_v2.py` (design stage) | Two-endpoint dose/pH design under point, marginal and joint rules |
| auxiliary | `cassia.pores` | pyGAPS 4.x adapter: BET area (Rouquerol criteria) from published N₂ isotherms |

**Reuse of existing packages.**
- The isotherm and kinetic equations are taken from the `adsorption` package (Amvro23, MIT licence). Its gas-phase concentration transform and its fixed-temperature Temkin form are bypassed for aqueous data.
- pyGAPS (Iacomi and Llewellyn, 2019) is used only for gas-isotherm descriptors.

**Release contents.** The release includes:
- the MIT licence, a `CITATION.cff` file and a changelog;
- a pinned environment (`requirements-lock.txt`, Python 3.13.17);
- SHA-256 checksums and `$Id` revision lines for the bundled databases;
- an end-to-end example workflow for user data (`examples/user_data_workflow.py`);
- 35 automated tests (85% line coverage), covering mass and electron conservation for every scenario, saturation-index statuses, interpolation accuracy, the conformal constructions and a negative control for the equivalence test.

Every number and figure in this paper is regenerated by `scripts/run_v2.py`, with seeds listed in the script.

### 2.3 Chromium chemistry layer

**Databases.** PHREEQC 3 was used with the three chromium-bearing databases distributed by the USGS (Parkhurst and Appelo, 2013):
- `minteq.v4.dat` (MINTEQA2 v4 lineage; Allison et al., 1991);
- `minteq.dat` (MINTEQA2 v3 lineage);
- `llnl.dat` (LLNL thermo.com.V8.R6.230).

Their revision lines and checksums are given in Table S2. The databases bundled with phreeqpython contain no chromium species and are rejected programmatically.

**Cr(VI) speciation.** Cr was entered as `Cr(6)` with one K⁺ per Cr, representing a K₂Cr₂O₇ stock, in NaCl at the stated ionic strength. pH was imposed by charge balance on Cl⁻ or, where infeasible, on Na⁺. Valence states entered separately are not re-equilibrated in a PHREEQC initial solution, so this represents a freshly prepared Cr(VI) solution. Fractions are per Cr atom.

**Speciation grid.** To represent the concentration dependence of dimerisation wherever it matters, a grid was computed with `minteq.v4.dat`:
- pH 1.5–10 in 0.25 steps (35 values);
- dissolved Cr(VI) 0.01–700 mg L⁻¹ (11 values);
- I 10⁻³–10⁻¹ M (5 values);
- T 288–318 K (4 values).

That is 7,700 PHREEQC solutions, interpolated multilinearly. The simulator evaluates speciation at the *current* dissolved Cr(VI), pH, I and T. The ML features are evaluated at each row's C₀, pH, I and T. Interpolation error was checked against 300 direct PHREEQC calls at random interior points.

**Redox scenarios.** Total Cr was entered without a valence state at an assumed pe, restricted to the water-stability field (pe ≤ 20.75 − pH). The Eh at which Cr(VI) and Cr(III) are equally abundant was obtained by interpolation.

**Cr(III) solubility.** The v0.1 analysis relied on convergence of the PHREEQC phase-boundary option, which the reviewer rightly criticised. The saturation index of the Cr(III) hydroxide (Cr(OH)₃(am) in minteq.v4; Cr(OH)₃(A) in minteq) is now evaluated directly for solutions of prescribed dissolved Cr(III). Saturation is located by bisection on log₁₀ C between 10⁻⁶ and 5,000 mg L⁻¹. Every evaluation returns a status:
- `bracketed` (with the SI at the solution);
- `undersaturated_to_cmax`: SI < 0 even at 5,000 mg L⁻¹, a direct thermodynamic statement;
- `supersaturated_at_cmin`;
- `speciation_failed`: reported, never interpreted.

`llnl.dat` contains no Cr(III) hydroxide, and its only Cr(III) solid (eskolaite) is written through CrO₄²⁻ and O₂. A Cr(III)-only solubility is therefore undefined there and is reported as not computable.

**External validation.** An independent check against the critically evaluated constants of Ball and Nordstrom (1998) and Rai et al. (1987) requires their tabulated values. These could not be retrieved from the build environment. CASSIA accepts a reference-constant table for this comparison (Section 4.7).

### 2.4 Multi-mechanism known-truth simulator

**General model.** One kinetic model contains five removal pathways (`cassia.simulate`). Its states are:
- dissolved Cr(VI), C_v;
- adsorbed Cr(VI), q₆;
- Cr(III) formed on the surface, q₃s;
- Cr(III) re-adsorbed from solution, q₃r;
- dissolved Cr(III), C₃;
- precipitated Cr(III), P₃;
- surface reductant, R_s (mmol e⁻ g⁻¹);
- leached dissolved reductant, L (mmol e⁻ L⁻¹);
- pH.

**Rate expressions.**
- *Adsorption:* linear-driving-force approach to a Langmuir equilibrium, whose affinity depends on the HCrO₄⁻ + Cr₂O₇²⁻ fraction from the speciation grid and on surface charge relative to pH_pzc.
- *Direct surface-mediated reduction* of dissolved Cr(VI) (pathway B): r_B = k_B h(pH) (R_s/R₀) C_v.
- *Adsorption-coupled reduction* of adsorbed Cr(VI) in place (pathway C): r_C = k_C h(pH) (R_s/R₀) q₆.
- *Dissolved-phase reduction* by reductants leached from the carbon (pathway D): r_D = k_D h(pH) L C_v.
- In all three, h(pH) = [1 + 10^(pH − pH₅₀)]⁻¹, with an Arrhenius temperature dependence.

**Fate of Cr(III).**
- Surface-formed Cr(III) is retained with a pH-dependent affinity θ₃(pH) and released at a rate k_rel (1 − θ₃) q₃s.
- Dissolved Cr(III) re-adsorbs towards a Langmuir equilibrium.
- Cr(III) precipitates or redissolves towards the PHREEQC Cr(OH)₃(am) solubility at the current pH and I.

**pH drift (scenario E).** pH follows a proton balance: reduction consumes protons, precipitation releases three per Cr, and proton exchange with the carbon relaxes pH towards pH_pzc. The total buffer capacity is that of the carbon plus 2.303[H⁺].

**Conservation.** The model conserves both chromium mass and electron equivalents:
- C_v + D(q₆ + q₃s + q₃r) + C₃ + P₃ = C₀;
- 3 × (Cr reduced) = D(R₀ − R_s) + (L_leached − L).

Both are verified by unit tests for every scenario.

**Scenarios.**

| Scenario | Pathways switched on |
|---|---|
| A | adsorption only (negative control) |
| B | adsorption + direct surface reduction |
| C | adsorption + adsorption-coupled reduction |
| D | adsorption + dissolved-phase reduction |
| E | all pathways + Cr(III) re-adsorption + precipitation + pH drift |

**Material and study populations.**
- *Families and family-level parameters.* There are six families: chemically and physically activated carbons, wood and manure biochars, hydrochars, and biomass-derived porous carbons. Each has physically plausible ranges of S_BET and pH_pzc. The family-level means of reductant capacity, leachable reductant and specific capacity are redrawn for every corpus from a common hyper-prior, so which family is most reductive is itself random across replicates. This removes the circularity identified by the reviewer, where the most reductive family was fixed by construction.
- *Experiments per study.* Each virtual study contains an 8-point isotherm (C₀ = 5–300 mg L⁻¹ at an optimum pH of 1.8–3.5), a 6-point pH edge (pH 2–8), an 8-point kinetic series and a 6-point dose series (0.5–8 g L⁻¹ at C₀ = 20 mg L⁻¹).
- *Measurement model.* Cr(VI) and total Cr carry 3% relative error, detection limits and a study-level calibration bias. XPS Cr(III) fractions carry an absolute error of 0.05.
- *Evidence classes and reporting.* Classes are assigned at study level (C 0.40, B 0.35, A 0.25) and then resolved per observation. Ionic strength and final pH are reported with class-dependent probabilities.

Every parameter, its distribution and the numerical settings (LSODA, rtol 10⁻⁶, atol 10⁻⁸) are listed in Table S1. The simulator is one family of plausible mechanisms. It is not calibrated to any material; its purpose is to make the ground truth available.

### 2.5 Mechanism benchmark (conventional isotherm analysis)

**Design.** For each scenario, 100 replicate corpora of 30 studies were generated (15,000 virtual studies in total). In every study, the isotherm series was analysed as in conventional practice:
- nonlinear Langmuir, Freundlich, Sips and Temkin fits (Temkin at the experimental T), with AICc/BIC and Akaike weights (Burnham and Anderson, 2002);
- Shapiro–Wilk and runs tests on the residuals;
- PFO/PSO kinetic fits.

In 20 replicates per scenario, 200-replicate residual-bootstrap intervals were computed for Langmuir q_max.

**Outcomes.** The outcomes are:
- the true shares of Cr(VI) disappearance attributable to adsorbed Cr(VI), Cr(III) on the solid, dissolved Cr(III) and precipitate (mass-weighted over the isotherm series);
- the fraction of fits with R² ≥ 0.95;
- the ratio of fitted to true q_max;
- bootstrap-interval coverage of the true q_max;
- for Class A/B studies, the aqueous Cr(III) index.

Across-replicate variability gives Monte Carlo standard errors (MCSE).

### 2.6 Mechanism-identifiability engine

Four cumulative measurement designs were defined:

| Design | Measurements | Evidence class |
|---|---|---|
| D1 | residual Cr(VI) at equilibrium (isotherm + pH edge) | C |
| D2 | D1 + Cr(VI) kinetics | C |
| D3 | D2 + total dissolved Cr at every point | B |
| D4 | D3 + XPS Cr(III) fraction on three solids + final pH at every point | A |

Four complementary analyses were applied.

1. **Local identifiability.** For each study, the sensitivity matrix J of all observables was computed with respect to the 14 log-transformed or pH-scale kinetic parameters plus the study bias, by forward differences. The Fisher information JᵀWJ (W = inverse measurement variance) was combined with a weakly informative Gaussian prior (SD 0.5 log₁₀ units). The Laplace-approximation posterior SD of the reduction share then follows by the delta method and was compared with its prior SD. This analysis assumes the mechanistic model structure is known, so it is an optimistic bound.
2. **Model-free detection tests.** Three diagnostics were each assessed over 20 noisy replicates per study: lack of isotherm fit (best R² < 0.95 or runs p < 0.05; D1), a one-sided z-test that dissolved Cr(III) = Σ(C_e,tot − C_e,VI) exceeds analytical noise (D3), and XPS Cr(III) fraction > 0.15 (D4).
3. **Structural equivalence test.**
   - *Question:* can an adsorption-only model reproduce the observations of each design? This is the reviewer's suggested test of whether different mechanism combinations give indistinguishable data.
   - *Model:* the flexible empirical model a conventional analyst would fit. It has a pH-dependent Sips equilibrium with a free logistic pH response, an Avrami-type kinetic curve with its own q_e, a free study bias and free pH relaxation. Dissolved and solid Cr(III) are zero by construction.
   - *Fit and test:* fitted by multi-start weighted least squares, then a χ² goodness-of-fit test, with p > 0.05 taken as 'indistinguishable'. Each study was paired with a negative control: the same material with reduction switched off. The test was run at 3% and 10% analytical noise.
4. **Profile likelihood and metadata bias.**
   - *Profile likelihood* (Raue et al., 2009): the reductant capacity R₀ was profiled, with all other parameters re-optimised, for one low-reduction and one high-reduction study under D1, D3 and D4. The reduction share implied along each profile was recorded.
   - *Metadata bias:* the first-order bias in the inferred reduction share was propagated for misreported pH (+0.5) and ionic strength (×10), via Δθ = Σ_post JᵀW(y_true − y_assumed).

Sample sizes: 300 studies (200 from scenario E and 100 from B–D) for analyses 1, 2 and 4; 175 studies, with paired controls, for analysis 3.

### 2.7 Machine-learning models and validation

**Targets and transforms.**
- *Targets:* q_app,VI (all rows), q_tot,aq (rows with total Cr) and, for the design layer, residual concentrations C_e,VI and C_e,tot.
- *Transform:* each is modelled as log₁₀(y + c), with c = 0.1 mg g⁻¹ for uptake and LOD/2 for concentrations.
- *Back-transformation:* uses Duan's (1983) smearing factor estimated on training residuals. The naive inverse is reported for comparison.

**Features.**
- *Base set:* C₀, pH, D, t, T, log₁₀ I, an I-imputed flag, S_BET, V_p, d_p and pH_pzc. Unreported I is imputed at 0.01 M and flagged.
- *Chem set:* adds the PHREEQC speciation features and pH − pH_pzc.

**Models.** Hyperparameters were fixed a priori (Table S3): RF (300 trees, min. leaf 2), XGBoost (500 rounds, depth 4, learning rate 0.05) and GPR (Matérn-5/2 ARD + white noise, standardised inputs). All seeds were fixed. Nested tuning by inner 4-fold grouped CV (RF: 9 configurations; XGBoost: 4) was evaluated on three corpora to quantify its effect.

**Validation design.**
- *Corpora:* 20 independent scenario-E corpora of 30 studies (six families; 840 observations each).
- *Schemes:* random 5-fold (diagnostic only), grouped 5-fold by study, leave-one-study-out (LOSO) and leave-one-family-out (LOFO, six held-out families).
- *External test:* models trained on each corpus were tested on an *external* corpus generated with different hyper-prior draws, so that population shift is included.
- *GPR:* evaluated on three corpora because of its cost. Its ARD kernel hyperparameters were optimised on a seeded 250-row subsample of each training fold and then held fixed while the GP was refitted on all training rows, so no test information is used.

**Reporting.** Each corpus reports pooled R², RMSE, MAE and median absolute error, with 1,000-replicate study-level (cluster) bootstrap intervals and per-study errors. Across corpora the distribution of each metric is reported, and model differences are paired by corpus. Between-study variability is quantified by a random-intercept decomposition of out-of-study residuals (τ², σ², ICC; Gelman and Hill, 2007). A learning curve (12–96 studies; four replicates each) estimates how many independent studies a transferable model needs.

### 2.8 Study-level conformal prediction

Three interval constructions were compared. All are calibrated on out-of-study residuals, with the study as the unit.

- **M0 – pooled residuals (v0.1 method).** The half-width is the ⌈(n+1)(1−α)⌉/n quantile of pooled absolute residuals from an inner grouped 10-fold CV. Residuals within a study are dependent, so M0 has no formal guarantee for a new study.
- **M1 – study-level CV+ with study-weighted scores.**
  - *Construction:* the CV+ interval of Barber et al. (2021), formed with study-level folds (K = 10). Each calibration study contributes total weight 1 (weight 1/n_j per observation), as in the weighted pooling of Dunn et al. (2023), and the test study carries one further unit of weight.
  - *Target:* marginal coverage of a random observation from a new study.
- **M2 – study-level CV+ with study scores.**
  - *Score:* s_j is the within-study 0.9-quantile of absolute residuals.
  - *Interval:* bounded by the ⌊α(J+1)⌋-th and ⌈(1−α)(J+1)⌉-th order statistics of the fold-model predictions ∓ s_j, over the J calibration studies.
  - *Target:* that at least 90% of a new study's observations are covered, with probability ≥ 1 − 2α under exchangeability of studies (the CV+ bound).
  - *Requirement:* M2 needs J ≥ (1 − α)/α studies; with fewer, it correctly returns an unbounded interval.

**Evaluation.** Each construction was evaluated by outer LOSO on 10 corpora (29 calibration studies per test study) for all four targets, at nominal levels 0.50–0.95. Reported for each:
- mean, median, 10th-percentile and worst per-study coverage;
- the fraction of studies covered at the nominal level;
- the observation-weighted coverage;
- the mean width;
- the number of calibration studies.

### 2.9 Treatment-design layer

**Models.** Residual Cr(VI) and residual total Cr are modelled directly (Section 2.7), each with its own M1 interval. Modelling effluent through q, as in v0.1, is ill-conditioned when removal approaches completeness.

**Design problem.**
- *Grid:* pH 2–7 in 0.5 steps and dose 0.5–8 g L⁻¹ (the training range), 24 h contact.
- *Rule:* choose the lowest dose meeting both limits.
- *Targets:* S1, C₀ = 20 mg L⁻¹ with limits Cr(VI) ≤ 0.5 and total Cr ≤ 2 mg L⁻¹; S2, C₀ = 5 mg L⁻¹ with limits 0.1 and 0.5 mg L⁻¹.

**Decision rules.**
- *Point:* point predictions.
- *Marginal 90%:* each endpoint's upper bound at α = 0.10.
- *Joint 90%:* Bonferroni, each endpoint at α = 0.05, so that the joint miscoverage is ≤ 0.10 by the union bound.
- *Joint 90% / 80% + 4 confirmatory batches:* four batches on the new material (pH 3, doses 0.5–4 g L⁻¹, 3% analytical noise) estimate its study offset for each endpoint. The centre is shifted by that offset. The half-width is the study-weighted (1 − α) quantile of offset-adjusted out-of-study residuals of the training studies, with α = 0.05 or 0.10 per endpoint (heuristic; it assumes the confirmatory batches are exchangeable with the study's other observations).
- *Point + verify (model-guided, test-confirmed):* the point model proposes the pH and the lowest dose. Duplicate confirmatory batches are run at that design, and the design is accepted when both mean measured concentrations are below 90% of their limits. Otherwise the dose is escalated along the grid, for at most four rounds.

**Cohorts.** Twenty independent cohorts were used, each with its own 30-study training corpus and 15 new materials drawn from the same hyper-priors (300 materials per target).

**Truth and metrics.** For every material, the simulator gave the true attainability: the lowest grid dose at which any pH meets both limits. For every rule and target the following were reported:
- unsafe recommendations (recommended but non-compliant in truth), with Wilson 95% intervals (Wilson, 1927);
- compliant designs found among attainable materials;
- infeasibility declarations for attainable materials;
- reliability among recommendations;
- the dose overhead relative to the true minimum dose.

**Sensitivity.** The operating-window sensitivity of point + verify designs was evaluated at pH ± 0.25 and ± 0.5 and dose × 0.8 and × 1.2. Predicted and true effluent were compared on a sub-grid.

### 2.10 Propagation of unreported chemistry

**Prediction level.** For studies that did not report ionic strength or final pH, 1,000 scenarios were drawn: I log-uniform on 10⁻³–10⁻¹ M and a pH-basis ambiguity uniform on ±0.5. They were propagated through the speciation features, the q_app,VI predictions and the M1 intervals. For each observation the following were recorded: the 90% range of f(HCrO₄⁻), the spread of the point prediction, and the union interval width against the nominal width. Grouped 5-fold models were used on four corpora.

**Design level.** For six new materials per corpus, the designer's assumed I was drawn 300 times from the same range, while the true I was 0.01 M. The point + verify design was re-run for each draw (confirmatory batches at the true conditions), and the true outcome of each accepted design was simulated. This gives distributions of recommended dose, residual total Cr and compliance probability.

**Mechanistic level.** The first-order bias in the inferred reduction share is given by the identifiability engine (Section 2.6).

### 2.11 Literature-data protocol and seed dataset

The extraction protocol (`data/PROTOCOL.md`) is PRISMA 2020-aligned (Page et al., 2021). It specifies:
- the databases: Scopus, Web of Science, PubMed, Zenodo, Mendeley Data, Figshare and Data in Brief;
- the exact search strings;
- the eligibility rules;
- table-first extraction with digitisation as a recorded fallback;
- unit normalisation and provenance fields;
- quality-appraisal fields (replicates, SD, blanks, mass balance, pH drift);
- observation-level evidence classification.

The seed dataset is Ramirez et al. (2020), Data in Brief 29, 105292: activated carbon from teak-wood residue, with equilibrium and kinetic data in the supplementary files. It is registered in `data/study_registry.csv`. Execution of the search and the extraction is pending (Section 3.8).

### 2.12 Use of AI tools in the research process

The CASSIA code, the analysis scripts and the manuscript drafts were developed with an AI coding assistant (Claude, Anthropic), working under the author's direction in a version-controlled repository.

**Safeguards.** The scientific design and the reviewer-driven changes were specified by the author. The correctness checks built into the work are:
- unit tests for conservation laws and numerical constructions;
- negative controls (scenario A);
- cross-checks of interpolated against direct PHREEQC calculations;
- full regeneration of every reported number from seeded scripts.

[The author should add a statement of their own verification of code and results here before submission.]

## 3. Results

Sections 3.1 and 3.7 report thermodynamic calculations. Sections 3.2–3.6 report **synthetic, known-truth benchmarks**. They are statements about what the methods can and cannot recover under the simulated mechanisms, not empirical findings about real carbons. Section 3.8 reports the status of the empirical application.

### 3.1 Chromium chemistry across three databases

**Cr(VI) speciation is robust to the database choice (Table 2, Fig. 1a,b).**
- The log K values for HCrO₄⁻ and Cr₂O₇²⁻ agree to within 0.016 and 0.041 log units across the three databases.
- Over pH 1.5–10 and 10–500 mg L⁻¹, species fractions differ from minteq.v4 by at most 0.031 (minteq; HCrO₄⁻ at pH 1.5, 500 mg L⁻¹) and 0.044 (llnl; CrO₄²⁻ at pH 10, 500 mg L⁻¹). The llnl difference reflects the alkali–chromate ion pairs (NaCrO₄⁻, KCrO₄⁻) that only the MINTEQ databases contain.
- The HCrO₄⁻/CrO₄²⁻ crossover lies at pH 6.37–6.38 in all three databases (10 mg L⁻¹, I = 0.01 M).
- The crossover shifts with ionic strength, from pH 6.46 (I = 0.001 M) to 6.24 (I = 0.1 M). The temperature effect is small: pH 6.37 at 288 K, 6.40 at 318 K.
- H₂CrO₄ is the only species with a large log K disagreement (6.42, 5.65 and 5.18). It never exceeds about 2% of Cr(VI) above pH 1.5.

**Concentration matters more than the database.**
- At pH 3 the dichromate share of Cr rises from 0.17% at 1 mg L⁻¹ to 1.7% (10), 7.4% (50), 13.1% (100), 28.2% (300) and 43.7% (700 mg L⁻¹) in minteq.v4.
- The three databases agree to within 0.02 at every concentration.
- A single isotherm spanning 5–300 mg L⁻¹ therefore samples materially different sorbate mixtures. This is why the revised simulator evaluates speciation at the current dissolved Cr(VI) (reviewer comment 3A).
- The four-dimensional speciation grid reproduced 300 direct PHREEQC solutions with a maximum absolute error of 0.011 and a mean error of 0.0011 in species fraction.

**Table 2.** Equilibrium constants (25 °C) for CrO₄²⁻ + n H⁺ ⇌ product.

| Product | minteq.v4 | minteq | llnl |
|---|---|---|---|
| HCrO₄⁻ | 6.510 | 6.509 | 6.494 |
| Cr₂O₇²⁻ (2 CrO₄²⁻ + 2 H⁺) | 14.560 | 14.557 | 14.519 |
| H₂CrO₄ | 6.419 | 5.651 | 5.175 |

**The redox boundary is database dependent (Fig. 1c).**
- The Eh at which Cr(VI) and Cr(III) are equally abundant (50 mg L⁻¹ total Cr) is 63–83 mV higher in minteq.v4 than in llnl over pH 2.5–10. The two MINTEQ versions agree to within −19 to +5 mV.
- At pH 2, neither MINTEQ database reaches 50% Cr(VI) anywhere inside the water-stability field. At the O₂/H₂O limit, only 12% Cr(VI) remains in minteq.v4.
- Under the conditions where Cr(VI) removal by carbons is usually optimal, Cr(VI) is thermodynamically unstable with respect to reduction at almost any attainable Eh. Its persistence is kinetic.
- The 63–83 mV spread is a *between-database* difference, not an independently validated uncertainty interval for real systems (reviewer comment 3C).

**Direct saturation indices replace phase-boundary convergence (Fig. 1d, Table S4).**
- Of 174 bisection evaluations (two databases × three ionic strengths × 29 pH values), 146 were bracketed and 28 returned 'undersaturated up to 5,000 mg L⁻¹'. None failed.
- A 625-point SI map over pH 3–9 and 10⁻³–10³ mg L⁻¹ Cr(III) converged everywhere.
- Undersaturation up to 5,000 mg L⁻¹ occurs at pH ≤ 4.0 in minteq.v4 and pH ≤ 3.5 in minteq. This is now a direct thermodynamic result, not an inference from non-convergence (reviewer comment 3B).
- At I = 0.01 M, the dissolved Cr(III) at saturation with Cr(OH)₃(am) is:

| pH | minteq.v4 (mg L⁻¹) | minteq (mg L⁻¹) |
|---|---|---|
| 4.0 | undersaturated up to 5,000 | 266 |
| 4.5 | 16 | 9.5 |
| 5.0 | 1.42 | 0.79 |
| 6.0 | 0.023 | 0.017 |
| 7.0 | 1.2 × 10⁻³ | 1.8 × 10⁻³ |
| 8.0 | 1.4 × 10⁻⁴ | 7.8 × 10⁻⁴ |

- The two MINTEQ versions therefore differ by more than an order of magnitude at pH 4, and by a factor of 5.5 at pH 8.
- llnl.dat cannot provide this quantity (Section 2.3).

![Fig. 1](figures/fig1_chemistry.png)

**Fig. 1.** PHREEQC chromium chemistry with three databases (thermodynamic calculations; no experimental data). (a, b) Cr(VI) species fractions per Cr atom at 10 and 500 mg L⁻¹, I = 0.01 M NaCl. (c) Eh at which Cr(VI) and Cr(III) are equally abundant (50 mg L⁻¹ total Cr), with the O₂/H₂O boundary. (d) Dissolved Cr(III) at saturation with Cr(OH)₃(am) by direct saturation-index bisection. Triangles mark pH values at which solutions remain undersaturated up to 5,000 mg L⁻¹ (filled: minteq.v4; open: minteq).

### 3.2 Mechanism benchmark: what conventional isotherm analysis recovers (synthetic)

**Scale and mechanism mix.** Across 15,000 virtual studies (5 scenarios × 100 corpora × 30 studies), the true share of Cr(VI) disappearance due to reduction (isotherm series, mass-weighted) had these medians:
- scenario B: 0.24 (5–95%: 0.02–0.92);
- scenario C: 0.22 (0.03–0.75);
- scenario D: 0.04 (0.002–0.54);
- scenario E: 0.38 (0.05–0.93).

In scenario E the removed chromium was distributed as follows (medians): 0.62 adsorbed Cr(VI), 0.24 dissolved Cr(III), 0.045 Cr(III) on the solid and essentially no precipitate (95th percentile 0.065), consistent with the Cr(OH)₃ undersaturation at the isotherm pH of 1.8–3.5 (Section 3.1). The distributions are shown in Fig. 2a,d.

**Goodness of fit does not discriminate mechanisms (Fig. 2b).**
- *Negative control.* In the pure-adsorption control (A), only 46.7% of best-isotherm fits reached R² ≥ 0.95 (MCSE 0.8%). For low-capacity families (manure biochar, hydrochar; median true q_max 1–4 mg g⁻¹), only 16% did, because their adsorption signal is close to analytical noise.
- *Reduction can improve fit quality.* The share of R² ≥ 0.95 fits rose to 58.8% (B) and 59.3% (E). For the low-capacity families it roughly doubled, to 31–34%, because reduction adds a smooth, saturating removal signal.
- *High-reduction studies still fit well.* Among studies whose reduction share exceeded 0.2, R² ≥ 0.95 was reached in 46% (B), 41% (C), 17% (D) and 52% (E).
- *Residual diagnostics did not help.* The runs test flagged 4.9–7.9% of Langmuir fits in every scenario, including the control.
- *Kinetic model choice did not help.* PSO was preferred over PFO in 25% of control series and 27–34% of reduction scenarios, so it carries no mechanistic information.

**Publishable fits overstate capacity in proportion to reduction, and the size depends on the mechanism (Fig. 2c, Table 3).** Fits whose q_max reached the fitting bound because the isotherm never saturated (0.6–3.1%) were excluded.

**Table 3.** Fitted/true Langmuir q_max for publishable fits (R² ≥ 0.95), by scenario (synthetic; 3,000 studies per scenario).

| Scenario | Median ratio (IQR) | Median ratio, reduction share > 0.3 (n) | Spearman ρ, reduction share vs ratio (mean ± SD over replicates) |
|---|---|---|---|
| A (control) | 1.02 (0.97–1.18) | – | – |
| B | 1.32 (1.07–2.60) | 4.88 (604) | 0.60 ± 0.19 |
| C | 1.15 (1.04–1.50) | 1.60 (446) | 0.32 ± 0.19 |
| D | 1.04 (0.98–1.21) | 13.7 (43; mostly very-low-capacity materials) | 0.13 ± 0.20 |
| E | 1.31 (1.08–2.29) | 2.25 (857) | 0.63 ± 0.17 |

- With reduction shares below 0.05, the median ratio was 1.02–1.03 in every scenario.
- The qualitative conclusion is that a saturating fit to apparent removal absorbs the reduction sink into q_max. It holds for direct surface reduction (B), adsorption-coupled reduction (C) and the combined mechanism (E).
- It is weak for dissolved-phase reduction (D), which rarely dominates removal at the simulated leachable-reductant levels.

**Bootstrap intervals are miscalibrated for reasons beyond mechanism.**
- The 95% residual-bootstrap interval of Langmuir q_max covered the true q_max in only 22% of publishable control fits (n = 289).
- Coverage fell further to 8.8–17.4% in the reduction scenarios.
- Eight-point isotherms that do not fully reach equilibrium produce intervals that are both precise and biased, even without reduction. Reduction makes this worse.

**A total-Cr mass balance tracks reduction where Cr(III) stays dissolved.** For Class A/B studies, the aqueous Cr(III) index ranked studies by reduction share with Spearman ρ = 0.73 (B), 0.52 (C), 0.33 (D) and 0.71 (E). It underestimated the share by a median of 23% in C and E, where Cr(III) is partly retained on the solid. It is therefore a lower bound, as intended.

**Family patterns follow from assumed capacity, not from assumed reductance.**
- Even with family-level reductant priors redrawn in every replicate, hydrochar had the highest mean reduction share in 83 of 100 scenario-E replicates.
- With hydrochar excluded, the most reductive family varied (manure biochar 56, wood biochar 39, physically activated carbon 5).
- Within scenario E the drawn family reductant prior correlated with the reduction share (ρ = 0.40), as intended. The dominance of low-surface-area chars, however, comes from the assumed proportionality between S_BET and adsorption capacity: a small adsorption sink inflates the reduction *share*.
- This is a property of the simulator and must not be read as an empirical statement about hydrochars (reviewer comment 2).

![Fig. 2](figures/fig2_benchmark.png)

**Fig. 2.** Mechanism benchmark (synthetic; 15,000 virtual studies). (a) True reduction share of Cr(VI) disappearance by scenario (box: IQR; whiskers: 5–95%). (b) Share of studies whose best isotherm fit reaches R² ≥ 0.95, against the true reduction share (bins with ≥ 20 studies); dotted line, pure-adsorption control. (c) Fitted/true Langmuir q_max for fits with R² ≥ 0.95 (median and IQR; log scale). (d) Fate of the removed Cr(VI) in scenario E. Scenario definitions are in Section 2.4 and all parameter distributions in Table S1.

### 3.3 Identifiability: which measurements establish the mechanism (synthetic)

**If the model were known (Fig. 3a).**
- With the mechanistic model structure assumed known, residual Cr(VI) at equilibrium (D1) already constrained the reduction share. The median Laplace posterior SD was 0.014, 0.024, 0.032 and 0.065 for true shares of < 0.05, 0.05–0.2, 0.2–0.4 and > 0.4, against prior SDs of 0.04–0.28.
- Adding total Cr (D3) roughly halved these values, to 0.008–0.035, and adding XPS and final pH (D4) reduced them to 0.007–0.021.
- This analysis is deliberately optimistic. It assumes that the analyst knows the exact rate laws, which no analyst of real data does.

**Model-free diagnostics (Fig. 3b).**
- *Lack of isotherm fit* is not specific. It flagged 35% of studies with < 5% reduction and 58% of studies with > 40%.
- *The dissolved-Cr(III) z-test* had a 8% detection rate at < 5% reduction, rising to 40%, 64% and 66% in the higher bins.
- *XPS* had rates of 12%, 51%, 59% and 75%.
- The two Cr(III)-based diagnostics are complementary. The z-test misses Cr(III) that stays on the solid; XPS misses Cr(III) released to solution.

**Structural equivalence: an adsorption-only model explains most Class C data (Fig. 3c, Table 4).** This is the decisive result.

**Table 4.** Structural equivalence test: share of studies in which a flexible adsorption-only model is not rejected (χ² p > 0.05).

| Design | Studies with reduction share > 0.3 where the adsorption-only model is *not* rejected, 10% noise | Same, 3% noise | False-rejection rate of paired pure-adsorption controls |
|---|---|---|---|
| D1 | 79% | 67% | 3.4% |
| D2 | 72% | 65% | 10–13% |
| D3 | 32% | 7% | 10–13% |
| D4 | 12.5% | 3% | 10–13% |

- The flexible adsorption-only model was fitted to data generated with the full mechanism (n = 175 studies plus 175 paired pure-adsorption controls).
- The control false-rejection rates for D2–D4 are slightly above nominal. This reflects imperfect flexibility of the empirical kinetic model.
- Residual Cr(VI) data therefore usually cannot exclude the hypothesis 'adsorption only', even when most of the Cr(VI) was reduced. The Cr(III) observations of Classes B and A usually can.
- Kinetic Cr(VI) data (D2) add little, because conventional kinetic models are fitted with their own free equilibrium uptake.

**Profile likelihood (Fig. 3d).** All profile points converged; the minimum −2 log L was consistent with the number of observations (3.4 for n = 14 under D1; 21–79 for n = 44–61 under D3/D4).

- *Low-reduction study (true share 0.15).* The reductant capacity R₀ was practically non-identifiable. The 95% profile region extended from at least 16-fold below to 1.6-fold above the true value under every design. The implied reduction share nevertheless stayed within 0.13–0.14 (D1) and 0.17–0.18 (D3, D4). A share can be well determined while the parameter that generates it is not, *provided the model structure is correct*.
- *High-reduction study (true share 0.60).* R₀ was sharply identified, within one grid step of 0.2 log units, under all designs. The share at the minimum was 0.601 in each case.
- *Caveat.* The profiles use the known model structure and include the true parameters among the optimiser starts. They therefore show identifiability under the correct model, not the ease of estimation. The equivalence test removes that assumption, and under it Class C data lose most of their discriminating power.

**Misreported metadata.**
- Analysing the data at a pH 0.5 units above the true value biased the inferred reduction share by a median of 0.032 (90th percentile 0.157) under D1, and by 0.024 (0.092) under D3.
- A ten-fold error in ionic strength caused much smaller biases: medians of 0.003 and 0.001, 90th percentiles of 0.025 and 0.015.
- The pH basis (initial vs final) is therefore the reporting omission that matters most for mechanism inference.

![Fig. 3](figures/fig3_identifiability.png)

**Fig. 3.** Mechanism identifiability by measurement design (synthetic). (a) Laplace-posterior SD of the reduction share assuming the true model structure (median per bin; dotted: prior). (b) Detection rates of three model-free diagnostics. (c) Fraction of studies in which a flexible adsorption-only model is not rejected (χ² p > 0.05), for paired pure-adsorption controls and for studies with moderate or high reduction, at 3% and 10% analytical noise. (d) Profile likelihood of the reductant capacity for one low-reduction (dashed; non-monotone, with a secondary basin above 6× R₀) and one high-reduction (solid) study under D1, D3 and D4. The legend gives the true reduction share and the share at the profile minimum; dotted line, 95% threshold (3.84).

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

## 4. Discussion

### 4.1 What can be established from which measurements

The revised analysis gives a mechanism-resolved answer to the question in the title. It rests on synthetic evidence and should be read with that qualification.

**Class C (residual Cr(VI) only).**
- Such data establish that Cr(VI) disappears under the stated conditions, and how fast.
- They do not establish why. A flexible adsorption-only model fits Class C data from studies in which most of the Cr(VI) was reduced (79% not rejected at 10% analytical noise).
- Conventional isotherm statistics do not repair this. Fit quality is weakly, and sometimes positively, related to reduction; residual diagnostics flag 5–8% of fits regardless of mechanism; bootstrap intervals for q_max are precise but rarely cover the true capacity.
- The fitted q_max is an *apparent Cr(VI) removal capacity*. Under surface or adsorption-coupled reduction it exceeds the adsorption capacity in proportion to the reduction share.

**Class B (adds a total-Cr analysis).**
- Adding total Cr changes the picture qualitatively. The adsorption-only hypothesis is then rejected in most high-reduction studies (68–93%, depending on noise).
- The aqueous Cr(III) index ranks studies by reduction share, as a lower bound.

**Class A (adds solid-phase speciation and final pH).**
- This closes most of the remaining gap: the adsorption-only model is rejected in 88–97% of high-reduction studies.
- XPS and the dissolved-Cr(III) test are complementary. One sees Cr(III) retained on the solid, the other Cr(III) released to solution.

**Recommendation.** This hierarchy is a quantitative argument for the minimum reporting set of §4.4. It also gives a criterion for pooling literature data: Class C observations may train a model of apparent removal, but they must not be pooled with Class A/B observations as evidence of adsorption capacity.

**Why the known-structure analyses look better.** The Fisher and profile analyses are much more optimistic than the equivalence test, and the contrast is instructive. If the analyst knew the exact rate laws, even Class C data would constrain the reduction share. Real analysts do not, and once the structure is unknown, the information about mechanism in Class C data largely disappears. Identifiability here is a property of the data *and* of what the analyst is entitled to assume.

### 4.2 Thermodynamic feasibility versus surface kinetics

**What thermodynamics says.** At the pH where carbons remove Cr(VI) most effectively (1.5–3.5), Cr(VI) is unstable with respect to Cr(III) at almost any Eh inside the water-stability field (§3.1). The two MINTEQ databases do not even reach a 50% Cr(VI) fraction at pH 2. Cr(III) is also undersaturated with respect to Cr(OH)₃(am) up to thousands of mg L⁻¹ below pH ≈ 4.

**What it cannot say.** Equilibrium thermodynamics therefore says that reduction is possible and that its product will stay dissolved at low pH. It says nothing about whether reduction happens on the time scale of a batch test. That depends on:
- the rate of electron transfer from the carbon;
- its finite reductant capacity;
- the competition between adsorption and reduction.

These are kinetic and material-specific properties that equilibrium isotherms cannot represent, and that the simulator represents only through assumed rate laws.

**Consequence.** The redox calculations in CASSIA are labelled *scenarios*. A PHREEQC equilibrium prediction is never presented as a prediction of reduction in a batch experiment.

**Caveat on Eh.** The 63–83 mV difference between databases shows that even the thermodynamic boundary is uncertain. Without measured Eh, the redox state of a published experiment cannot be reconstructed.

### 4.3 Machine learning across studies

The random-versus-study gap (R² 0.89 against 0.63) reproduces, in this domain, what is known from ecology and materials informatics (Roberts et al., 2017). Three aspects matter for practice.

1. **One corpus cannot rank learners.** LOSO R² varied widely between corpora of identical design (95% range 0.41–0.82), and within-corpus cluster-bootstrap intervals were about ±0.15. A single literature corpus of 20–30 studies therefore cannot establish small differences between learners; the apparent RF–XGBoost advantage (0.04) was not significant.
2. **The error is a study-level offset.** Nearly half of the out-of-study error is a study-level offset (ICC 0.45). Neither more rows from the same studies nor hyperparameter tuning removes it; the learning curve shows that independent studies are the limiting resource. The offset also explains why study-level intervals must be wide (§3.5), and why confirmatory measurements on the new material are so valuable (§3.6).
3. **Chemistry features do not raise accuracy.** Chemistry-informed features did not improve accuracy even though, in the simulator, adsorption affinity depends on speciation by construction. Tree ensembles reconstruct smooth functions of pH, C₀, I and T from the raw inputs. In CASSIA, the role of these features is to make chemistry explicit and its uncertainty calculable, not to raise R².

### 4.4 A minimum reporting standard for mechanistic claims

The results support five requirements for any paper that claims Cr(VI) *adsorption* by a carbonaceous material:
1. total Cr alongside Cr(VI) for at least the isotherm and pH-edge series;
2. solid-phase Cr speciation (XPS or XANES) for at least the lowest- and highest-loading solids;
3. final as well as initial pH;
4. the background electrolyte and ionic strength;
5. replicate analytical SDs.

**Why these.** Items 1–2 move a study from Class C to Class A, which is the step that makes mechanisms distinguishable (§3.3). Item 3 matters most among the metadata: a 0.5-unit pH-basis ambiguity biases the inferred reduction share by up to 0.16 (90th percentile), whereas a ten-fold error in ionic strength biases it by at most 0.03.

**Using existing literature.** Studies that do not meet the standard remain useful, provided their endpoint is labelled as apparent removal. CASSIA enforces this label automatically.

### 4.5 Positioning

**What CASSIA does not offer.** It is not a new chromium-uptake predictor, and the paper does not report a better R².

**What it adds.**
- **Provenance and evidence at observation level.** Every number is tied to its source and to the measurements that support its interpretation.
- **A mechanism-identifiability layer.** It tells the user which conclusions a given set of measurements can support, and which additional measurement would change that.
- **Study-level uncertainty and decision support.** Prediction intervals and treatment designs are calibrated at the level that matters for a new material, the study.

Existing adsorption databases and ML studies of Cr(VI) removal typically pool residual-Cr(VI) observations without these distinctions. A systematic comparison with them belongs with the empirical application (§3.8).

### 4.6 Scope decisions: pyGAPS and RASPA3

**pyGAPS.** It recomputes BET areas, with Rouquerol criteria, from published N₂ isotherms (a unit-tested adapter), so that textural descriptors are consistent across studies. It is never applied to aqueous data.

**RASPA3.** It was audited and not used. It offers mature Monte Carlo and molecular-dynamics capabilities, but no model of aqueous Cr(VI)–carbon surface reactions, and its reaction-ensemble functions are experimental. A narrowly defined optional module, for example the pore volume accessible to hydrated HCrO₄⁻ versus Cr₂O₇²⁻ in idealised slit pores, could inform accessibility descriptors. It could not inform reduction, which this paper shows to be the decisive and unidentifiable process.

### 4.7 Limitations

1. **No empirical validation.** All findings in §3.2–3.7 come from simulations. The systematic search and extraction, including the seed dataset of Ramirez et al. (2020), have not been carried out because the publisher and repository hosts were not reachable from the computational environment. This is the main limitation of the paper, and the next step of the work.
2. **One model family.** The simulator spans five mechanisms and randomised family-level priors, but all are variants of one kinetic model family. Results that depend on assumed scaling, such as the dominance of low-surface-area chars in reduction share, are properties of the simulator. Mechanisms not represented include electron shuttling by dissolved organic matter, Cr(V) intermediates, competitive anion effects and surface-site heterogeneity.
3. **Thermodynamic validation.** Only three USGS databases were compared. An independent validation against critically evaluated constants (Ball and Nordstrom, 1998; Rai et al., 1987) is outstanding. The background electrolyte was NaCl, and major-ion effects beyond it were not explored.
4. **Simplified measurement model.** XPS, for example, is treated as an unbiased bulk measurement of the Cr(III) fraction. Real XPS probes only a surface layer and can itself reduce Cr(VI) under the beam.
5. **Conformal guarantees.** The conformal guarantees assume exchangeable studies. A new material from an unrepresented family violates this assumption, as the leave-one-family-out results show.
6. **Approximate procedures.** The identifiability analyses use local (Fisher) approximations and profiles computed with the true parameters among the optimiser starts. The equivalence test uses one flexible empirical model class; a more flexible class would make Class C data even less informative.

## 5. Conclusions

### 5.1 Demonstrated computational findings (thermodynamic calculations and known-truth simulations)

1. Three PHREEQC databases agree on Cr(VI) speciation to within 0.044 in species fraction. They differ by 63–83 mV on the Cr(VI)/Cr(III) boundary and by more than an order of magnitude on Cr(OH)₃(am) solubility at pH 4. Direct saturation indices show Cr(III) undersaturated up to 5,000 mg L⁻¹ at pH ≤ 3.5–4.0.
2. Across 15,000 virtual studies under five mechanisms, isotherm fits of residual Cr(VI) with R² ≥ 0.95 occurred in 41–52% of studies with more than 20% reduction under surface, adsorption-coupled and combined reduction. Fitted q_max then exceeded the adsorption capacity in proportion to the reduction share. The effect was weak when reduction occurred only in solution.
3. A flexible adsorption-only model could not be rejected for 79% of high-reduction studies observed through residual Cr(VI) alone (10% noise). Adding total Cr reduced this to 32%, and adding XPS and final pH to 12.5%.
4. Random cross-validation overstated transferability: random-forest R² was 0.89 under random splits, 0.63 under leave-one-study-out, 0.56 under leave-one-family-out and 0.56 on external corpora. Speciation features and nested tuning did not improve prediction. Forty-five per cent of out-of-study variance was between-study, and leave-one-study-out R² rose with the number of independent studies (0.47 at 12, 0.75 at 96).
5. Pooled-residual intervals reached nominal *average* coverage but covered only 65–75% of new studies at 90%. Study-level CV+ intervals with one score per study covered 88–93% of studies, at widths of 2.1–3.1 log₁₀ units. 95% intervals for total-Cr endpoints require at least 19 calibration studies with total-Cr data.
6. Point predictions gave unsafe treatment designs for 21% (S1) and 7% (S2) of new materials, while calibrated bounds certified almost nothing. A model-guided, test-confirmed protocol found compliant designs for 54–62% of attainable materials with no unsafe recommendation, at a median dose twice the true minimum.

### 5.2 Proposed scientific and engineering implications

1. Residual-Cr(VI) isotherm constants should be reported, and pooled, as apparent removal capacities. Mechanistic claims require Class A evidence (§4.4).
2. ML models of Cr(VI) removal should be validated by leave-one-study-out (or family-out) splits, with study-level uncertainty, before they are used for new materials.
3. Treatment designs for a new carbon should be based on calibrated upper bounds for both residual Cr(VI) and total Cr, and should be confirmed experimentally.
4. CASSIA's schema, evidence classes and extraction protocol are proposed as a basis for a shared, provenance-aware chromium-removal dataset. The empirical application will be reported separately.

## Data and code availability

All code, the three PHREEQC databases (with checksums), the virtual corpora, the result tables and the figure scripts are in the `cassia/` directory of the repository `karthice1/adsorption` (branch `claude/exciting-hawking-ayadsw`; CASSIA 0.2.0). Every number in this paper is regenerated by `python scripts/run_v2.py --stage all`. The master-dataset template, schema, study registry, search log and PRISMA flow templates are in `cassia/data/`. A versioned release with a Zenodo DOI will be created before publication [DOI to be inserted]. No experimental data were generated, and no published experimental data were analysed, in this study.

## Funding

This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors. [Author to confirm.]

## CRediT authorship contribution statement

**Karthik [Surname]:** Conceptualization, Methodology, Software, Formal analysis, Validation, Visualization, Writing – original draft, Writing – review & editing. [Add co-authors and roles as appropriate.]

## Declaration of competing interest

The author declares no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

## Declaration of generative AI and AI-assisted technologies in the manuscript preparation process

During the preparation of this work the author used Claude (Anthropic) to develop and run analysis code and to draft and revise the manuscript text; its use in the research process is described in Section 2.12. After using this tool, the author reviewed and edited the content as needed and takes full responsibility for the content of the published article. [The author should confirm this wording against the target journal's current policy.]

## References

Allison, J.D., Brown, D.S., Novo-Gradac, K.J., 1991. MINTEQA2/PRODEFA2, A Geochemical Assessment Model for Environmental Systems: Version 3.0 User's Manual. EPA/600/3-91/021. U.S. Environmental Protection Agency, Athens, GA.

Amvro23, 2022. adsorption: a Python package for isotherm, kinetic and dynamic adsorption models (MIT licence). https://github.com/amvro23/adsorption [accessed 30 October 2026].

Ball, J.W., Nordstrom, D.K., 1998. Critical evaluation and selection of standard state thermodynamic properties for chromium metal and its aqueous ions, hydrolysis species, oxides, and hydroxides. J. Chem. Eng. Data 43, 895–918. https://doi.org/10.1021/je980080a

Barber, R.F., Candès, E.J., Ramdas, A., Tibshirani, R.J., 2021. Predictive inference with the jackknife+. Ann. Stat. 49, 486–507.

Breiman, L., 2001. Random forests. Mach. Learn. 45, 5–32.

Burnham, K.P., Anderson, D.R., 2002. Model Selection and Multimodel Inference: A Practical Information-Theoretic Approach, 2nd ed. Springer, New York.

Chen, T., Guestrin, C., 2016. XGBoost: A scalable tree boosting system. In: Proc. 22nd ACM SIGKDD Int. Conf. Knowledge Discovery and Data Mining, pp. 785–794.

Dong, X., Ma, L.Q., Li, Y., 2011. Characteristics and mechanisms of hexavalent chromium removal by biochar from sugar beet tailing. J. Hazard. Mater. 190, 909–915.

Duan, N., 1983. Smearing estimate: a nonparametric retransformation method. J. Am. Stat. Assoc. 78, 605–610.

Dunn, R., Wasserman, L., Ramdas, A., 2023. Distribution-free prediction sets for two-layer hierarchical models. J. Am. Stat. Assoc. 118, 2491–2502. [volume/pages to be verified]

Gelman, A., Hill, J., 2007. Data Analysis Using Regression and Multilevel/Hierarchical Models. Cambridge University Press, Cambridge.

Iacomi, P., Llewellyn, P.L., 2019. pyGAPS: a Python-based framework for adsorption isotherm processing and material characterisation. Adsorption 25, 1533–1542.

Kotaś, J., Stasicka, Z., 2000. Chromium occurrence in the environment and methods of its speciation. Environ. Pollut. 107, 263–283.

Mohan, D., Pittman, C.U., 2006. Activated carbons and low cost adsorbents for remediation of tri- and hexavalent chromium from water. J. Hazard. Mater. 137, 762–811.

Mohan, D., Rajput, S., Singh, V.K., Steele, P.H., Pittman, C.U., 2011. Modeling and evaluation of chromium remediation from water using low cost bio-char, a green adsorbent. J. Hazard. Mater. 188, 319–333.

Page, M.J., McKenzie, J.E., Bossuyt, P.M., et al., 2021. The PRISMA 2020 statement: an updated guideline for reporting systematic reviews. BMJ 372, n71.

Park, D., Yun, Y.-S., Park, J.M., 2004. Reduction of hexavalent chromium with the brown seaweed *Ecklonia* biomass. Environ. Sci. Technol. 38, 4860–4864.

Parkhurst, D.L., Appelo, C.A.J., 2013. Description of input and examples for PHREEQC version 3 — A computer program for speciation, batch-reaction, one-dimensional transport, and inverse geochemical calculations. U.S. Geological Survey Techniques and Methods, book 6, chap. A43.

Pedregosa, F., Varoquaux, G., Gramfort, A., et al., 2011. Scikit-learn: Machine learning in Python. J. Mach. Learn. Res. 12, 2825–2830.

Rai, D., Sass, B.M., Moore, D.A., 1987. Chromium(III) hydrolysis constants and solubility of chromium(III) hydroxide. Inorg. Chem. 26, 345–349.

Ramirez, A., et al., 2020. Experimental and computational data set on adsorption of Cr(VI) from water using an activated carbon. Data Brief 29, 105292. https://doi.org/10.1016/j.dib.2020.105292 [full author list to be verified]

Ran, Y.A., Sharma, S., Balestra, S.R.G., Li, Z., Calero, S., Vlugt, T.J.H., Snurr, R.Q., Dubbeldam, D., 2024. RASPA3: A Monte Carlo code for computing adsorption and diffusion in nanoporous materials and thermodynamics properties of fluids. J. Chem. Phys. 161, 114106.

Rasmussen, C.E., Williams, C.K.I., 2006. Gaussian Processes for Machine Learning. MIT Press, Cambridge, MA.

Raue, A., Kreutz, C., Maiwald, T., Bachmann, J., Schilling, M., Klingmüller, U., Timmer, J., 2009. Structural and practical identifiability analysis of partially observed dynamical models by exploiting the profile likelihood. Bioinformatics 25, 1923–1929.

Roberts, D.R., Bahn, V., Ciuti, S., et al., 2017. Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. Ecography 40, 913–929.

Simonin, J.-P., 2016. On the comparison of pseudo-first order and pseudo-second order rate laws in the modeling of adsorption kinetics. Chem. Eng. J. 300, 254–263.

Tran, H.N., You, S.-J., Hosseini-Bandegharaei, A., Chao, H.-P., 2017. Mistakes and inconsistencies regarding adsorption of contaminants from aqueous solutions: A critical review. Water Res. 120, 88–116.

Wang, J., Guo, X., 2020. Adsorption kinetic models: Physical meanings, applications, and solving methods. J. Hazard. Mater. 390, 122156.

Wilson, E.B., 1927. Probable inference, the law of succession, and statistical inference. J. Am. Stat. Assoc. 22, 209–212.

Software versions: CASSIA 0.2.0; phreeqpython 1.6.2 (Vitens; https://github.com/Vitens/phreeqpython); PHREEQC databases minteq.v4.dat, minteq.dat and llnl.dat (USGS; revisions in Table S2); adsorption 0.1.0 (Amvro23); pyGAPS 4.3.1.dev693; NumPy 2.5.3; SciPy 1.18.1; pandas 3.0.5; scikit-learn 1.9.1; XGBoost 3.4.1; matplotlib 3.11.2; Python 3.13.17.
