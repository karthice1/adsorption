# CASSIA: A Provenance-Aware Computational Framework for Evidence-Constrained Analysis of Cr(VI) Removal by Carbonaceous Materials

## Highlights

<!--HIGHLIGHTS-->

## Abstract

<!--ABSTRACT-->

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
- 34 automated tests (85% line coverage), covering mass and electron conservation for every scenario, saturation-index statuses, interpolation accuracy, the conformal constructions and a negative control for the equivalence test.

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
- *GPR:* evaluated on five corpora because of its cost.

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

**Cohorts.** Twenty independent cohorts were used, each with its own 30-study training corpus and 15 new materials drawn from the same hyper-priors (300 materials per target).

**Truth and metrics.** For every material, the simulator gave the true attainability: the lowest grid dose at which any pH meets both limits. For every rule and target the following were reported:
- unsafe recommendations (recommended but non-compliant in truth), with Wilson 95% intervals (Wilson, 1927);
- compliant designs found among attainable materials;
- infeasibility declarations for attainable materials;
- reliability among recommendations;
- the dose overhead relative to the true minimum dose.

**Sensitivity.** The operating-window sensitivity of joint-rule designs was evaluated at pH ± 0.25 and ± 0.5 and dose × 0.8 and × 1.2. Predicted and true effluent were compared on a sub-grid.

### 2.10 Propagation of unreported chemistry

**Prediction level.** For studies that did not report ionic strength or final pH, 1,000 scenarios were drawn: I log-uniform on 10⁻³–10⁻¹ M and a pH-basis ambiguity uniform on ±0.5. They were propagated through the speciation features, the q_app,VI predictions and the M1 intervals. For each observation the following were recorded: the 90% range of f(HCrO₄⁻), the spread of the point prediction, and the union interval width against the nominal width. Grouped 5-fold models were used on four corpora.

**Design level.** For six new materials per corpus, the designer's assumed I was drawn 300 times from the same range, while the true I was 0.01 M. The joint-rule design was re-optimised for each draw, and the true outcome of each recommended design was simulated. This gives distributions of recommended dose, residual total Cr and compliance probability.

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
