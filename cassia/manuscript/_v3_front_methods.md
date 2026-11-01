# CASSIA: A Provenance-Aware Computational Framework and Known-Truth Benchmark for Evidence-Constrained Analysis of Cr(VI) Removal by Carbonaceous Materials

## Highlights

<!--HIGHLIGHTS-->

## Abstract

<!--ABSTRACT-->

**Keywords:** hexavalent chromium; carbonaceous adsorbents; known-truth simulation benchmark; PHREEQC; mechanism identifiability; leave-one-study-out validation; conformal prediction; data provenance

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

**Scope: a methodological benchmark, not an experimentally validated framework.** This paper is a computational-methodology study. It establishes, against a known truth, what the methods in CASSIA can and cannot recover and how reliable their outputs are. It does not establish that CASSIA interprets real chromium-removal experiments reliably; that requires the empirical application, which is defined (schema, evidence classes, extraction protocol, the open-access seed dataset of Ramirez et al., 2020) but has not been carried out (Section 3.8). No result in this paper is an empirical finding about a real carbon, and the framework should not be described as experimentally validated until that application is reported.

## 2. Methods

### 2.1 Endpoints, evidence classes and terminology

Three endpoints are computed and never merged:
- **Apparent Cr(VI) uptake**, q_app,VI = (C₀ − C_e,VI)/D. It is always computable.
- **Total-aqueous-Cr uptake**, q_tot,aq = (C₀ − C_e,tot)/D. It requires a total-Cr analysis (AAS or ICP). It is not adsorption either, because it includes Cr(III) retained or precipitated on the solid.
- **Aqueous Cr(III) index**, Σ(C_e,tot − C_e,VI)/Σ(C₀ − C_e,VI). It is a lower bound on the share of Cr(VI) disappearance due to reduction.

Evidence classes are assigned **per observation** (`cassia.dataset.classify`), because a single paper may report total Cr or solid-phase data for some experiments only. Oxidation-state identification and chromium mass balance are treated as different analytical tasks:
- **Class C:** residual Cr(VI) only. Solid-phase data without a total-aqueous-Cr analysis also count as C.
- **Class B:** Cr(VI) and total aqueous Cr.
- **Class A1:** B plus solid oxidation-state evidence, a final or controlled pH, and a reported electrolyte. XPS and XANES are recorded separately because their sampling depths and artefacts differ.
- **Class A2:** A1 plus a solid Cr inventory (digestion and total-Cr analysis) and phase or coordination evidence distinguishing retained from precipitated Cr(III).

A1 establishes that Cr(III) is present on the solid; only A2 quantifies how much Cr sits where.

A study-level label (the best class present) is kept for reporting only. Class C observations may be labelled only as *apparent* Cr(VI) removal.

Negative apparent uptake (C_e,VI > C₀ within analytical error) is retained in the data, flagged, and censored at q = 0 for model targets. Values below the limit of detection (LOD) enter as LOD/2 and are flagged.

### 2.2 Software architecture and release

CASSIA 0.3.2 is organised in six layers (Table 1). Analyses are run by scripts that are listed, with their outputs, in Table 2.

**Table 1.** CASSIA layers, modules and responsibilities.

| Layer | Module | Responsibility |
|---|---|---|
| 1 Data | `cassia.dataset` | 37-field schema with DOI and source location per row; validation and QC flags; endpoints; observation-level evidence classes C, B, A1, A2 |
| 2 Chemistry | `cassia.chemistry` | PHREEQC adapter (phreeqpython 1.6.2): Cr(VI) speciation, the 4-D speciation grid, redox scenarios, Cr(III) saturation indices, background electrolytes, database provenance |
| 3 Simulator | `cassia.simulate`, `cassia.measurement` | Known-truth multi-mechanism batch model (scenarios A–E), mass and electron balances; analytical and solid-phase measurement models (XPS, XANES, inventory, phase evidence); parameter table (Table S1) |
| 4 Identifiability | `cassia.identifiability` | Fisher information by measurement design (D1–D5), model-free detection tests, calibrated rejection and equivalence tests, profile likelihood, metadata-bias propagation |
| 5 Models and uncertainty | `cassia.fitting`, `cassia.validation` | Nonlinear isotherm/kinetic selection; RF, XGBoost and GPR under four validation schemes; study-level conformal prediction (M0, M1, M2, M2s) |
| 6 Decision support | `scripts/run_v2.py` (design stage), `scripts/run_v3.py` (seqver stage) | Two-endpoint dose/pH design under point, marginal and joint rules; sequential verification with simple and error-controlled acceptance |
| auxiliary | `cassia.pores` | pyGAPS 4.x adapter: BET area (Rouquerol criteria) from published N₂ isotherms |

**Table 2.** Scripts that regenerate the results (run in this order by `scripts/reproduce_all.sh`).

| Script | Outputs |
|---|---|
| `scripts/run_v2.py --stage all` | `results/v2/`: chemistry (Tables S2–S4), mechanism benchmark (Table 4, Fig. 2), ML validation (Table 6, Fig. 4), interval methods M0–M2 (Table 7, Table S5), treatment design (Table 8), missing chemistry (Fig. 8), Fig. 1 |
| `scripts/run_v3.py --stage chem3 … seqver` | `results/v3/`: electrolyte sensitivity (Table S10), identifiability and equivalence tests (Table 5, Table S8), M2s (Table 7b), ML across scenarios (Table S6), censoring (Table S7), learning curve, sequential verification (Table 8b) |
| `scripts/run_v3.py --stage figures` | Figs. 3, 5, 6, 7 (must run after `run_v2.py`, whose older versions of these files it replaces) |
| `scripts/typeI_validation.py` | independent negative-control cohort; type-I errors with Wilson intervals (Table S14) |
| `scripts/seqver_oc.py` | operating characteristic of the verification rules (Table S9) |
| `scripts/bench_partition.py` | true partition by scenario; Monte Carlo precision (Tables S11, S12) |
| `scripts/endpoint_validation.py` | out-of-study accuracy of the four endpoint models (Table S13) |
| `scripts/m2s_family.py` | M2s by material family (Table S5b) |
| `scripts/verification_table.py` | accepted designs and unsafe rates with explicit denominators (Table 8c) |
| `scripts/make_table_s1.py` | Table S1, generated from the implementation |
| `scripts/make_supp_figs.py` | Supplementary Fig. S1 |
| `manuscript/make_supplement_v3.py`, `manuscript/assemble_v3.py` | Supplementary Material and manuscript text with the code commit |
| `scripts/make_release_manifest.py` | `release/MANIFEST.sha256`, `release/RELEASE.md` |

**Reuse of existing packages.**
- The isotherm and kinetic equations are taken from the `adsorption` package (Amvro23, MIT licence). Its gas-phase concentration transform and its fixed-temperature Temkin form are bypassed for aqueous data.
- pyGAPS (Iacomi and Llewellyn, 2019) is used only for gas-isotherm descriptors.

**Release contents.** The release includes:
- the MIT licence, a `CITATION.cff` file and a changelog;
- a pinned environment (`requirements-lock.txt`, Python 3.13.17);
- SHA-256 checksums and `$Id` revision lines for the bundled databases;
- an end-to-end example workflow for user data (`examples/user_data_workflow.py`);
- a release record (`release/RELEASE.md`, `release/MANIFEST.sha256`) with SHA-256 checksums of every code, database, data, result and figure file, the code commit hash and the exact commands;
- 40 automated tests, covering mass and electron conservation for every scenario, saturation-index statuses, interpolation accuracy, the conformal constructions and a negative control for the equivalence test, and agreement of Table S1 with the constants used by the code.

Every number, table and figure in this paper is regenerated by the scripts of Table 2, run in order by `scripts/reproduce_all.sh` (about 9–11 h on 4 cores). All random seeds are fixed in the scripts.

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

**Cr(III) solubility.** Non-convergence of the PHREEQC phase-boundary option does not show that no solubility constraint exists, so it is not used. Instead, the saturation index of the Cr(III) hydroxide (Cr(OH)₃(am) in minteq.v4; Cr(OH)₃(A) in minteq) is evaluated directly for solutions of prescribed dissolved Cr(III). Saturation is located by bisection on log₁₀ C between 10⁻⁶ and 5,000 mg L⁻¹. Every evaluation returns a status:
- `bracketed` (with the SI at the solution);
- `undersaturated_to_cmax`: SI < 0 even at 5,000 mg L⁻¹, a direct thermodynamic statement;
- `supersaturated_at_cmin`;
- `speciation_failed`: reported, never interpreted.

`llnl.dat` contains no Cr(III) hydroxide, and its only Cr(III) solid (eskolaite) is written through CrO₄²⁻ and O₂. A Cr(III)-only solubility is therefore undefined there and is reported as not computable.

**External validation.** An independent check against the critically evaluated constants of Ball and Nordstrom (1998) and Rai et al. (1987) requires their tabulated values. These could not be retrieved from the build environment. CASSIA accepts a reference-constant table for this comparison (Section 4.7).

**Background electrolyte.** The electrolyte is NaCl unless stated otherwise. Its effect on Cr(VI) speciation was checked against NaNO₃ and Na₂SO₄ (I = 0.01 and 0.1 M; §3.1). Adsorption competition by these anions on the carbon surface is not represented.

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
- *Families and family-level parameters.* There are six families: chemically and physically activated carbons, wood and manure biochars, hydrochars, and biomass-derived porous carbons. Each has physically plausible ranges of S_BET and pH_pzc. The family-level means of reductant capacity, leachable reductant and specific capacity are redrawn for every corpus from a common hyper-prior, so which family is most reductive is itself random across replicates. No family is therefore most reductive by construction, and family rankings of reduction share in the outputs are not imposed by fixed inputs.
- *Experiments per study.* Each virtual study contains an 8-point isotherm (C₀ = 5–300 mg L⁻¹ at an optimum pH of 1.8–3.5), a 6-point pH edge (pH 2–8), an 8-point kinetic series and a 6-point dose series (0.5–8 g L⁻¹ at C₀ = 20 mg L⁻¹).
- *Measurement model.* Cr(VI) and total Cr carry 3% relative error, detection limits and a study-level calibration bias.
- *Solid-phase measurements* follow technique-specific models (`cassia.measurement`; Table S1):
  - **XPS** is surface-sensitive. The reported Cr(III) fraction is f_s + b(1 − f_s), with surface fraction f_s = e·f/(e·f + 1 − f), where f is the true fraction, e a per-study surface-enrichment factor (log-normal, SD 0.4) and b ~ U(0, 0.15) beam- or vacuum-induced reduction of Cr(VI) during analysis; the fitting noise is 0.05.
  - **XANES** is bulk-sensitive: no enrichment, b ~ U(0, 0.03), noise 0.03.
  - **Solid Cr inventory** is measured with 5% relative error.
  - **Phase evidence** gives the precipitated share of solid Cr(III) with an absolute error of 0.10.

  These are assumptions encoding qualitative behaviour, not instrument specifications.
- *Evidence classes and reporting.* Classes are assigned at study level (C 0.40, B 0.30, A1 0.20, A2 0.10; A1 split equally between XPS and XANES) and then resolved per observation. Ionic strength and final pH are reported with class-dependent probabilities.

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

**Number of replicates.** The 100 replicate corpora per scenario were checked by a Monte Carlo precision analysis (`scripts/bench_partition.py`, Table S12). For each headline rate, the between-replicate SD s gives MCSE = s/√R and the number of replicates R = (1.96 s/h)² needed for a 95% half-width h. With R = 100, the MCSE was ≤ 0.010 for unconditional fit rates, ≤ 0.018 for rates conditional on reduction share > 0.2, ≤ 0.006 for mean reduction shares and ≤ 0.020 for Spearman correlations. The contrasts on which conclusions rest are 0.1–0.4. A half-width of 0.01 on the unconditional fit rates would need 228–357 replicates; this precision is not needed for any conclusion and was not pursued.

### 2.6 Mechanism-identifiability engine

Six cumulative measurement designs were defined:

| Design | Measurements | Evidence class |
|---|---|---|
| D1 | residual Cr(VI) at equilibrium (isotherm + pH edge) | C |
| D2 | D1 + Cr(VI) kinetics | C |
| D3 | D2 + total dissolved Cr at every point | B |
| D4x | D3 + XPS Cr(III) fraction on three solids + final pH | A1 (XPS) |
| D4n | D3 + XANES Cr(III) fraction on three solids + final pH | A1 (XANES) |
| D5 | D4n + solid Cr inventory + phase evidence on the same solids | A2 |

**Which D5 measurements enter which analysis.** The Fisher analysis and the precipitation-share estimate use all D5 measurements. The adsorption-only goodness-of-fit, rejection and equivalence tests use D5 *without* phase evidence (D4n + inventory): the phase measurement is the precipitated share of solid Cr(III), which is undefined when, as under the adsorption-only null, the solid holds no Cr(III). Phase evidence is therefore interpreted separately, as a discriminator of precipitation, and the reported rejection power and equivalence rates for D5 do not include it.

(A further design, D4u, with an unbiased solid oxidation state, is used only for the profile likelihoods in Supplementary Fig. S1.)

Four complementary analyses were applied.

1. **Local identifiability.**
   - *Sensitivities.* For each study, the sensitivity matrix J of all observables was computed by forward differences. It covers the log-transformed or pH-scale kinetic parameters, the study bias and the technique nuisances: the XPS enrichment and beam reduction, or the XANES beam reduction.
   - *Posterior.* The Fisher information JᵀWJ, where W is the inverse measurement variance, was combined with a weakly informative Gaussian prior: SD 0.5 log₁₀ units for kinetic parameters and the uniform-equivalent SDs for the nuisances. The delta method then gives Laplace-approximation posterior SDs of the reduction share and of the precipitation share.
   - *Caveat.* This analysis assumes the mechanistic model structure is known, so it is an optimistic bound.
2. **Model-free detection tests.** Each was assessed over 20 noisy replicates per study, and in the paired pure-adsorption controls to measure false alarms:
   - a one-sided z-test that dissolved Cr(III) = Σ(C_e,tot − C_e,VI) exceeds analytical noise (B);
   - XPS Cr(III) fraction > 0.20 (A1, XPS);
   - XANES Cr(III) fraction > 0.10 (A1, XANES).
3. **Adsorption-only goodness of fit and equivalence (Appendix B).**
   - *Model and controls.* A flexible empirical adsorption-only model was fitted to each design's observations by weighted least squares: a pH-dependent Sips equilibrium, a free kinetic curve, study bias, pH relaxation and, for A1/A2, a bounded measurement-induced apparent Cr(III) fraction. Every study was paired with a negative control: the same material with reduction switched off. Both arms were run at 3% and 10% analytical noise.
   - *Three outcomes per study.*
     - **rejected**, when χ² exceeds a threshold calibrated on the negative controls (their 95th percentile; type-I error estimated by split-half cross-validation);
     - **equivalent within ε**, when a noncentral-χ² equivalence test, the goodness-of-fit analogue of two one-sided tests (Schuirmann, 1987; Wellek, 2010), shows that the misfit non-centrality is below a margin λ₀(ε). The margin corresponds to a systematic discrepancy of ε = 5, 10 or 20% in every concentration (0.05 absolute for oxidation-state fractions);
     - **inconclusive**, otherwise.
   - *Terminology.* "Not rejected" is never interpreted as "equivalent".
   - *Operating characteristics.* Type-I error, power and false-equivalence rate were evaluated against the noise-free misfit of each study.
4. **Metadata bias.** The first-order bias in the inferred reduction share was propagated for misreported pH (+0.5) and ionic strength (×10), via Δθ = Σ_post JᵀW(y_true − y_assumed). These analyses use designs D1 and D3.

**Sample sizes.**
- Analyses 1 and 2: 300 studies (200 from scenario E and 100 from scenarios B–D).
- Analysis 3: 120 studies (60 from E, 20 each from B–D), each with a paired control, × 6 designs × 2 noise levels.

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
- *GPR:* evaluated on three corpora because of its cost, and treated as exploratory. Its ARD kernel hyperparameters were optimised on a seeded 250-row subsample of each training fold and then held fixed while the GP was refitted on all training rows, so no test information is used.

**Reporting.** Each corpus reports pooled R², RMSE, MAE and median absolute error, with 1,000-replicate study-level (cluster) bootstrap intervals and per-study errors. Across corpora the distribution of each metric is reported, and model differences are paired by corpus. Between-study variability is quantified by a random-intercept decomposition of out-of-study residuals (τ², σ², ICC; Gelman and Hill, 2007). A learning curve (12–96 studies; eight replicates each) estimates how many independent studies a transferable model needs.

**Additional analyses.**
- *Mechanism scope.* The main ML analysis uses scenario E. To test whether its conclusions depend on the mechanism, RF and XGBoost (random, LOSO, LOFO) and the interval methods were repeated on five corpora of each scenario A–E.
- *Cross-mechanism transfer.* Models trained on scenario-E corpora were also tested on corpora generated under scenarios A–D.
- *Censoring and LOD treatment.* Four treatments were compared on five corpora (RF, LOSO): censoring at zero with LOD/2 (baseline), dropping censored rows, an inverse-hyperbolic-sine transform that keeps negative values, and LOD instead of LOD/2.

### 2.8 Study-level conformal prediction

Four interval constructions were compared. All are calibrated on out-of-study residuals, with the study as the unit of resampling. They differ in the unit of the conformity score and in their guarantee:

| Method | Conformity score | Procedure class | Coverage statement |
|---|---|---|---|
| M0 | one score per observation, pooled | inner grouped K-fold residuals (cross-validation residual method; neither split nor cross-conformal in the formal sense) | none for a new study |
| M1 | one score per observation, weighted 1/n_j per study | study-level CV+ (Barber et al., 2021) with group weights | none claimed; targets coverage of a random observation of a new study |
| M2 | one score per study (within-study 0.9-quantile) | study-level CV+ | none claimed (Appendix A, Remark 5) |
| M2s | one score per study (within-study 0.9-quantile) | split conformal over studies | exact: P(≥ ⌈0.9 n⌉ of a new study's n observations covered) ≥ 1 − α (Appendix A) |

Only M2s carries a proven study-level guarantee. The coverage of M0–M2 for unseen studies is established empirically, by outer leave-one-study-out evaluation, and reported per study.


- **M0 – pooled residuals (the conventional approach).** The half-width is the ⌈(n+1)(1−α)⌉/n quantile of pooled absolute residuals from an inner grouped 10-fold CV. Residuals within a study are dependent, so M0 has no formal guarantee for a new study.
- **M1 – study-level CV+ with study-weighted scores.**
  - *Construction:* the CV+ interval of Barber et al. (2021), formed with study-level folds (K = 10). Each calibration study contributes total weight 1 (weight 1/n_j per observation), as in the weighted pooling of Dunn et al. (2023), and the test study carries one further unit of weight.
  - *Target:* marginal coverage of a random observation from a new study.
- **M2 – study-level CV+ with study scores (empirical).**
  - *Score:* s_j is the ⌈0.9 n_j⌉-th smallest absolute residual of study j.
  - *Interval:* bounded by the ⌊α(J+1)⌋-th and ⌈(1−α)(J+1)⌉-th order statistics of the fold-model predictions ∓ s_j.
  - *Status:* no coverage guarantee is claimed. The CV+ bound is proven for scalar scores of exchangeable points, and its extension to within-study quantile scores is not established (Appendix A).
- **M2s – study-level split conformal (exact).**
  - *Construction:* the training studies are split at random into halves, one for model fitting and J for calibration. Each calibration study is scored by s_j as above, and the half-width is the ⌈(1−α)(J+1)⌉-th smallest score.
  - *Guarantee (Appendix A):* under exchangeability of studies, with probability ≥ 1 − α at least ⌈0.9 n⌉ of a new study's n observations are covered. The guarantee is marginal over the study population (families, study sizes), not conditional on a family.
  - *Requirement:* J ≥ (1 − α)/α calibration studies. With fewer, the interval is unbounded.

**Evaluation.** Each construction was evaluated by outer LOSO on 10 corpora for all four targets, at nominal levels 0.50–0.95. There were 29 training studies per test study; for M2s, about 14 of these were used for calibration. Reported for each:
- mean, median, 10th-percentile and worst per-study coverage;
- the fraction of studies covered at the nominal level;
- the observation-weighted coverage;
- the mean width;
- the number of calibration studies.

### 2.9 Treatment-design layer

**Models.** Residual Cr(VI) and residual total Cr are modelled directly (Section 2.7), each with its own M1 interval. Modelling effluent indirectly through q is ill-conditioned when removal approaches completeness.

**Design problem.**
- *Grid:* pH 2–7 in 0.5 steps and dose 0.5–8 g L⁻¹ (the training range), 24 h contact.
- *Rule:* choose the lowest dose meeting both limits.
- *Targets:* S1, C₀ = 20 mg L⁻¹ with limits Cr(VI) ≤ 0.5 and total Cr ≤ 2 mg L⁻¹; S2, C₀ = 5 mg L⁻¹ with limits 0.1 and 0.5 mg L⁻¹.

**Decision rules.**
- *Point:* point predictions.
- *Marginal 90%:* each endpoint's upper bound at α = 0.10.
- *Joint 90%:* Bonferroni, each endpoint at α = 0.05, so that the joint miscoverage is ≤ 0.10 by the union bound.
- *Joint 90% / 80% + 4 confirmatory batches:* four batches on the new material (pH 3, doses 0.5–4 g L⁻¹, 3% analytical noise) estimate its study offset for each endpoint. The centre is shifted by that offset. The half-width is the study-weighted (1 − α) quantile of offset-adjusted out-of-study residuals of the training studies, with α = 0.05 or 0.10 per endpoint (heuristic; it assumes the confirmatory batches are exchangeable with the study's other observations).
- *Point + verify (model-guided, test-confirmed):* the point model proposes the pH and the lowest dose. Confirmatory batches are run at that design, and, if they fail, the dose is escalated along the grid for at most four rounds. Two acceptance rules were compared:
  - **(i) Simple rule:** duplicates; accept when both mean measured concentrations are ≤ 90% of their limits.
  - **(ii) error-controlled rule:** triplicates; accept when the one-sided upper confidence bound of each endpoint is ≤ its limit. The bound is mean × exp[z·√(σ_r²/n + σ_b²)], with σ_r = 0.03 random and σ_b = 0.05 calibration uncertainty. z is the 1 − α/(2R) normal quantile, with α = 0.05 and R = 4 rounds (Bonferroni over two endpoints and four looks).

  The design cohort uses the simple rule with one realised verification sequence per material. Both rules were then compared on an independent verification cohort of 300 materials per target, drawn from the same hyper-priors: the probability of accepting a non-compliant design was estimated with 1,000 simulated verification sequences per material, target and rule. Two analytical-error models were used: calibration bias shared by both endpoints with correlated replicate errors (ρ = 0.5), and independent biases and errors.

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

**Safeguards.** The scientific design and all methodological choices were specified by the author. The correctness checks built into the work are:
- unit tests for conservation laws and numerical constructions;
- negative controls (scenario A);
- cross-checks of interpolated against direct PHREEQC calculations;
- full regeneration of every reported number from seeded scripts.

[The author should add a statement of their own verification of code and results here before submission.]
