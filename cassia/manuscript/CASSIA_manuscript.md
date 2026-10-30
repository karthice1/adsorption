# CASSIA: A Data-Provenance-Aware, Speciation-Constrained Computational Framework for Distinguishing Cr(VI) Adsorption, Reduction and Precipitation in Carbonaceous Wastewater Treatment Systems

**Author:** Karthik [surname, affiliation and ORCID to be completed by the author]
**Corresponding author:** karthikresearchworks@gmail.com

**Article type:** Research article (computational methods)

---

## Highlights

- A decrease in measured Cr(VI) is not evidence of adsorption. CASSIA keeps the two apart through provenance-tagged endpoints and evidence classes.
- Two PHREEQC databases agree on Cr(VI) speciation but differ in redox and Cr(III)-solid thermodynamics.
- In known-truth simulations, isotherm fits with R² ≥ 0.95 conceal reduction and overstate Langmuir q_max.
- Leave-one-study-out validation shows how much random splits overstate the transferability of machine-learning models.
- Study-aware conformal intervals restore nominal coverage where model-internal intervals under-cover.

## Abstract

<!--ABSTRACT-->

**Keywords:** hexavalent chromium; activated carbon; biochar; PHREEQC; adsorption isotherm identifiability; leave-one-study-out validation; conformal prediction; data provenance

---

## 1. Introduction

Hexavalent chromium is a priority contaminant in effluents from electroplating, leather tanning, pigment and wood-preservation industries. Its toxicity and mobility differ sharply from those of Cr(III) (Kotaś and Stasicka, 2000). Carbonaceous materials are the most studied low-cost sorbents for Cr(VI), and the literature now contains many batch studies on activated carbons, biochars and biomass-derived porous carbons (Mohan and Pittman, 2006; Mohan et al., 2011; Dong et al., 2011). Most of these studies report isotherm and kinetic constants, and a growing number train machine-learning models on the pooled values.

The difficulty is that Cr(VI) on carbon is a redox-active system. At low pH, oxygen-containing functional groups and the carbon matrix reduce Cr(VI) to Cr(III). Part of the Cr(III) is released to solution, part stays bound to the surface, and at higher pH it can precipitate as Cr(OH)₃ (Park et al., 2004; Dong et al., 2011). The quantity most often measured, residual Cr(VI) by the 1,5-diphenylcarbazide method, therefore records the *disappearance* of Cr(VI), not its adsorption. When q_e = (C₀ − C_e)/m is computed from residual Cr(VI) and fitted with Langmuir, Freundlich or pseudo-second-order equations, the fitted constants mix adsorption, reduction and precipitation in proportions the analysis cannot see. Critical reviews have documented related problems in aqueous adsorption modelling: linearised fitting, R² as the only selection criterion, and misread kinetic models (Tran et al., 2017; Simonin, 2016; Wang and Guo, 2020).

Machine-learning studies add a further risk. Rows from one publication share an adsorbent, a laboratory and an analytical protocol. Random train/test splits therefore leak study identity into the test set and inflate apparent accuracy (Roberts et al., 2017). A model can reach a high random-split R² while being unable to predict an unseen material.

This work builds and tests the infrastructure needed to ask whether aqueous speciation, adsorption-model selection and uncertainty-aware learning can separate genuine adsorption from apparent Cr(VI) removal. Four research questions are posed:

- **RQ1.** Does PHREEQC-derived Cr(VI) speciation improve predictive performance over raw pH and composition descriptors?
- **RQ2.** Can conventional isotherm and kinetic fitting distinguish adsorption from apparent removal?
- **RQ3.** Do models generalise across studies and material families, and how much do random splits overstate this?
- **RQ4.** How much uncertainty comes from unreported solution chemistry and from between-study heterogeneity, and can prediction intervals be made reliable?

RQ2–RQ4 cannot be answered with published data alone, because no published study reports the true adsorption/reduction split for every observation. We therefore pair a PHREEQC chemistry layer with a *known-truth* simulation study. Virtual batch experiments are generated from a mechanistic model that couples speciation-dependent adsorption, acid-promoted surface reduction and Cr(III) partitioning. They are then analysed exactly as published data would be. This design tests whether a method *can* recover the mechanism when the truth is known. It is a necessary condition for trusting the method on real data, but it does not replace real data. The multi-study literature dataset is defined here (schema, evidence classes, extraction protocol), but its extraction is reported as ongoing, and no empirical conclusion about real carbons is drawn in this paper.

The contribution is deliberately methodological. CASSIA does not claim a new adsorption mechanism, and the existing packages it uses (phreeqpython, adsorption, pyGAPS) are used for what they were designed for.

## 2. Methods

### 2.1 Software audit and roles

Four open-source code bases were audited for fitness of purpose (Table 1).

**Table 1.** Roles of the audited code bases in CASSIA.

| Code base | Licence | Role in CASSIA | Audit finding that shaped the design |
|---|---|---|---|
| phreeqpython (PHREEQC 3) (Parkhurst and Appelo, 2013) | Apache-2.0 | Aqueous Cr chemistry (core) | Bundled databases (`phreeqc.dat`, `wateq4f_PWN.dat`, `vitens.dat`, `pitzer.dat`) contain no Cr species. Cr-bearing databases from the USGS distribution must be supplied and validated. The `change_ph` helper relies on a `Fix_pH` phase absent from `minteq.v4.dat`. |
| adsorption (amvro23) | MIT | Isotherm and kinetic equations (core) | `Isotherms.set_inlet` applies a gas-phase transform (P·M_r/RT) to x. The Temkin equation fixes T = 298.15 K. Both are bypassed for aqueous data, and the equations are reused directly. |
| pyGAPS (Iacomi and Llewellyn, 2019) | MIT | BET recalculation from published N₂ isotherms (secondary) | Designed for gas/vapour isotherms. Used only for textural descriptors, never for aqueous Cr data. `area_BET_raw` expects the cross-section in nm². |
| RASPA3 (Ran et al., 2024) | see upstream | Not used in this release | No aqueous Cr(VI)–carbon surface-reaction model, and its reaction-ensemble functionality is experimental. Deferred (Section 4.4). |

### 2.2 Chromium chemistry layer

PHREEQC 3 was driven through phreeqpython 1.6.2 with two databases distributed by the USGS: `minteq.v4.dat` (MINTEQA2 v4 lineage; Allison et al., 1991) and `llnl.dat` (LLNL `thermo.com.V8.R6.230`). Other databases were rejected programmatically. Both contain the Cr(VI) protonation and dimerisation equilibria, Cr(III) hydrolysis species and the Cr(VI)/Cr(III) redox couple.

*Cr(VI) speciation.* Cr was entered as `Cr(6)` with one K⁺ per Cr, representing K₂Cr₂O₇ stock. PHREEQC does not re-equilibrate valence states entered separately in an initial solution, so this represents a freshly prepared Cr(VI) solution regardless of the default pe. Background electrolyte was NaCl at the stated ionic strength. pH was imposed by charge balance on Cl⁻ (acid addition) or, where the chromate charge exceeded the cation pool, on Na⁺ (base addition). Fractions are reported per Cr atom, so Cr₂O₇²⁻ counts twice.

*Redox scenarios.* Because batch papers almost never report Eh, the Cr(VI)/Cr(III) partition was computed as an explicit scenario: total Cr entered without a valence state, at an assumed pe. The grid was restricted to the stability field of water (pe ≤ 20.75 − pH at 25 °C). Outside it, PHREEQC does not converge and the results have no physical meaning. For each pH, the pe at which Cr(VI) and Cr(III) are equally abundant (pe₅₀) was obtained by interpolation. Eh = 0.05916·pe V.

*Cr(III) solubility.* The dissolved Cr(III) in equilibrium with Cr(OH)₃(am) or Cr₂O₃ was computed with the SOLUTION phase-boundary option (concentration adjusted to SI = 0 at fixed pH) in `minteq.v4.dat`. In `llnl.dat` the only Cr(III) solid (eskolaite) is written in terms of CrO₄²⁻ and O₂(aq), so a Cr(III)-only solubility is undefined without fixing the redox state. This asymmetry is reported as a database finding, not hidden. Below pH ≈ 4.5 the phase-boundary iteration does not converge for Cr(OH)₃(am), which indicates no solubility constraint at realistic concentrations.

*Speciation features for learning.* For each experimental condition (pH, C₀, I), CASSIA computes f(HCrO₄⁻), f(CrO₄²⁻), f(Cr₂O₇²⁻), the mean Cr(VI) charge per Cr atom, and the PHREEQC ionic strength. pH − pH_pzc is added as a surface-charge proxy.

### 2.3 Data provenance, endpoints and evidence classes

The master dataset schema (`data/cassia_schema.csv`) has 31 fields. Every row carries a DOI, a source location (table, figure or SI file) and an extraction method. Three endpoints are kept separate:

- q_app,VI = (C₀ − C_e,VI)/D, apparent Cr(VI) disappearance, always computable;
- q_total = (C₀ − C_e,total)/D, total-Cr uptake, only when total Cr was measured;
- the aqueous Cr(III) index, Σ(C_e,total − C_e,VI)/Σ(C₀ − C_e,VI), the share of Cr(VI) disappearance recovered as dissolved Cr(III).

Studies are assigned to evidence classes automatically (`cassia.dataset.classify`). **Class A** needs total and Cr(VI)-specific analyses, solid-phase Cr speciation (XPS or XANES), a final or controlled pH, and a reported background electrolyte. **Class B** needs either total Cr or solid-phase speciation. **Class C** has residual Cr(VI) only. For Class C rows the only permitted endpoint label is `apparent_CrVI_disappearance`. Unreported Eh, ionic strength, alkalinity and major ions are never imputed into the data file. They enter only as labelled scenario ranges (Section 2.7). The PRISMA 2020-aligned search, screening and extraction protocol (Page et al., 2021) is in `data/PROTOCOL.md`.

### 2.4 Known-truth virtual batch experiments

Each virtual study draws one material from one of three families (activated carbon, biochar, biomass-derived porous carbon). Family-specific ranges cover BET area, point of zero charge and reductant capacity, the latter being largest for biochars, consistent with their higher content of oxidisable organic matter. The batch model is a set of ordinary differential equations integrated with LSODA:

- Cr(VI) uptake follows a linear-driving-force approach to a Langmuir equilibrium q_eq = q_max K C/(1 + K C). The affinity is K = K₀ [0.05 + f(HCrO₄⁻ + Cr₂O₇²⁻)] [0.1 + (1 + exp((pH − pH_pzc)/0.8))⁻¹]. The speciation fractions are interpolated from PHREEQC (`minteq.v4.dat`, 50 mg L⁻¹, I = 0.01 M).
- Reduction is first order in dissolved Cr(VI) and proportional to the remaining surface reductant R_s (finite capacity R₀), with acid promotion (1 + 10^(pH − pH₅₀))⁻¹ and Arrhenius temperature dependence (E_a = 35 kJ mol⁻¹).
- Reduced chromium is partitioned into a surface-bound fraction (logistic in pH) and a dissolved fraction. The dissolved fraction is capped at the PHREEQC Cr(OH)₃(am) solubility, and any excess precipitates.

Measured concentrations carry 3% relative analytical noise, a 0.02 mg L⁻¹ detection limit and a study-level multiplicative calibration bias (log-normal, σ = 0.05). Each study contains an 8-point isotherm (C₀ = 10–400 mg L⁻¹ at the study's optimum pH, 1.8–3.5), a 6-point pH edge (pH 2–8, C₀ = 100 mg L⁻¹) and an 8-point kinetic series. Dose, temperature, ionic strength and equilibration time vary between studies. The training corpus has 24 studies (8 per family, 528 observations, seed 2026). An independent external corpus of 9 studies uses a different seed (777). Mass conservation is checked by unit test.

The simulator encodes one plausible mechanistic picture, chosen to be qualitatively consistent with the literature. It is not calibrated to any material. Its purpose is to make the ground truth available, not to predict real carbons.

### 2.5 Isotherm and kinetic model benchmarking (Model A)

The Langmuir, Freundlich, Sips and Temkin isotherms, and the pseudo-first-order (PFO) and pseudo-second-order (PSO) kinetic equations, were fitted in their nonlinear forms by bounded least squares. The equations were taken from `adsorption.models`, with Temkin evaluated at the experimental temperature. Models were ranked by AICc and Akaike weights (Burnham and Anderson, 2002), with BIC also reported. Residuals were tested with Shapiro–Wilk (normality) and a Wald–Wolfowitz runs test (systematic misfit). Parameter 95% confidence intervals came from a 300-replicate residual bootstrap. For each virtual study the fitted Langmuir q_max was compared with the true adsorption capacity, and the true share of Cr(VI) disappearance caused by reduction was computed from the simulator state.

### 2.6 Machine-learning models and validation (Models B and C)

The target was q_app,VI, modelled as log₁₀(q + 0.1 mg g⁻¹) and back-transformed for metrics. Two feature sets were compared. The *base* set (C₀, pH, D, t, T, I, S_BET, V_p, d_p, pH_pzc) is what a paper typically reports. The *chem* set adds the PHREEQC speciation features (Section 2.2). Three learners were used: random forest (Breiman, 2001; 400 trees), XGBoost (Chen and Guestrin, 2016; 500 rounds, depth 4) and Gaussian-process regression with a Matérn-5/2 ARD kernel plus white noise on standardised features (Rasmussen and Williams, 2006), all via scikit-learn (Pedregosa et al., 2011).

Four validation schemes were applied to every model–feature combination: random 5-fold (diagnostic only), grouped 5-fold by study, leave-one-study-out (LOSO) and leave-one-material-family-out (LOFO). The models were then refitted on all 24 studies and scored on the 9-study external corpus. RMSE, MAE, median absolute error and R² are reported in mg g⁻¹.

*Uncertainty (Model C).* Model-internal intervals (GPR predictive SD; RF across-tree spread) were compared with study-aware split-conformal intervals (Lei et al., 2018). In each outer LOSO fold, absolute residuals were collected from an inner 3-fold grouped CV over the training studies, and the ⌈(n+1)α⌉/n quantile set the interval half-width. Because calibration residuals are always out-of-study, the intervals absorb between-study variability. Coverage was evaluated at nominal levels 0.5–0.95, with PICP₉₀ (coverage at 90%) and the average coverage error (ACE, mean |empirical − nominal|). Between-study variability was quantified by a method-of-moments random-intercept decomposition of the out-of-study residuals (τ², σ², ICC; Gelman and Hill, 2007).

### 2.7 Missing-chemistry scenarios and robust design

For the pH-edge conditions, the effect of unreported chemistry on f(HCrO₄⁻) was propagated by Monte Carlo (40 draws per condition). Ionic strength was drawn log-uniformly over 10⁻³–10⁻¹ M, and an initial-versus-final pH ambiguity was drawn uniformly over ±0.5 units.

A treatment-design exercise compared nominal and robust decision rules. For 15 new virtual materials, pH (2.0–7.0) and dose (0.5–15 g L⁻¹) were chosen to minimise dose subject to residual Cr(VI) ≤ 0.5 mg L⁻¹ and dissolved total Cr ≤ 2 mg L⁻¹, starting from C₀ = 20 mg L⁻¹. The nominal rule used random-forest point predictions of q_app,VI and q_total. The robust rule used the lower 90% study-aware conformal bound on uptake. Each recommended design was then checked against the simulator's true outcome.

### 2.8 Reproducibility

All numbers, tables and figures are produced by `scripts/run_all.py` (seeded). The chemistry layer, endpoint logic, fitting and simulator mass balance are covered by 14 unit tests. Code, databases and outputs are in the `cassia/` directory of the accompanying repository.

<!--RESULTS-->
