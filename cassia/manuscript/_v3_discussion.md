## 4. Discussion

### 4.1 What can be established from which measurements

The revised analysis gives a mechanism-resolved answer to the question in the title. It rests on synthetic evidence and should be read with that qualification.

**Terminology.** Three outcomes are distinguished throughout (Appendix B). 'Rejected' means the calibrated χ² test rejected the adsorption-only model. 'Equivalent within ε' means the equivalence test demonstrated that the adsorption-only misfit is smaller than a stated margin. 'Not rejected' alone is inconclusive and is never read as structural equivalence.

**Class C (residual Cr(VI) only).**
- Such data establish that Cr(VI) disappears under the stated conditions, and how fast.
- They do not establish why. In the benchmark, adsorption-only was *not rejected by the specified test* for 63–76% of high-reduction studies, and it was *demonstrated equivalent within 10%* for 43–65% of them (Table 4). For 90% of these studies the noise-free misfit was truly within 10%: the non-rejection reflects genuine structural equivalence of the Class C observations, not merely low power.
- Conventional isotherm statistics do not repair this. Fit quality is weakly, and sometimes positively, related to reduction; residual diagnostics flag 5–8% of fits regardless of mechanism; bootstrap intervals for q_max are precise but rarely cover the true capacity.
- The fitted q_max is an *apparent Cr(VI) removal capacity*. Under surface or adsorption-coupled reduction it exceeds the adsorption capacity in proportion to the reduction share.

**Class B (adds a total-Cr analysis).**
- Adding total Cr changes the picture qualitatively. The calibrated test then rejected adsorption-only in 80–88% of high-reduction studies, with cross-validated type-I error of 3–5%, and equivalence within 10% was demonstrated for only 12–20%.
- The aqueous Cr(III) index ranks studies by reduction share, as a lower bound.

**Class A1 and A2 (add solid-phase evidence).**
- A1 (oxidation state of the solid and final pH) raised calibrated power to 90–96%. XPS and XANES gave similar power in the fitted test, but their model-free diagnostics differed: the assumed beam reduction and surface enrichment made a simple XPS threshold false-alarm in 7.1% of pure-adsorption controls, against 1.0% for XANES (§3.3). An XPS Cr(III) signal is therefore not, on its own, evidence of chemical reduction.
- A2 (adds a solid Cr inventory and phase evidence) added little to rejection power over A1, but it was the only design that constrained the *precipitated* share (posterior SD 0.026 against 0.037–0.14). Its benefit is the partitioning of the reduced Cr, not the detection of reduction.
- These detectability statements are conditional on the assumed measurement models; with larger real XPS biases, the A1/A2 advantage would shrink.

**Recommendation.** This hierarchy is a quantitative argument, within the benchmark, for the minimum reporting set of §4.4. It also gives a criterion for pooling literature data: Class C observations may train a model of apparent removal, but they must not be pooled with Class A/B observations as evidence of adsorption capacity.

**Why the known-structure analyses look better.** The Fisher and profile analyses are much more optimistic than the equivalence test, and the contrast is instructive. If the analyst knew the exact rate laws, even Class C data would constrain the reduction share. Real analysts do not, and once the structure is unknown, the information about mechanism in Class C data largely disappears: the equivalence test shows that a flexible adsorption-only model reproduces them within a few per cent. Identifiability here is a property of the data *and* of what the analyst is entitled to assume.

### 4.2 Thermodynamic feasibility versus surface kinetics

**What thermodynamics says.** At the pH where carbons remove Cr(VI) most effectively (1.5–3.5), Cr(VI) is unstable with respect to Cr(III) at almost any Eh inside the water-stability field (§3.1). The two MINTEQ databases do not even reach a 50% Cr(VI) fraction at pH 2. Cr(III) is also undersaturated with respect to Cr(OH)₃(am) up to thousands of mg L⁻¹ below pH 4.25 (minteq.v4, I ≤ 0.01 M; 3.5 with minteq.dat).

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

4. **Scope across mechanisms.** The validation ordering (random ≫ LOSO > LOFO) held in all five mechanism scenarios, but a model trained under one mechanism lost 0.05–0.10 in R² when applied to corpora generated under solution or adsorption-coupled reduction. A literature model is therefore conditional on the mechanism mix of its training studies, which Class C data cannot reveal.
5. **The Gaussian process is exploratory.** It was run on 3 corpora only, its kernel was not tuned per corpus, and it collapsed outside random splits. No conclusion about kernel methods in general is drawn from it.

### 4.4 A minimum reporting standard for mechanistic claims

The results support five requirements for any paper that claims Cr(VI) *adsorption* by a carbonaceous material:
1. total Cr alongside Cr(VI) for at least the isotherm and pH-edge series;
2. solid-phase Cr speciation (XPS or XANES) for at least the lowest- and highest-loading solids;
3. final as well as initial pH;
4. the background electrolyte and ionic strength;
5. replicate analytical SDs.

**Why these.** Items 1–2 move a study from Class C to Class A, which is the step that made mechanisms distinguishable in the benchmark (§3.3). In the benchmark, items 1–2 raised the calibrated power to reject adsorption-only from about 35% to 80–96%. XANES is preferable to XPS where available, and XPS spectra should be acquired with short exposure and reported with beam-damage controls. Item 3 matters most among the metadata: a 0.5-unit pH-basis ambiguity biases the inferred reduction share by up to 0.16 (90th percentile), whereas a ten-fold error in ionic strength biases it by at most 0.03.

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

**RASPA3.** It was audited and not used. It offers mature Monte Carlo and molecular-dynamics capabilities, but no model of aqueous Cr(VI)–carbon surface reactions, and its reaction-ensemble functions are experimental. A narrowly defined optional module, for example the pore volume accessible to hydrated HCrO₄⁻ versus Cr₂O₇²⁻ in idealised slit pores, could inform accessibility descriptors. It could not inform reduction, which in the benchmark was the decisive process and was not identifiable from Class C data.

### 4.7 Limitations

1. **No empirical validation.** All findings in §3.2–3.7 come from simulations. The extraction of the seed dataset (Ramirez et al., 2020) and of 20–30 further studies was not carried out, because the publisher and repository hosts were not reachable from the computational environment (§3.8). This is the main limitation of the paper: the empirical application is still outstanding.
2. **One model family.** The simulator spans five mechanisms and randomised family-level priors, but all are variants of one kinetic model family. Results that depend on assumed scaling, such as the dominance of low-surface-area chars in reduction share, are properties of the simulator. Mechanisms not represented include electron shuttling by dissolved organic matter, Cr(V) intermediates, competitive anion effects and surface-site heterogeneity. Cross-mechanism transfer (§3.4) was tested only within this family.
3. **Thermodynamic validation.** Three USGS databases were compared. An independent validation against critically evaluated constants (Ball and Nordstrom, 1998; Rai et al., 1987) is outstanding: those constants could not be retrieved and verified in this revision. NaCl, NaNO₃ and Na₂SO₄ backgrounds were compared (§3.1); Ca²⁺, Mg²⁺, carbonate and dissolved organic matter were not.
4. **Assumed measurement models.** XPS is modelled with log-normal surface enrichment (SD 0.4 in log), beam reduction up to 15% and noise 0.05; XANES with beam reduction up to 3% and noise 0.03; inventories with 5% error. These distributions are plausible but not calibrated against round-robin data. The A1/A2 results depend on them.
5. **Equivalence margins are conventions.** The margins (5, 10, 20% systematic discrepancy) are expressed in units of the analytical SD and are choices, not physical constants. Equivalence is demonstrated relative to one flexible empirical model class; a more flexible class would make Class C data even less informative.
6. **Conformal guarantees.** The M2s guarantee is exact but marginal over exchangeable studies (Appendix A). A new material from an unrepresented family violates exchangeability, as the leave-one-family-out results show. M2 has no proven guarantee.
7. **Approximate procedures.** The Fisher analysis assumes known rate laws, and the profiles include the true parameters among optimiser starts. Sequential-verification error rates hold under the assumed analytical-error model (3% random, 5% calibration uncertainty).
8. **Censoring.** Negative apparent uptakes (about 5% of observations) were censored at zero in the main analysis. The conclusions were unchanged under LOD-based censoring, an asinh transform or dropping censored values (LOSO R² 0.58–0.68; mean M1 coverage 0.94; Table S7).

## 5. Conclusions

### 5.1 Demonstrated computational findings (thermodynamic calculations and known-truth simulations)

1. Three PHREEQC databases agree on Cr(VI) speciation to within 0.044 in species fraction. They differ by 63–83 mV on the Cr(VI)/Cr(III) boundary and by more than an order of magnitude on Cr(OH)₃(am) solubility at pH 4. Direct saturation indices show Cr(III) undersaturated up to 5,000 mg L⁻¹ at pH ≤ 4.25 (minteq.v4, I ≤ 0.01 M; 4.5 at I = 0.1 M) and pH ≤ 3.5 (minteq), on a 0.25-pH grid (Table S4). The constants have not been validated independently.
2. Across 15,000 virtual studies under five mechanisms, isotherm fits of residual Cr(VI) with R² ≥ 0.95 occurred in 41–52% of studies with more than 20% reduction under surface, adsorption-coupled and combined reduction. Fitted q_max then exceeded the adsorption capacity in proportion to the reduction share. The effect was weak when reduction occurred only in solution.
3. Observed through residual Cr(VI) alone, adsorption-only was not rejected by the calibrated test for 63–76% of high-reduction studies, and was demonstrated equivalent within a 10% margin for 43–65%. With total Cr, the calibrated test rejected it in 80–88% (type-I error 3–5%); with solid-phase evidence, in 90–96%. Only the A2 design constrained the precipitated share.
4. LOSO R² for apparent uptake (0.48–0.65 across the five scenarios) was far below random-split R² (0.86–0.89). Speciation features did not improve it, and a corpus of 20–30 studies cannot rank learners. With 8 replicates per size, LOSO R² rose from 0.31 (12 studies) to 0.66 (48) and then levelled off (0.66 at 96; §3.4).
5. Pooled conformal intervals reached nominal average coverage but covered only 65–75% of studies. The exact study-level split-conformal interval M2s met its guarantee, but at 3.3 log₁₀ units at 90% and unbounded at 95% with literature-scale numbers of calibration studies.
6. Point-prediction designs were unsafe for 21% (S1) and 7% (S2) of new materials, and calibrated-bound designs were rarely available. Model-guided designs with sequential confirmatory testing were accepted for 98–105 of 300 materials per target; no unsafe accepted design was observed in the simulated cohorts (Wilson 95% upper bound 1.3%), and the error-controlled rule kept the per-material false-acceptance probability ≤ 0.7%.

### 5.2 Proposed scientific and engineering implications

1. Residual-Cr(VI) isotherm constants should be reported, and pooled, as apparent removal capacities. Mechanistic claims require Class A evidence (§4.4).
2. ML models of Cr(VI) removal should be validated by leave-one-study-out (or family-out) splits, with study-level uncertainty, before they are used for new materials.
3. For a new carbon, prediction should propose a candidate design, not certify it. Compliance with both the Cr(VI) and total-Cr limits should be established by sequential confirmatory testing with an error-controlled acceptance rule. Pure prediction was unsafe, and calibrated intervals were too wide to be useful, in the benchmark.
4. CASSIA's schema, evidence classes and extraction protocol are proposed as a basis for a shared, provenance-aware chromium-removal dataset. Their value for real data depends on the empirical application, which is the next step of this work.
