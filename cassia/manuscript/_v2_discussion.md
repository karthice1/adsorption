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

[ML DISCUSSION]

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
4. [ML CONCLUSION]
5. [UQ CONCLUSION]
6. [DESIGN CONCLUSION]

### 5.2 Proposed scientific and engineering implications

1. Residual-Cr(VI) isotherm constants should be reported, and pooled, as apparent removal capacities. Mechanistic claims require Class A evidence (§4.4).
2. ML models of Cr(VI) removal should be validated by leave-one-study-out (or family-out) splits, with study-level uncertainty, before they are used for new materials.
3. Treatment designs for a new carbon should be based on calibrated upper bounds for both residual Cr(VI) and total Cr, and should be confirmed experimentally.
4. CASSIA's schema, evidence classes and extraction protocol are proposed as a basis for a shared, provenance-aware chromium-removal dataset. The empirical application will be reported separately.
