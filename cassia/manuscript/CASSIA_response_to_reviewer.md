# Response to the Reviewer

**Manuscript:** CASSIA: A Provenance-Aware Computational Framework for Evidence-Constrained Analysis of Cr(VI) Removal by Carbonaceous Materials
*(previous title: "CASSIA: A Data-Provenance-Aware, Speciation-Constrained Computational Framework for Distinguishing Cr(VI) Adsorption, Reduction and Precipitation in Carbonaceous Wastewater Treatment Systems")*

We thank the reviewer for an unusually thorough and constructive report. We agree with its central diagnosis: the original submission showed what CASSIA does in one synthetic environment of our own construction, and it did not show that the framework works on real experiments.

We have rebuilt the research rather than editing the text. The software is version 0.2.0. Every simulation, analysis and figure has been regenerated: 15,000 benchmark studies under five mechanisms, 350 identifiability studies, 20 replicate ML corpora, and 20 design cohorts of 15 new materials each.

**One requirement is not yet met: empirical validation (Major Comment 1).** We state this plainly here and in the manuscript (Sections 1, 3.8 and 4.7). The literature extraction requires access to publisher and repository hosts, which were not reachable from the computational environment in which this revision was produced. We have not substituted invented or mocked data. The manuscript is therefore now explicitly a methods paper, and we would understand if the editor prefers to see the empirical application before further consideration.

Reviewer text is in *italics*; our response follows each item. "§" refers to sections of the revised manuscript.

---

## Major Comment 1 – Absence of empirical validation

*"…the central claims … are based on artificial observations … Extract the Ramirez et al. dataset … at least 20 additional studies … run CASSIA against the actual observations … clearly separate empirical findings from synthetic benchmarking."*

**Response.** We agree that this is the most important gap, and it remains open in this revision. What we have done and not done:

1. **Not done.** The Ramirez et al. (2020) data (Data in Brief 29, 105292) and the further studies have not been extracted. The required hosts (pmc.ncbi.nlm.nih.gov, sciencedirect.com, zenodo.org, data.mendeley.com, figshare.com, the Crossref and Europe PMC APIs) refused connections from the build environment. We did not transcribe values from memory or fabricate them. The registry records this (`data/study_registry.csv`).
2. **Done: separation of evidence.** Synthetic and thermodynamic results are separated from any empirical claim:
   - every synthetic section and figure caption is labelled "synthetic, known-truth";
   - the Abstract, Highlights and Conclusions state that no empirical data were analysed;
   - a dedicated §3.8 reports the empirical status.
3. **Done: the empirical pipeline is ready to execute.** It consists of:
   - the schema (34 fields, row-level DOI and source location);
   - observation-level evidence classification;
   - QC checks (e.g. Ce,VI > Ce,tot, unit errors);
   - an end-to-end workflow script for user data (`examples/user_data_workflow.py`);
   - empty, dated execution records (`search_log.csv`, `prisma_flow.csv`).
   Running the Ramirez benchmark requires only the supplementary files.
4. **Changed: positioning.** The title, abstract and conclusions now present CASSIA as a framework for *evidence-constrained* interpretation (see §4.1 of the report). The empirical application is announced as the next step.

---

## Major Comment 2 – Circular validation by the simulator

*"…observed conclusions are dependent on the assumptions embedded in the simulator … build a multi-mechanism synthetic benchmark … at least 100 independent simulation replications … report the distributions of true adsorption, reduction and precipitation shares … The central conclusion should survive changes in the data-generating mechanism."*

**Response.** Implemented in full (§2.4–2.5, §3.2, Table S1, Fig. 2).

- **Five mechanisms.** One general kinetic model now contains five pathways, run as scenarios A–E:
  - A: adsorption only (negative control);
  - B: direct surface-mediated reduction;
  - C: adsorption-coupled reduction of *adsorbed* Cr(VI), which the reviewer identified as central;
  - D: dissolved-phase reduction by leached reductants followed by partitioning;
  - E: all pathways, plus Cr(III) re-adsorption, Cr(OH)₃(am) precipitation/dissolution and proton-balance pH drift.
- **Conservation.** Electron-equivalent and mass balances are exact and unit-tested for every scenario.
- **Replication.** 100 replicate corpora of 30 studies per scenario (15,000 studies). Monte Carlo SEs are reported (e.g. MCSE 0.8–1.0 percentage points for the R² ≥ 0.95 fractions).
- **Circularity.** The family-level means of reductant capacity, leachable reductant and specific capacity are redrawn from a common hyper-prior in every corpus. Remaining structure is reported transparently: hydrochars still had the highest reduction share in 83/100 replicates, which we traced to the assumed S_BET–capacity scaling (a small adsorption sink inflates the reduction share), not to a reductant assumption. §3.2 states that this must not be read as an empirical claim.
- **Distributions** of the true adsorbed, Cr(III)-on-solid, dissolved-Cr(III) and precipitated shares are reported for every scenario (Fig. 2a,d; `results/v2/bench.json`).
- **Robustness of the central conclusion.**
  - *Holds across mechanisms.* Among studies with > 20% reduction, isotherm fits still reached R² ≥ 0.95 in 46% (B), 41% (C), 17% (D) and 52% (E).
  - *Size depends on mechanism.* Publishable fits overstated q_max increasingly with reduction (Spearman ρ 0.60, 0.32, 0.13, 0.63). The effect is strong for B, C and E and weak for D, and §3.2 says so explicitly.
  - *Corrected statement.* We have replaced the original claim that the conclusions "do not depend on simulator details" with this mechanism-resolved statement.

---

## Major Comment 3 – PHREEQC chemistry

**A. Fixed-concentration speciation look-up.** *Corrected.*
- A PHREEQC grid over pH × dissolved Cr(VI) (0.01–700 mg L⁻¹) × I × T (7,700 solutions) is interpolated at the *current* dissolved Cr(VI), pH, I and T inside the simulator, and at each row's conditions for the ML features.
- Interpolation error against 300 direct PHREEQC calls: maximum 0.011, mean 0.0011 (§2.3, §3.1).
- Major-ion composition beyond the background electrolyte is not represented, because the schema records major ions but the benchmark does not vary them. This is stated as a limitation.

**B. Cr(III) solubility.** *Corrected.*
- Non-convergence is no longer interpreted. Saturation indices are evaluated directly, and solubility is located by bisection with an explicit status per evaluation.
- 174 evaluations: 146 bracketed, 28 undersaturated up to 5,000 mg L⁻¹, 0 failures. A 625-point SI map is fully converged.
- "No constraint below pH ≈ 4" is now a direct thermodynamic result (SI < 0 at 5,000 mg L⁻¹), and it is database dependent (pH ≤ 4.0 in minteq.v4, ≤ 3.5 in minteq).

**C. Database validation.** *Partially addressed.*
- A third Cr-bearing database (minteq.dat) was added, with checksums and revision provenance.
- The between-database differences are characterised (log K; speciation; redox boundary 63–83 mV; Cr(III) solubility differing by > 10× at pH 4).
- They are explicitly labelled as between-database differences, not validated uncertainty intervals.
- We could not retrieve the tabulated constants of Ball and Nordstrom (1998) or Rai et al. (1987) from the build environment, and we have not quoted values we could not verify. CASSIA accepts a reference-constant table for this comparison; §4.7 lists it as outstanding.

**D. Uncertainty propagation.** *Implemented through the complete workflow (§2.10, §3.7).* Using 1,000 Monte Carlo draws per condition (previously 40), the unreported ionic strength and pH basis are propagated to:
- the speciation feature;
- ML predictions and their M1 intervals;
- the recommended dose and the true residual total Cr (300 re-optimised designs per material);
- the inferred reduction share (first-order bias through the identifiability engine).

Results (§3.7):
- **Speciation feature.** Unknown ionic strength and pH basis widened the 90% range of f(HCrO₄⁻) to 0.24 at pH 5–6.5, but only to 0.035–0.08 elsewhere.
- **Predictions.** The spread of the point prediction was 0.08–0.14 log₁₀ units, and the union of the M1 intervals was only 6–14% wider than the nominal interval.
- **Design.** Under the test-confirmed design protocol, unknown ionic strength changed the accepted dose for 3 of the 10 feasible materials, and changed the probability of reaching an accepted design within four rounds (6–100% across draws for some materials). It never produced a non-compliant accepted design.
- **Mechanistic level.** A 0.5-unit pH error biased the inferred reduction share by up to 0.16 (90th percentile), against ≤ 0.025 for a ten-fold error in ionic strength.

Conclusion: missing chemistry is second-order for prediction and design and first-order (pH basis) for mechanism inference.

---

## Major Comment 4 – Treatment-design methodology

*"…how the second predictive model is trained, validated or calibrated … uncertainty bounds for both endpoints … simultaneous feasibility … sensitivity … predicted vs actual … repeated independent cohorts and a joint … safety criterion."*

**Response.** Rebuilt (§2.9, §3.6, Fig. 6).
1. **Separate models.** The design no longer infers effluent from uptake, which was ill-conditioned near complete removal. Residual Cr(VI) and residual total Cr are modelled directly, each by its own RF model with its own study-level CV+ interval. The total-Cr model is trained on the observations that carry a total-Cr analysis (Classes A/B).
2. **Validation.** Both endpoint models (and q_tot,aq) are validated by LOSO with per-study conformal coverage (Table 6).
3. **Joint criterion.** A Bonferroni rule (each endpoint at α = 0.05; joint miscoverage ≤ 0.10 by the union bound) is compared with point and marginal-90% rules.
4. **Replicated cohorts.** 20 independent cohorts × 15 new materials (300 per target, two target sets).
5. **Truth-based metrics.** True attainability is computed from the simulator for every material, giving unsafe recommendations, missed attainable designs and dose overhead, with Wilson intervals.
6. **Sensitivity and calibration.** Operating-window sensitivity (pH ± 0.25/0.5, dose × 0.8/1.2) and predicted-vs-true effluent are reported.

Results (§3.6, Table 7, Fig. 6; S1: C₀ = 20 mg L⁻¹, limits 0.5/2 mg L⁻¹; 169 of 300 materials truly attainable):

| Rule | Unsafe recommendations | Compliant designs found (of attainable) | Notes |
|---|---|---|---|
| Point predictions | 62/300 (20.7%; Wilson 16.5–25.6%) | 25% | Unsafe |
| Marginal and joint 90% (Bonferroni) | 0 | 0 | Recommended nothing |
| Joint 80% + 4 confirmatory batches | 2 | 7% | |
| Model-guided, test-confirmed | 0 | 62% | Mean 1.4 verification rounds; median dose twice the true minimum |

S2 gave the same picture.
- We therefore report the reviewer's concern as a finding: study-level uncertainty from 30 heterogeneous studies is too large to certify low effluent concentrations for a new material.
- Joint 95%-per-endpoint bounds for total Cr are unbounded with the ~17 studies that report total Cr; at least 19 are needed.
- The defensible engineering use of the models is to propose designs that are then confirmed by testing.

We agree that a marginal interval does not guarantee joint compliance after optimisation over a grid. The selected design is chosen *because* its bound is favourable, so §3.6 reports the realised compliance rather than the nominal one.

---

## Major Comment 5 – Study-aware conformal prediction

*"…the unit of calibration … whether the conformity score is calculated for each observation or for each study … report coverage per study, median, worst-study, width, coverage against nominal, number of calibration studies … clarify split, cross-conformal or other."*

**Response.** Clarified and replaced (§2.8, §3.5, Fig. 5, Table 6).
- **The v0.1 method is renamed M0** (pooled-residual) and described accurately. It is a cross-validation-calibrated heuristic without a formal guarantee for new studies.
- **Two cluster-level methods were added**, both CV+ constructions (Barber et al., 2021) with the *study* as the exchangeable unit:
  - **M1** weights each calibration study equally (observation weight 1/n_j; cf. Dunn et al., 2023) and targets coverage of a random observation of a new study;
  - **M2** uses one score per study (within-study 0.9-quantile of absolute residuals) and targets 'at least 90% of a new study's observations covered' with probability ≥ 1 − 2α.
- **Reporting.** For outer LOSO on 10 corpora (29 calibration studies per test study; four endpoints; nominal 0.50–0.95) we report per-study coverage distributions, the median, 10th-percentile and worst study, the fraction of studies covered at the nominal level, the width and the number of calibration studies.

Results (§3.5, Table 6, Fig. 5):
- **M0** reached nominal average coverage (0.89–0.90 at 90%), but only 65–75% of test studies were covered at ≥ 90%, and the worst study had 11–29% coverage.
- **M1** raised the share of adequately covered studies to 76–87%; **M2** raised it to 88–93%, above its 1 − 2α = 0.80 guarantee.
- **Width.** M2 intervals spanned 2.1–3.1 log₁₀ units.
- **Family-conditional coverage** was lower for hydrochars (0.77–0.91), so the guarantees are marginal over families, not conditional.

---

## Major Comment 6 – ML study design

- **A. LOFO with three families.** Six families are now used. LOFO now evaluates six held-out families per corpus, over 20 corpora. We agree this still does not represent the universe of carbonaceous materials, and say so.
- **B. Uncertainty of metrics.** Every metric now carries a 1,000-replicate study-level bootstrap interval, its distribution over 20 independent corpora, and paired model comparisons (RF − XGBoost per corpus, with a sign test). Result: RF exceeded XGBoost under LOSO by 0.038 ± 0.090 R², better in 14/20 corpora (sign test p = 0.12). The difference is not significant and is no longer presented as model superiority.
- **C. Model-development details.** Now given (§2.7, Table S3): fixed a-priori hyperparameters, seeds, scaling (GPR only), missing-data handling (I imputed at 0.01 M with an indicator) and fold-level per-study errors. Nested grouped tuning was evaluated on three corpora: it changed LOSO R² by +0.001 (RF) and −0.040 (XGBoost), so the fixed a-priori settings are retained.
- **D. Log target.** Negative apparent uptake is censored at q = 0 (flagged and counted), so q + 0.1 > 0 always. Values below the LOD enter as LOD/2. Back-transformation uses Duan's smearing estimator; the mean bias was −1.86 mg g⁻¹ with the naive inverse and −1.30 mg g⁻¹ with smearing. About 4.8% of observations were censored.
- **Effective sample size.** The learning curve (12–96 studies) shows how LOSO performance depends on the number of *studies*, which is the relevant sample size: LOSO R² was 0.47 ± 0.14 with 12 studies and 0.75 ± 0.08 with 96, while random-split R² stayed at 0.80–0.90. The between-study ICC was 0.45 ± 0.14.

---

## Major Comment 7 – Literature-extraction protocol

**Response.** We agree that a protocol is not an executed search, and we no longer describe the dataset as provenance-aware *data*; it is a provenance-aware *schema and protocol*.
- `data/PROTOCOL.md` now contains:
  - the exact search strings and databases;
  - dated execution logs (`search_log.csv`) and a PRISMA 2020 flow template (`prisma_flow.csv`), both empty because the search has not been run;
  - a duplicate-removal step and double-extraction QC with re-digitisation checks.
- Evidence classification is now **observation-level**, as recommended, with a study-level summary for reporting only.

---

## Major Comment 8 – Software contribution

| Requirement | Revision |
|---|---|
| Software architecture | Table 1, README (six layers, module responsibilities) |
| Installation instructions | README; `pyproject.toml` with extras `[ml,pores,dev]` |
| User-facing API | README "Minimal API"; `examples/user_data_workflow.py` |
| Versioned release | CASSIA 0.2.0, CHANGELOG |
| Licence | MIT (`LICENSE`) |
| Automated tests | 35 tests |
| Coverage report | 85% line coverage (`results/coverage_report.txt`) |
| Real-data workflow | Workflow script ready; demonstrated only on a labelled synthetic file, because no real data could be retrieved |
| Permanent DOI | **Pending.** We will mint it via the Zenodo–GitHub integration at acceptance or resubmission; `CITATION.cff` carries the placeholder |
| Reproducible environment | `requirements-lock.txt` (Python 3.13.17); database SHA-256 checksums and revisions (`databases/SHA256SUMS`, `PROVENANCE.md`) |

---

## Section-by-section comments

- **4.1 Title.** Changed to the reviewer's suggested wording.
- **4.2 Introduction.** Expanded on redox chemistry on carbon, thermodynamic data sources and group-aware uncertainty quantification. The contribution is restated as three questions (what can be established, what remains unidentifiable, how much confidence is warranted). We did not add citations we could not verify from the build environment; a comparison with existing chromium-data platforms is flagged for the author to complete.
- **4.3 Methods.** Every simulator parameter, distribution, consumption equation, partitioning law, measurement model and LSODA setting is in Table S1, so the simulator can be rebuilt without the source code.
- **4.4 Results.** All synthetic figures are labelled as such. Parameter distributions are referenced in every caption. Fig. 4 now carries across-corpus 95% intervals.
- **4.5 Discussion.** The claim that conclusions do not depend on simulator details is replaced by the mechanism-resolved evidence of §3.2. A new §4.2 separates thermodynamic feasibility from the kinetics of surface-mediated reduction.
- **4.6 Conclusions.** Divided into "Demonstrated computational findings" and "Proposed scientific and engineering implications".

---

## Minor comments

| # | Issue | Action |
|---|---|---|
| 1 | CASSIA not expanded | Defined: *Chromium Adsorption, Speciation and Surface-reduction Identifiability Analysis* |
| 2 | Author name incomplete | **Author to complete** (placeholder retained) |
| 3 | Affiliation missing | **Author to complete** |
| 4 | ORCID placeholder | **Author to complete** |
| 5 | CRediT placeholders | Roles listed; **author to confirm and add co-authors** |
| 6 | Ramirez et al. reference | Volume and article number added (29, 105292); full author list **to be verified by the author** |
| 7 | Software attribution | Versions given (phreeqpython 1.6.2, adsorption package by Amvro23 (MIT), pyGAPS 4.x, scikit-learn 1.9.1, XGBoost 3.4.1) |
| 8 | Database provenance | Revision lines, source URL, retrieval date and SHA-256 (Table S2) |
| 9 | q_total terminology | Renamed q_tot,aq ("total-aqueous-Cr uptake") and explicitly distinguished from adsorption |
| 10 | Negative uptake / LOD | Defined (§2.1, §2.7) |
| 11 | Tables | Tables 1 and 3 restructured; long tables moved to the Supplement |
| 12 | Figures | Fonts enlarged (≥ 9 pt), panel labels added, synthetic status in captions |
| 13 | Page layout | Footer line numbering suppressed and spacer paragraphs removed in the Word file |
| 14 | Funding | Statement added; **author to confirm** |
| 15 | AI declaration | Methods §2.12 describes the AI tool's role in code development and the safeguards. The declaration follows the current Elsevier template; the author must add a personal verification statement and check it against the target journal's current policy |

---

## The reviewer's five non-negotiable requirements

| # | Requirement | Status |
|---|---|---|
| 1 | Empirical literature validation | **Not met** (no data access); pipeline ready (Major Comment 1) |
| 2 | Multi-mechanism synthetic benchmark | **Met** (§3.2) |
| 3 | Correct statistical validation | **Met** for cluster-level conformal, total-Cr prediction and joint design (§3.4–3.6) |
| 4 | Strengthened PHREEQC validation | **Partly met**: concentration dependence and Cr(III) saturation met; independent constant validation outstanding (Major Comment 3C) |
| 5 | Publishable software release | **Met except the DOI**, which requires the author's Zenodo account |
