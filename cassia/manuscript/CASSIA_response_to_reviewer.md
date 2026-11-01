# Response to the Reviewer (third round)

**Manuscript:** CASSIA: A Provenance-Aware Computational Framework for Evidence-Constrained Analysis of Cr(VI) Removal by Carbonaceous Materials (revision 3; software CASSIA 0.3.1)

We thank the reviewer for the report. The page and line numbers it cites (for example, the fixed 50 mg L⁻¹ speciation lookup, the 40 Monte Carlo draws, the 15 virtual materials, the inner three-fold CV and Table 4 with 85.4% coverage) belong to the first version of the manuscript. Most of these points were addressed in the two previous revisions. Below, for each comment, we say what the current manuscript contains, where to find it, and what was added in this revision. Where a request is still not met, we say so.

**New in this revision.**
- A Monte Carlo precision analysis for the number of benchmark replicates, and the full true partition of Cr(VI) disappearance for all five scenarios (§2.5, §3.2, Tables S11–S12).
- An out-of-study validation of all four endpoint models used by the design layer, including the two total-Cr models (§3.6, Table S13).
- A table in §2.8 stating, for every interval method, the unit of the conformity score, the conformal procedure class and the coverage statement.
- Uncertainty propagation now also reports residual total Cr at accepted designs and explains the effect on fitted capacity (§3.7).
- The manuscript no longer refers to its own revision history; it reads as a stand-alone article.
- Code commit `0528c60cc3dcd567b246e9e883f41b5ef20de0ef`, regenerated SHA-256 manifest and test log.

---

## Major Comment 1 – Empirical validation

**Status: still not met.** We agree that this is the decisive weakness, and that the paper should contain both synthetic and empirical results.
- *The source.* We confirmed the reviewer's information: Ramirez et al. (2020) is open access on PubMed Central (PMC7044665), with raw isotherm and kinetic files as supplementary material. The companion article is J. Environ. Chem. Eng. 8, 103702. This is now stated in §3.8 and in the study registry (`data/study_registry.csv`).
- *Why it was not extracted.* On 1 November 2026 we tried again. The network policy of the computational environment rejected every host that serves the article or its files: www.ncbi.nlm.nih.gov, pmc.ncbi.nlm.nih.gov, doi.org, ScienceDirect, Crossref, Europe PMC, Mendeley Data and the institutional repository (repository.udem.edu.co). The same applies to the 20 further studies.
- *What we did not do.* We did not transcribe values from memory, estimate them, or digitise them from figures. Any of these would produce unverifiable numbers presented as empirical data.
- *What is in the manuscript.* §3.8 now gives the reviewer's seven-step plan as an ordered protocol: extraction, condition reconstruction, classification into C/B/A1/A2, identification of total-Cr and XPS/XANES studies, application of §3.3–3.5, and separate reporting of empirical findings. The pipeline that carries out these steps is tested and needs no code change. The Abstract, §3.8 and §4.7 state that no published data were analysed.

We recognise that the reviewer may recommend rejection on this point alone. The empirical application can be completed as soon as the source files are supplied by the author or the hosts are allowed.

## Major Comment 2 – Circular validation by the simulator

**Status: addressed (§2.4, §2.5, §3.2–3.6, Tables S1, S6, S11, S12; Fig. 5).**
- **Multi-mechanism benchmark.** The simulator implements the five scenarios the reviewer proposed:
  - A, adsorption only;
  - B, direct surface-mediated reduction;
  - C, adsorption-coupled reduction of *adsorbed* Cr(VI);
  - D, dissolved-phase reduction by leached reductants;
  - E, all pathways with Cr(III) re-adsorption, precipitation and pH drift.

  Chromium mass and electron equivalents are conserved, and unit tests check this for every scenario.
- **Biochar circularity.** Family-level reductant capacity, leachable reductant and specific capacity are redrawn from a common hyper-prior for every corpus, so no family is most reductive by construction. Even so, hydrochar had the highest reduction share in 83 of 100 scenario-E replicates. We traced this to the assumed proportionality between S_BET and adsorption capacity: a small adsorption sink inflates the reduction *share*. §3.2 reports this as a property of the simulator, not as evidence about hydrochars.
- **Replicates and precision.** There are 100 independent replicate corpora per scenario (15,000 virtual studies). The new precision analysis (Table S12) gives MCSE ≤ 0.010 for unconditional fit rates, ≤ 0.018 for conditional rates, ≤ 0.006 for mean reduction shares and ≤ 0.020 for correlations. The contrasts the conclusions rest on are 0.1–0.4. A ±0.01 half-width on the unconditional rates would need 228–357 replicates; no conclusion requires it.
- **True shares.** Table S11 gives median and 5–95% ranges of adsorbed Cr(VI), dissolved Cr(III), Cr(III) on the solid and precipitate for every scenario. For example, the median reduction share was 0.24 in B, 0.22 in C, 0.04 in D and 0.38 in E.
- **Robustness across mechanisms.** The central conclusions were re-tested under every scenario:
  - Isotherm fits overstated capacity in proportion to the reduction share in B, C and E, weakly in D (Table 3).
  - Random splits overstated transferability in all scenarios: R² 0.86–0.89 random against 0.48–0.65 LOSO (Fig. 5a, Table S6).
  - The study-level interval guarantee held in all scenarios (Fig. 5c).
- **Limits.** All scenarios belong to one kinetic model family (§4.7, limitation 2).

## Major Comment 3 – PHREEQC chemistry

**A. Fixed-concentration lookup. Status: addressed (§2.3, §3.1).**
- The lookup was replaced by a four-dimensional minteq.v4 grid of 7,700 PHREEQC solutions:
  - pH 1.5–10 in 0.25 steps;
  - dissolved Cr(VI) 0.01–700 mg L⁻¹;
  - I 10⁻³–10⁻¹ M;
  - T 288–318 K.
- The simulator evaluates speciation at the *current* dissolved Cr(VI), and the ML features at each row's C₀, pH, I and T.
- Against 300 direct PHREEQC solutions, the maximum interpolation error was 0.011 in species fraction.
- Major ions: NaCl, NaNO₃ and Na₂SO₄ backgrounds were compared (Table S10). Species fractions differed by ≤ 0.024. Ca, Mg, carbonate and DOM are stated as not explored.

**B. Cr(III) solubility. Status: addressed (§2.3, §3.1, Table S4).**
- Non-convergence is no longer interpreted. Saturation indices are computed directly, and saturation is located by bisection between 10⁻⁶ and 5,000 mg L⁻¹. Every evaluation carries a convergence status: bracketed, undersaturated to the maximum, supersaturated at the minimum, or speciation failed.
- Undersaturation up to 5,000 mg L⁻¹ extends to pH 4.25 (minteq.v4, I ≤ 0.01 M) and pH 3.5 (minteq).

**C. Database validation. Status: not met; stated as a limitation.**
- §3.1 now says explicitly that the 63–83 mV spread is a between-database difference, not a validated uncertainty interval.
- The tabulated constants of Ball and Nordstrom (1998) and Rai et al. (1987) could not be retrieved, for the same network reason as in Major Comment 1. CASSIA accepts a reference-constant table, so the comparison can be run as soon as the values are available.

**D. Uncertainty propagation. Status: addressed (§2.10, §3.7, Fig. 8).** Unreported ionic strength (10⁻³–10⁻¹ M) and pH basis (±0.5) are propagated with 1,000 draws per observation (previously 40) through the complete workflow:

| Quantity | Result |
|---|---|
| HCrO₄⁻ fraction | 90% range 0.035–0.24, largest at pH 5–6.5 |
| ML prediction | spread 0.08–0.14 log₁₀, against an interval width of 1.2–1.3 |
| Prediction interval | union interval 6–14% wider than nominal |
| Recommended dose | changed for 3 of 10 materials with accepted designs (300 draws each) |
| Residual total Cr at accepted designs | 95th percentile 0.08–1.38 mg L⁻¹ (limit 2); no non-compliant accepted design |
| Reduction share | bias up to 0.16 (90th percentile) for a 0.5-unit pH-basis error |
| Fitted capacity | q_max itself is unaffected, because isotherm fitting uses neither I nor pH basis; its interpretation changes via the reduction share |

## Major Comment 4 – Treatment-design methodology

**Status: addressed (§2.7, §2.9, §3.6, Tables 6, 7, 7b, S9, S13; Fig. 7).**
1. **Separate total-Cr models.** Residual Cr(VI) and residual total Cr are modelled directly as effluent concentrations. Each has its own RF model with the same feature set, trained only on rows with the corresponding measurement; for total Cr these are about 18 Class A/B studies per corpus.
2. **Validation (new, Table S13).** Leave-one-study-out R² (log scale) was 0.74 for C_e,VI and 0.65 for C_e,tot, against 0.94 and 0.91 under random splits.
3. **Uncertainty bounds for both endpoints.** Study-level intervals for all four endpoints are in Table 6. For total Cr, only about 17 training studies have total-Cr data. The CV+ intervals (M1, M2) are therefore effectively unbounded at 95%, because a study-level 95% quantile needs at least 19 calibration studies. The exact M2s interval, which calibrates on half of them (about 9), is unbounded at 95% and for 63% of test studies at 90%. This is reported as a data requirement.
4. **Simultaneous feasibility.** The joint criterion uses Bonferroni bounds (α/2 per endpoint) over the grid. Feasibility is evaluated against the simulator's true attainability for every material.
5. **Sensitivity.** pH ± 0.25/± 0.5 and dose × 0.8/× 1.2 were evaluated.
6. **Predicted against true effluent.** Spearman 0.65 for Cr(VI) and 0.55 for total Cr; the 95% upper bounds held for 97.9% and 100% of grid points.

**Scale of the evaluation.** The design analysis uses 20 independent cohorts of 15 materials (300 materials per target), each with its own training corpus. The first-version claim of a "seven-fold improvement" from 15 materials has been removed.

**Results.**
- Point predictions were unsafe for 20.7% (Wilson 16.5–25.6%) and 7.0% of materials for the two targets.
- Joint 90% bounds recommended almost nothing.
- Model-guided designs with sequential confirmatory testing were accepted for 98–105 materials per target, with no unsafe accepted design observed (Wilson upper bound 1.3%).
- The false-acceptance risk of the sequential verification was simulated separately (Table 7b). Under the error-controlled rule the worst per-material probability was ≤ 0.7%.

## Major Comment 5 – Study-aware conformal prediction

**Status: addressed (§2.8, §3.5, Tables 6, 6b, S5; Fig. 6; Appendix A).**
- **Unit of calibration.** A new table in §2.8 states, for each method:
  - M0: pooled observation-level scores from an inner grouped K-fold. This is a cross-validation residual method, not split or cross-conformal in the formal sense, and it has no study-level guarantee.
  - M1: observation-level scores weighted 1/n_j per study, in a study-level CV+ construction. No guarantee is claimed.
  - M2: one score per study, the within-study 0.9-quantile, in study-level CV+. No guarantee is claimed.
  - M2s: one score per study, in a split-conformal procedure over studies, with an exact guarantee proven in Appendix A: with probability ≥ 1 − α, at least ⌈0.9n⌉ of a new study's n observations are covered.
- **Requested statistics.** These are reported for all four endpoints and nominal levels 0.50–0.95 (10 corpora, outer leave-one-study-out):
  - coverage per study (Fig. 6b);
  - mean, median, 10th-percentile and worst-study coverage;
  - the fraction of studies covered;
  - mean width;
  - coverage against nominal level (Fig. 6a);
  - the number of independent calibration studies.
- **Results.**
  - Pooled M0 reached nominal *average* coverage but covered only 65–75% of studies at 90%, with worst studies at 11–29%.
  - M2s met its guarantee: 0.93 of studies at nominal 0.90. The price was 3.3 log₁₀ units of width, and the intervals were unbounded at 95%.
  - The guarantee is marginal, not family-conditional: hydrochar studies were covered at 84% (Table S5).
- **Earlier claims.** The first-version Table 4 (RF 85.4%, GPR 85.6%, conformal 90.3–94.9%) was based on pooled residuals and has been replaced.

---

## Release and reproducibility

| Item | Status |
|---|---|
| Version | CASSIA 0.3.1 |
| Code commit | `0528c60cc3dcd567b246e9e883f41b5ef20de0ef` |
| New scripts | `scripts/bench_partition.py` (Tables S11–S12), `scripts/endpoint_validation.py` (Table S13) |
| Tests | 39 tests pass; line coverage 81% |
| Manifest | `release/MANIFEST.sha256`, `release/RELEASE.md` |
| Tagged release and Zenodo DOI | **Author action** (cannot be created from the build environment) |

## Items that remain open

| Item | Reason | Remedy |
|---|---|---|
| Empirical application (Ramirez et al. + ≥ 20 studies) | Source hosts blocked by the environment's network policy | Author supplies the files, or the hosts are allowed |
| Validation against Ball and Nordstrom (1998) / Rai et al. (1987) | Same | Same; CASSIA accepts a reference-constant table |
| Release tag and DOI | Requires the author's GitHub/Zenodo account | Author action |
| Author name, affiliation, ORCID, CRediT, funding | Author-specific | Author to complete |
