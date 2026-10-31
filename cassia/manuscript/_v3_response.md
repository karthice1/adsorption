# Response to the Reviewer (second round)

**Manuscript:** CASSIA: A Provenance-Aware Computational Framework for Evidence-Constrained Analysis of Cr(VI) Removal by Carbonaceous Materials (revision 2; software CASSIA 0.3.0)

We thank the reviewer for a second, detailed report. Every major comment led to new analyses, and the manuscript text was revised to match their results. Where a request could not be met, we say so plainly and explain why. Section, table and figure numbers refer to the revised manuscript; "S" numbers refer to the Supplementary Material.

**Summary of what changed.**
- A formal equivalence test (noncentral χ², with margins of 5/10/20%) was added. The rejection test is now calibrated, and its type-I error and power are reported with cross-validation (§2.6, §3.3, Table 4, Table S8, Appendix B).
- Evidence classes are now C, B, A1 and A2. XPS and XANES are simulated with separate, explicitly stated bias models, and solid Cr inventory and phase evidence are simulated (§2.1, §2.4, §2.6).
- An exact study-level split-conformal interval (M2s) is given with a theorem and proof (Appendix A). The earlier guarantee for M2 is withdrawn.
- The false-acceptance risk of sequential verification was simulated for two acceptance rules (Table 7b, Table S9).
- The ML conclusions were re-tested across mechanism scenarios A–E (Table S6, Fig. 5). The learning curve was repeated with 8 replicates, and censoring sensitivity was added (Table S7).
- The PHREEQC pH-threshold discrepancy was corrected, and an electrolyte sensitivity analysis was added (Table S10).
- Release: code commit `{COMMIT}`, SHA-256 manifest, test log and coverage report (`release/`).
- **Not met:** the empirical application (Major Comment 1). This is explained below.

---

## Major Comment 1 – Empirical validation with Ramirez et al. (2020) and 20–30 studies

**Response: not met in this revision.** We could not retrieve the Data in Brief article of Ramirez et al. (2020; 29, 105292), its Mendeley Data supplement, or the primary articles of a further 20–30 studies.
- *Why.* The computational environment used for this revision allows outbound access only to an approved list of hosts. Every publisher and repository host needed returned HTTP 403: ScienceDirect, PubMed Central, Europe PMC, Crossref, Mendeley Data, Zenodo and Figshare.
- *What we did not do.* We did not transcribe values from memory and did not estimate them, because that would produce unverifiable data presented as empirical.
- *What is ready.* The extraction protocol (`data/PROTOCOL.md`, now with the C/B/A1/A2 criteria), the schema, the study registry (with the full Ramirez et al. author list and article number), the class-assignment code and the user-data workflow are complete and tested. Loading the Ramirez data requires no code change.
- *What the manuscript says.* §3.8 and §4.7 (limitation 1) state that the requirement remains unmet. Every quantitative claim is labelled synthetic or thermodynamic. The Abstract ends with the same statement.

We recognise that the reviewer may consider this decisive. If the editor agrees, the empirical application can be completed as soon as the source files are available, either supplied by the author or retrieved from a network-enabled environment.

## Major Comment 2 – A formal equivalence test instead of "could not be rejected"

**Response: done (§2.6, §3.3, Table 4, Table S8, Appendix B).**
- **Margins.** The margin is a systematic discrepancy of relative size ε in every concentration, with ε = 5, 10 and 20%, plus 0.05 in the oxidation-state fraction. It is expressed as a non-centrality λ₀(ε) = Σ(m_i/σ_i)².
- **Test.** Equivalence is declared when χ² ≤ the α-quantile of the noncentral χ²(ν, λ₀), with α = 0.05. This is the goodness-of-fit form of an equivalence test (Wellek, 2010), and, as in TOST (Schuirmann, 1987), equivalence must be demonstrated.
- **Calibrated type-I error.** The nominal χ² test over-rejected (6.7–15%), because the empirical model only approximates the generating model. Rejection thresholds are now the 95th percentiles of 120 paired negative controls. Their split-half cross-validated type-I error is 3.3–12.5%, mostly 3–7%.
- **Power.** For high-reduction studies, calibrated power was 24–37% for Class C designs, 80–88% with total Cr and 90–96% with solid-phase evidence.
- **Operating characteristics.** The simulator provides the noise-free misfit of every study, so the false-equivalence rate (0–4.3%, ≤ α) and the equivalence power (96–99% at 3% noise) are reported against the truth.
- **GOF details.** Appendix B states the observation vector, the error model, the treatment of LOD and shared calibration error, the parameters of the adsorption-only model, the fitting procedure, the degrees of freedom and the decision rule.
- **Terminology.** The manuscript now distinguishes three outcomes: "rejected by the calibrated test", "equivalent within ε" and "inconclusive". The phrase "adsorption-only was not rejected by the specified statistical test" is used wherever only non-rejection is shown. "Not rejected" is never read as structural equivalence.

**What the new test showed.** For high-reduction studies observed through residual Cr(VI) only, adsorption-only was *demonstrated equivalent within 10%* for 43–65% of studies. For 84–92% of these high-reduction studies the noise-free misfit was truly within that margin, so the test recovered a genuine structural equivalence. The earlier statement therefore holds in a stronger, formally tested form, for those studies.

## Major Comment 3 – Evidence classes and solid-phase measurements

**Response: done (§2.1, §2.4, §2.6, §3.3, §4.1, §4.4; `cassia/measurement.py`, `cassia/dataset.py`).**
- **Classes.**
  - C: residual Cr(VI) only.
  - B: adds total aqueous Cr.
  - A1: adds the solid oxidation state, the pH basis and the ionic strength.
  - A2: adds a solid Cr inventory and phase evidence as well.
  The class is assigned automatically per observation. The protocol table in `data/PROTOCOL.md` was updated to match.
- **XPS and XANES are modelled separately.**
  - XPS: log-normal surface enrichment of Cr(III) (SD 0.4 in log), beam/vacuum reduction of Cr(VI) uniform on 0–15%, and noise SD 0.05.
  - XANES: beam reduction 0–3% and noise SD 0.03.
  - Solid inventory: 5% relative error. Phase fraction: 0.10 absolute error.
- **Designs.** D4x (A1 by XPS), D4n (A1 by XANES) and D5 (A2) are evaluated separately. The empirical adsorption-only model may explain an apparent solid Cr(III) signal by beam reduction, up to the technique's assumed maximum.
- **Results.**
  - A simple XPS threshold false-alarmed in 7.1% of pure-adsorption controls, against 1.0% for XANES.
  - In the fitted test, XPS and XANES gave similar power (90–96%).
  - Only A2 constrained the precipitated share (posterior SD 0.026, against 0.037 for D4n).
  - All detectability statements are explicitly conditional on these assumed measurement models (§3.3, §4.7 limitation 4).

## Major Comment 4 – False-acceptance control in sequential verification; conclusions aligned with results

**Response: done (§2.9, §3.6, Table 7b, Table S9, §5).**
- **Simulation.** 1,000 verification sequences were simulated per material and target, over 300 materials per target. Each sequence had up to four dose rounds and both endpoints, a shared or an independent calibration bias, and correlated replicate errors.
- **Rules compared.** The v0.2 rule (duplicates, mean ≤ 0.9·limit) was compared with an error-controlled rule: triplicates, with a one-sided upper confidence bound ≤ limit, Bonferroni-corrected over two endpoints and four rounds.
- **False-acceptance probability.** The worst per-material value was 1.4–1.7% for the v0.2 rule and 0.1–0.7% for the controlled rule. No material exceeded 5% under either rule.
- **Operating characteristic** (`scripts/seqver_oc.py`). At a true concentration equal to the limit, the acceptance probability over four rounds was 3.7% for the v0.2 rule and 0.8% for the controlled rule.
- **Wilson intervals.** All design rates now carry Wilson 95% intervals. For the accepted designs, "no unsafe recommendation" is replaced by "no unsafe accepted designs were observed in the simulated test cohorts", with a Wilson 95% upper bound of 1.3% of materials per target.
- **Conclusions** (§5.1 item 6, §5.2 item 3). The conclusions now say what the results support:
  - the model proposes a candidate design;
  - compliance is established by sequential confirmatory testing;
  - pure prediction was unsafe;
  - calibrated intervals were too wide to be useful.

## Major Comment 5 – Exact theorem and proof for the study-level conformal method

**Response: done, with a correction (§2.8, §3.5, Table 6b, Appendix A).**
- **Correction.** We could not prove the 1 − 2α guarantee claimed in v0.2 for M2 (CV+ with within-study quantile scores). The claim is withdrawn, and M2 is reported as an empirical method only.
- **New method M2s.** We added an exact study-level split-conformal method, M2s. Appendix A states the setting, the assumption (exchangeable studies conditional on the proper-training set), the theorem and the proof. With probability ≥ 1 − α, at least ⌈βn⌉ of a new study's n observations are covered, and the probability is ≤ 1 − α + 1/(J+1) for distinct scores.
- **Remarks in Appendix A.** These cover marginal versus family-conditional validity, the conservatism for small studies, the minimum number of calibration studies (J ≥ (1−α)/α) and the relation to Dunn et al. (2023).
- **Empirical check.** M2s tracked its nominal level (0.93, 0.80, 0.72 and 0.50 at nominal 0.90, 0.80, 0.70 and 0.50), within the theoretical bounds, in all five scenarios. It is unbounded at 95% with literature-scale numbers of studies, and this is reported as a result.

## Major Comment 6 – ML across scenarios, GPR, learning curve, censoring

**Response: done (§3.4, §4.3, Table S6, Table S7, Fig. 5).**
- **Scenarios A–E.** The ML conclusions were re-tested on five corpora per mechanism scenario.
  - Random-split R² was 0.86–0.89, LOSO 0.48–0.65 and LOFO 0.40–0.63, so the random-split inflation is not specific to scenario E.
  - Speciation features lowered LOSO R² in every scenario.
  - Models trained on scenario E and transferred to scenarios B and D lost 0.05–0.10 in R².
  - These claims are restricted to the simulator's model family.
- **GPR.** It is labelled exploratory (3 corpora, untuned kernel), and no general conclusion is drawn from it (§4.3 item 5).
- **Learning curve.** It was repeated with 8 replicates per corpus size; SE and the Spearman trend are reported. With 8 replicates, the mean LOSO R² was 0.31, 0.58, 0.66 and 0.66 for 12, 24, 48 and 96 studies (SE 0.09–0.02; Spearman ρ = 0.58, p < 0.001). The curve is monotonic in its means and levels off beyond about 48 studies. The non-monotonic v0.2 curve (4 replicates) is explained in the text as replicate noise.
- **Censoring.** The main analysis censors negative apparent uptake (about 5% of observations) at zero. The alternatives were LOD-based censoring, LOD/2, dropping, and an asinh transform without censoring. They gave LOSO R² 0.58–0.68 and mean M1 coverage 0.94 in every case (Table S7), so the conclusions are unchanged.

## Major Comment 7 – PHREEQC validation, pH threshold, background electrolyte

**Response: partly done (§3.1, §4.7 limitation 3, Table S4, Table S10).**
- **pH threshold corrected.** The text now matches Table S4. On the 0.25-pH grid, Cr(III) is undersaturated up to 5,000 mg L⁻¹ to pH 4.25 in minteq.v4 at I = 0.001–0.01 M (4.5 at I = 0.1 M), and to pH 3.5 in minteq.
- **Electrolyte.** NaCl, NaNO₃ and Na₂SO₄ backgrounds were compared. Cr(VI) species fractions differed by ≤ 0.024, the HCrO₄⁻/CrO₄²⁻ crossover by ≤ 0.004 pH units, and the CrO₃SO₄²⁻ ion pair stayed ≤ 1%. Ca, Mg, carbonate and DOM were not explored, and this is stated.
- **Independent validation: outstanding.** The critically evaluated constants (Ball and Nordstrom, 1998; Rai et al., 1987) could not be retrieved and verified from the environment, for the same reason as in Major Comment 1. The manuscript states this as a limitation and does not claim the validation.

---

## Software, release and reproducibility

| Item | Status |
|---|---|
| Version | CASSIA 0.3.0 (`pyproject.toml`, `CITATION.cff`, CHANGELOG) |
| Code commit | `{COMMIT}` (branch `claude/exciting-hawking-ayadsw`, repository `karthice1/adsorption`) |
| Content manifest | `release/MANIFEST.sha256` (SHA-256 of all code, data, databases, results and figures); `release/RELEASE.md` |
| Run commands | `scripts/run_v2.py --stage all`; `scripts/run_v3.py --stage all`; `scripts/seqver_oc.py`; `scripts/make_supp_figs.py` |
| Logs | `release/test_log.txt`; stage logs in `results/` |
| Tests | {NTESTS} tests pass; line coverage {COV} (`results/coverage_report.txt`) |
| Tagged release and Zenodo DOI | **Author action.** The tag `cassia-v0.3.0` and the Zenodo archive must be created from the GitHub account; this could not be done from the build environment. The commit and manifest identify the content until then |

## Status of the response document items (revised)

| Reviewer request | Status |
|---|---|
| Empirical validation (Ramirez + 20–30 studies) | **Not met** – source hosts unreachable; pipeline ready |
| Formal equivalence test, margins, calibrated type-I error and power, GOF details | Done |
| Terminology "not rejected by the specified statistical test" | Done |
| Evidence classes C/B/A1/A2; XPS/XANES modelled separately; solid inventory | Done |
| Sequential verification false-acceptance control | Done (simulated; assumed error model) |
| Wilson CIs; qualified "no unsafe accepted designs observed" | Done |
| Exact theorem and proof for the study-level conformal method | Done for M2s; M2 guarantee withdrawn |
| ML across scenarios A–E or restricted claim | Done (and claims restricted to the model family) |
| GPR exploratory; learning curve qualified; censoring sensitivity | Done |
| PHREEQC pH threshold vs Table S4 | Corrected |
| Independent PHREEQC validation | **Outstanding** – constants not retrievable |
| NaCl limitation | Addressed for NaCl/NaNO₃/Na₂SO₄; other ions stated as limitation |
| Tagged release, DOI | Commit and manifest done; **tag/DOI: author action** |
| Highlights qualified | Done |
| Placeholders removed | Done, except author-specific items (below) |
| Dunn et al. reference | Corrected: J. Am. Stat. Assoc. 118(544), 2491–2502 |
| Last-page page break | Fixed in the Word build |

## Highlights

The highlights were rewritten so that every quantitative claim is qualified as simulated:
{HIGHLIGHTS}

## Minor comments

| # | Issue | Action |
|---|---|---|
| 1 | Placeholders | `[COMMIT]` and `[DOI]` placeholders removed. Commit inserted; the DOI sentence now states that the tag and archive will be created by the author |
| 2 | Author name, affiliation, ORCID, CRediT, funding | **Author to complete.** These cannot be supplied by us |
| 3 | AI-use declaration | Retained; **author to confirm** against the journal's current policy |
| 4 | Dunn et al. (2023) | Volume, issue and pages corrected |
| 5 | Ramirez et al. (2020) | Full author list and article number |
| 6 | Software versions | Moved to Data and code availability |
| 7 | Profile likelihood figure | Moved to Fig. S1; non-smooth profile points disclosed as local optima and the profiles described as indicative only |
| 8 | Final page break | Removed from the Word build |
