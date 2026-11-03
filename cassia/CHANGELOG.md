# Changelog

## 0.3.7 (2026-11-03) - ninth revision
- Table S5c: the family-conditional cluster bootstrap now recomputes the pooled rate (the reported point estimate) in
  each resample, as Table S5b does; the estimand is defined in the caption, and the equal-corpus mean is written to the
  CSV for comparison (differs by at most 0.024).
- Reproducibility wording (Section 2.2, Section 2.12, Data availability) states which stages were rerun and that a
  single uninterrupted `reproduce_all.sh` run from the archive is still outstanding; archive naming gives the code
  and documents commits.
- Word tables: rows no longer split across pages.

## 0.3.6 (2026-11-03) - eighth revision
- Abstract and Conclusions: the M2s guarantee is stated as a finite-corpus randomisation guarantee for a randomly
  held-out study of the same realised corpus, not for a study from an independent literature population.
- `scripts/m2s_exchangeability.py --summarise` rebuilds Tables 7b/S5c from the saved per-study file and records the
  number of corpora behind each family-conditional bootstrap interval (`family_*_ncorpora`; 9 for hydrochar on the
  total-Cr endpoints); Table S5b/S5c captions say which corpora are resampled.
- Interval widths are stated on the log10(y + c) modelling scale (Section 3.5, Tables 7, 7b, S5b, S5d).
- Shortened for a computational-methodology journal: Appendix B -> Supplementary Note S1; full Sections 2.10/3.7 and
  former Fig. 8 -> Note S2 (Fig. S2); sequential-verification details and former Table 8b -> Note S3 (Table S15).

## 0.3.5 (2026-11-02) - seventh revision
- M2s: the reported evaluation is now the random role assignment of Appendix A (`scripts/m2s_exchangeability.py`,
  50 proper-training sets per corpus, all nominal levels; Table 7b, Table S5d, Fig. 6), with corpus-cluster bootstrap
  intervals for family-conditional rates (Table S5c). The original fixed-seed LOSO run is reported as a single realisation.
- Manuscript: guarantee stated for "a randomly selected held-out study from the same realised corpus"; extra-binomial
  type-I results marked nominal and exploratory; §4.4 items required per evidence class; §4.5 separates results from
  mathematics, the kinetic simulator and assumed measurement-error models; PHREEQC database agreement not presented as
  accuracy; electrolyte sensitivity limited to the electrolytes tested; Fig. 6 caption corrected.
- `requirements-lock.txt` records the interpreter used (CPython 3.13.18) with unchanged package versions.

## 0.3.4 (2026-11-02) - sixth revision
- M2s exchangeability under the hierarchical simulator (fixed family rotation, shared hyper-priors): Appendix A proves
  the required exchangeability from the random role assignment within a corpus and states the scope of the guarantee;
  `scripts/m2s_exchangeability.py` (Table S5c) enumerates role assignments on the 10 m2s corpora and reports
  within-corpus, family-conditional and new-corpus coverage; new unit test (42 tests).
- Manuscript: study-pooled versus corpus-level calibration named consistently; novelty claim narrowed; Highlights state
  that coverage is marginal; hydrochar Wilson interval corrected to 71-92%; test count 42 throughout.

## 0.3.3 (2026-11-02) - fifth revision
- Type-I error of the rejection test: new corpus-level analysis from the existing 360 controls, no refitting
  (`scripts/typeI_corpus.py`, Table S14b). Includes a corpus-level split-conformal threshold
  (`cassia.validation.cluster_conformal_threshold`), which is unbounded at 5% with 18 corpora, plus
  between-corpus heterogeneity and calibration-set sensitivity (explains the 6-corpus D1 rates of 14.2% and 10.8%).
  The study-pooled threshold is now reported as empirical average cross-corpus performance, without a finite-sample guarantee.
- Table S5b: corpus-cluster bootstrap intervals alongside the descriptive Wilson intervals; the reported
  interval is the wider envelope of the two (`scripts/m2s_family.py`).
- `scripts/reproduce_all.sh` runs without git metadata. The commit is taken from the argument, then
  `CASSIA_COMMIT`, then git, then `release/RELEASE.md`, and is resolved before any computation.
  `make_release_manifest.py` uses the same rule.
- README updated to 0.3.3 (reproduction, release checklist); new unit test for the cluster conformal threshold (41 tests).
- Fixed: `|` in supplementary table headers broke the Word table parser.

## 0.3.2 (2026-11-01) - fourth revision
- Rejection test: conformal threshold estimated across corpora; 240 new negative controls in 12 independent
  corpora; leave-one-corpus-out type-I error with corpus-bootstrap intervals (`scripts/typeI_validation.py`,
  Table S14). The previous 6-corpus calibration gave 14-15% type-I error for D1 on independent controls.
- Table S1 generated from the implementation (`scripts/make_table_s1.py`); measurement and evidence-class rows
  built from the code constants; unit test for agreement.
- M2s by material family with Wilson intervals (`scripts/m2s_family.py`, Table S5b).
- Verification denominators and both acceptance rules in one table (`scripts/verification_table.py`, Table 8c);
  the rule-comparison cohort documented as independent of the design cohort.
- `scripts/reproduce_all.sh` runs every script in dependency order; manuscript Table 2 lists outputs.
- Manuscript repositioned as a known-truth methodological benchmark.

## 0.3.1 (2026-11-01) - third revision
- Monte Carlo precision analysis of the mechanism benchmark and true partition of Cr(VI)
  disappearance for all scenarios (`scripts/bench_partition.py`; Tables S11-S12).
- Out-of-study validation of the four endpoint models used by the design layer, including the
  total-Cr models (`scripts/endpoint_validation.py`; Table S13).
- Operating characteristic of sequential verification as a stand-alone script (`scripts/seqver_oc.py`).
- Supplementary Fig. S1 (profile likelihoods; `scripts/make_supp_figs.py`).
- Manuscript: unit-of-calibration table for interval methods; revision-history language removed;
  empirical-application protocol expanded; Word build renders appendix mathematics.

## 0.3.0 (2026-10-31) - second revision
- Evidence classes C/B/A1/A2; separate measurement models for XPS (surface enrichment, beam/vacuum
  reduction), XANES (bulk), solid Cr inventory and phase evidence (`cassia.measurement`).
- Identifiability designs D1-D5 with technique nuisance parameters; precipitation share added.
- Adsorption-only test: calibrated rejection threshold (negative controls, split-half type-I), TOST-type
  noncentral-chi-square equivalence test with margins 5/10/20%, operating characteristics against
  the noise-free truth; three outcomes (rejected / equivalent / inconclusive).
- M2s: exact study-level split conformal with proof (Appendix A); M2 (CV+) guarantee withdrawn.
- ML and intervals for all five mechanism scenarios + cross-mechanism transfer; censoring/LOD
  sensitivity; learning curve with 8 replicates.
- Error-controlled sequential verification (upper confidence bound, Bonferroni over rounds and
  endpoints) with false-acceptance probabilities under shared and independent calibration bias.
- Background-electrolyte sensitivity (NaCl, NaNO3, Na2SO4).
- Release manifest (SHA-256 of all files, commit hash, commands, expected outputs).

## 0.2.0 (2026-10-30) - revision in response to peer review
- Chemistry: third Cr-bearing database (minteq.dat); speciation grid over pH x Cr(VI) x I x T with
  validation against direct PHREEQC; Cr(III) solubility by direct saturation-index bisection with
  convergence status; database provenance (SHA-256, $Id revision).
- Simulator: multi-mechanism model (scenarios A-E: adsorption, direct surface reduction,
  adsorption-coupled reduction, dissolved-phase reduction, combined with Cr(III) re-adsorption,
  precipitation and proton-balance pH drift); electron-equivalent and mass balances; six families with
  per-corpus hyper-priors; observation-level evidence classes; full parameter table.
- New identifiability engine: Fisher information/Laplace posterior by measurement design, model-free
  detection tests, adsorption-only equivalence tests with paired negative controls, profile likelihood,
  first-order metadata-bias propagation.
- Validation: >= 20 replicate corpora, cluster-bootstrap CIs, per-study metrics, nested grouped tuning,
  Duan smearing, learning curve, external corpora with different hyper-priors.
- Uncertainty: study-level CV+ conformal intervals (M1 weighted, M2 study-score) vs pooled residuals.
- Design: separate q_tot,aq model, joint (Bonferroni) two-endpoint constraint, 20 cohorts x 15 materials,
  operating-window sensitivity, predicted-vs-true effluent.
- Missing chemistry: 1000-draw propagation to predictions, intervals, dose and outcomes.
- Release files: LICENSE, CITATION.cff, requirements-lock.txt, databases/SHA256SUMS, coverage report.

## 0.1.0 (2026-10-30)
- Initial framework and manuscript draft (superseded; outputs archived in results/v0_1/).
