# Changelog

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
