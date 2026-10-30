# Changelog

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
