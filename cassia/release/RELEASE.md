# CASSIA 0.3.2 release record

- Date: 2026-11-01
- Code commit at manifest generation (parent of the commit that adds this file): `72a0b175c07b676093673cbf2b1e8c976de87142`
- Content identity: `release/MANIFEST.sha256` (120 files); verify with `sha256sum -c release/MANIFEST.sha256`
- Environment: `requirements-lock.txt` (Python 3.13.17); PHREEQC databases with checksums and provenance in `databases/`
- Tests: `release/test_log.txt`; coverage: `results/coverage_report.txt`
- Permanent archive: to be created by the author (GitHub release `cassia-v0.3.2` + Zenodo DOI); not available from the build environment

## Reproduce

```
pip install -e ..                      # adsorption equation package (repository root)
pip install -e .[ml,pores,dev]
pytest -q tests                        # expected: 40 tests pass (see release/test_log.txt)
bash scripts/reproduce_all.sh          # every table, figure and number, in dependency order (~9-11 h, 4 cores)
```

## Expected output tables (70)

- `results/v2/T_S1_simulator_parameters.csv`
- `results/v2/T_S3_hyperparameters.csv`
- `results/v2/T_cr3_solubility_SI.csv`
- `results/v2/T_db_provenance.csv`
- `results/v2/T_logK.csv`
- `results/v2/T_ml_summary.csv`
- `results/v2/T_redox_boundary.csv`
- `results/v2/T_uq_summary.csv`
- `results/v2/bench.json`
- `results/v2/bench_studies.csv.gz`
- `results/v2/chem.json`
- `results/v2/design.json`
- `results/v2/design_outcomes.csv`
- `results/v2/design_pred_vs_true.csv`
- `results/v2/design_sensitivity.csv`
- `results/v2/ident.json`
- `results/v2/ident_equivalence.csv`
- `results/v2/ident_fim_detection.csv`
- `results/v2/ident_profiles.csv`
- `results/v2/missing.json`
- `results/v2/missing_design.csv`
- `results/v2/missing_prediction.csv`
- `results/v2/ml.json`
- `results/v2/ml_learning_curve.csv`
- `results/v2/ml_metrics.csv`
- `results/v2/ml_per_study_LOSO_RF.csv`
- `results/v2/ml_tuning.csv`
- `results/v2/speciation_C0_10.csv`
- `results/v2/speciation_C0_500.csv`
- `results/v2/timing.json`
- `results/v2/uq.json`
- `results/v2/uq_per_study.csv`
- `results/v3/T_S1_simulator_parameters.csv`
- `results/v3/T_bench_mcse.csv`
- `results/v3/T_bench_partition.csv`
- `results/v3/T_censor_summary.csv`
- `results/v3/T_electrolyte_sensitivity.csv`
- `results/v3/T_endpoint_validation.csv`
- `results/v3/T_equiv3_summary.csv`
- `results/v3/T_lc3.csv`
- `results/v3/T_m2s_family.csv`
- `results/v3/T_m2s_family_full.csv`
- `results/v3/T_m2s_summary.csv`
- `results/v3/T_mlscen_summary.csv`
- `results/v3/T_mlscen_uq.csv`
- `results/v3/T_seqver_oc.csv`
- `results/v3/T_seqver_summary.csv`
- `results/v3/T_typeI_validation.csv`
- `results/v3/T_verification_denominators.csv`
- `results/v3/censor.json`
- `results/v3/censor_sensitivity.csv`
- `results/v3/chem3.json`
- `results/v3/endpoint_validation.csv`
- `results/v3/equiv3.csv.gz`
- `results/v3/equiv3.json`
- `results/v3/equiv3_part1.csv.gz`
- `results/v3/equiv3_part2.csv.gz`
- `results/v3/ident3.json`
- `results/v3/ident3_fim_detection.csv`
- `results/v3/lc3.csv`
- `results/v3/lc3.json`
- `results/v3/m2s.json`
- `results/v3/m2s_per_study.csv`
- `results/v3/mlscen.json`
- `results/v3/mlscen_metrics.csv`
- `results/v3/mlscen_uq.csv`
- `results/v3/seqver.csv`
- `results/v3/seqver.json`
- `results/v3/typeI_ctrl_part1.csv.gz`
- `results/v3/typeI_ctrl_part2.csv.gz`
