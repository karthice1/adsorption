# CASSIA 0.3.0 release record

- Date: 2026-10-31
- Code commit at manifest generation (parent of the commit that adds this file): `9d5a37bd6012b6abf2cacbca284a97aada468f68`
- Content identity: `release/MANIFEST.sha256` (103 files); verify with `sha256sum -c release/MANIFEST.sha256`
- Environment: `requirements-lock.txt` (Python 3.13.17); PHREEQC databases with checksums and provenance in `databases/`
- Tests: `release/test_log.txt`; coverage: `results/coverage_report.txt`
- Permanent archive: to be created by the author (GitHub release `cassia-v0.3.0` + Zenodo DOI); not available from the build environment

## Reproduce

```
pip install -e ..                      # adsorption equation package (repository root)
pip install -e .[ml,pores,dev]
pytest -q tests                        # expected: all tests pass (see release/test_log.txt)
python scripts/run_v2.py --stage all   # v0.2 analyses (~4-5 h, 4 cores) -> results/v2/
python scripts/run_v3.py --stage all   # v0.3 analyses (~2-3 h, 4 cores) -> results/v3/, figures
```

## Expected output tables (60)

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
- `results/v3/T_censor_summary.csv`
- `results/v3/T_electrolyte_sensitivity.csv`
- `results/v3/T_equiv3_summary.csv`
- `results/v3/T_lc3.csv`
- `results/v3/T_m2s_family.csv`
- `results/v3/T_m2s_summary.csv`
- `results/v3/T_mlscen_summary.csv`
- `results/v3/T_mlscen_uq.csv`
- `results/v3/T_seqver_oc.csv`
- `results/v3/T_seqver_summary.csv`
- `results/v3/censor.json`
- `results/v3/censor_sensitivity.csv`
- `results/v3/chem3.json`
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
