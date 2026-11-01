#!/usr/bin/env bash
# Regenerate every table, figure and number of the CASSIA 0.3.x paper, in dependency order.
# Run from the cassia/ directory after `pip install -e ..` and `pip install -e .[ml,pores,dev]`.
# Total time about 9-11 h on 4 cores; each line can also be run on its own.
set -euo pipefail
export PYTHONPATH=.
pytest -q tests | tee release/test_log.txt                     # 40 tests
python scripts/run_v2.py --stage all                           # results/v2: chemistry, benchmark, ML, UQ, design, missing chemistry; Figs 1, 2, 4, 8
python scripts/run_v3.py --stage chem3                         # electrolyte sensitivity (Table S10), thresholds (Table S4)
python scripts/run_v3.py --stage ident3                        # Fisher analysis, model-free diagnostics (Sect. 3.3)
EQUIV_PART=1 python scripts/run_v3.py --stage equiv3           # equivalence test, scenario E (~45 min)
EQUIV_PART=2 python scripts/run_v3.py --stage equiv3           # equivalence test, scenarios B-D (~45 min)
python scripts/run_v3.py --stage equiv3                        # summary: Table 5, Table S8
python scripts/run_v3.py --stage m2s                           # M2s (Table 7b)
python scripts/run_v3.py --stage mlscen                        # ML across scenarios (Table S6)
python scripts/run_v3.py --stage censor                        # censoring sensitivity (Table S7)
python scripts/run_v3.py --stage lc3                           # learning curve (Fig. 5b)
python scripts/run_v3.py --stage seqver                        # sequential verification (Table 8b)
python scripts/run_v3.py --stage figures                       # Figs 3, 5, 6, 7 (overwrites the run_v2 versions)
CTRL_PART=1 python scripts/typeI_validation.py                 # independent negative controls, part 1 (~45 min)
CTRL_PART=2 python scripts/typeI_validation.py                 # part 2 (~45 min)
python scripts/typeI_validation.py                             # Table S14: type-I error with Wilson intervals
python scripts/seqver_oc.py                                    # Table S9
python scripts/bench_partition.py                              # Tables S11-S12
python scripts/endpoint_validation.py                          # Table S13 (~15 min)
python scripts/m2s_family.py                                   # Table S5b
python scripts/verification_table.py                           # Table 8c
python scripts/make_table_s1.py                                # Table S1 from the implementation
python scripts/make_supp_figs.py                               # Fig. S1
(cd manuscript && python make_supplement_v3.py && python assemble_v3.py "$(git rev-parse HEAD)")
python scripts/make_release_manifest.py "$(python -c 'import cassia; print(cassia.__version__)')"
