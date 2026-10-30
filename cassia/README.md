# CASSIA 0.2

**CASSIA** (*Chromium Adsorption, Speciation and Surface-reduction Identifiability Analysis*) is an
open Python framework for **evidence-constrained** analysis of Cr(VI) removal by carbonaceous
materials. It keeps apparent Cr(VI) disappearance, total-aqueous-Cr disappearance and adsorption
apart, tells you which removal mechanisms a given set of measurements can and cannot identify, and
gives study-level calibrated prediction intervals and treatment designs.

The revised manuscript is `manuscript/CASSIA_manuscript.md` (Word version `CASSIA_manuscript.docx`);
the point-by-point reply to the reviewer is `manuscript/CASSIA_response_to_reviewer.md` (`.docx`), and the supplement `manuscript/CASSIA_supplementary.md` (`.docx`).

## Architecture

| Layer | Module | What it does |
|---|---|---|
| 1 Data | `cassia.dataset` | 34-field schema with DOI/source provenance, validation and QC flags, endpoints (`q_app_VI`, `q_tot_aq`), **observation-level** evidence classes A/B/C |
| 2 Chemistry | `cassia.chemistry` | phreeqpython/PHREEQC adapter for three Cr-bearing databases: Cr(VI) speciation, `SpeciationGrid` over pH x Cr(VI) x I x T, redox scenarios, Cr(III) saturation by direct SI bisection, database provenance |
| 3 Simulator | `cassia.simulate` | Known-truth multi-mechanism batch model (scenarios A-E), mass and electron balances, six material families with per-corpus hyper-priors, measurement model |
| 4 Identifiability | `cassia.identifiability` | Fisher information / Laplace posterior per measurement design, model-free detection tests, adsorption-only equivalence test, profile likelihood, metadata-bias propagation |
| 5A Models | `cassia.fitting`, `cassia.validation` | Nonlinear isotherm/kinetic selection (AICc, diagnostics, bootstrap); RF/XGBoost/GPR with random, grouped, leave-one-study-out and leave-one-family-out validation, cluster bootstrap CIs |
| 5B Uncertainty | `cassia.validation.StudyCVPlus` | Study-level CV+ conformal intervals (M1, M2) and the pooled-residual baseline (M0) |
| 6 Design | `scripts/run_v2.py` (`stage_design`) | Two-endpoint (Cr(VI), total Cr) dose/pH design with point, marginal and joint rules |
| aux | `cassia.pores` | pyGAPS adapter: BET (Rouquerol criteria) from published N2 isotherms |

## Install

```bash
git clone https://github.com/karthice1/adsorption && cd adsorption
git checkout claude/exciting-hawking-ayadsw
pip install -e .                      # the 'adsorption' equation package (repository root)
pip install -e cassia[ml,pores,dev]   # CASSIA, xgboost, pyGAPS, pytest
pytest -q cassia/tests                # 35 tests
```

Python 3.10+; the exact environment used for the paper is pinned in `requirements-lock.txt`
(Python 3.13.17). PHREEQC databases ship in `databases/` with checksums (`sha256sum -c databases/SHA256SUMS`).

## Minimal API

```python
from cassia import chemistry as ch, dataset as ds, fitting as fit, validation as val

ch.speciate_cr6(pH=3.0, cr_mgL=100, database="minteq.v4.dat").fractions   # Cr(VI) species per Cr atom
ch.cr3_solubility_si(pH=6.0)            # {'S_mgL': 0.0225, 'status': 'bracketed', ...}
df = ds.classify(ds.endpoints(df))      # endpoints + observation-level evidence classes
df = ch.speciation_features(df)         # PHREEQC features at each row's pH, C0, I, T
fits = fit.fit_isotherms(Ce, q, T=298.15); fit.selection_table(fits)
cp = val.StudyCVPlus("RF", val.CHEM).fit(train_df, "q_app_VI")
lo, hi, centre = cp.interval(new_df, alpha=0.1, method="M1")   # log10 scale; val.from_target() to mg/g
```

A complete workflow for user data in the schema is `examples/user_data_workflow.py`
(`--demo` runs it on a clearly labelled synthetic file).

## Reproduce the paper

```bash
cd cassia
python scripts/run_v2.py --stage all      # ~2 h on 4 cores; --quick for a smoke test
```

Stages `chem, bench, ident, ml, uq, design, missing, figures` can be run separately. Outputs go to
`results/v2/` (tables `T_*.csv`, per-stage JSON summaries) and `manuscript/figures/`.

## Scope and status

- Results in the paper are (i) PHREEQC thermodynamic calculations and (ii) known-truth simulation
  benchmarks. **No empirical literature data have been processed yet**: the publisher and data-repository
  hosts were not reachable from the build environment. `data/cassia_master_dataset.csv` is the empty
  template, `data/study_registry.csv` lists the seed study (Ramirez et al., 2020), and
  `data/PROTOCOL.md` the registered extraction protocol.
- RASPA3 is not used (manuscript Section 4.6).

## Release checklist (for the author)

1. Fill in author name/ORCID in `CITATION.cff` and `LICENSE`.
2. Tag the release commit `cassia-v0.2.0` and create a GitHub release.
3. Enable the repository in Zenodo (GitHub integration) so the release receives a DOI; put the DOI in
   `CITATION.cff` and in the manuscript's Data availability section.
