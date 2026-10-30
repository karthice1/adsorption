# CASSIA

**C**hromium **A**dsorption, **S**peciation and **S**urface-reduction **I**dentifiability **A**nalysis:
a data-provenance-aware, speciation-constrained computational framework for telling
genuine Cr(VI) adsorption apart from apparent removal by reduction or precipitation on
carbonaceous adsorbents.

The manuscript is in [`manuscript/CASSIA_manuscript.md`](manuscript/CASSIA_manuscript.md).

## Layout

| Path | Content |
|---|---|
| `cassia/chemistry.py` | phreeqpython/PHREEQC Cr adapter: Cr(VI) speciation, redox scenarios, Cr(III) solubility, ML speciation features |
| `cassia/fitting.py` | Nonlinear isotherm/kinetic fits (equations from `adsorption.models`), AICc/BIC, Akaike weights, residual diagnostics, bootstrap CIs |
| `cassia/dataset.py` | Master-dataset schema, validation, endpoint definitions, evidence classes A/B/C |
| `cassia/simulate.py` | Known-truth virtual batch experiments (adsorption + reduction + Cr(III) partition) |
| `cassia/validation.py` | Random / grouped / leave-one-study-out / leave-one-family-out CV, ablation, conformal intervals, between-study variance |
| `cassia/pores.py` | pyGAPS adapter for published N2 isotherms (BET with Rouquerol criteria) |
| `databases/` | `minteq.v4.dat` and `llnl.dat` from the USGS PHREEQC distribution (the two Cr-bearing databases validated here) |
| `data/` | Empty master-dataset template, schema, study registry, extraction protocol |
| `scripts/run_all.py` | Regenerates every table, figure and number in the manuscript |
| `results/` | Outputs of `run_all.py` |

## Reproduce

```bash
pip install -e ..            # the adsorption package (repository root)
pip install -e .[ml,pores,dev]
pytest -q tests
python scripts/run_all.py    # full run; --quick for a smoke test
```

## What is and is not in this release

- All chemistry results come from PHREEQC with the two bundled databases.
- All model-benchmarking results come from a **simulated** corpus with known ground
  truth (`cassia/simulate.py`). They test whether the methods can recover the
  adsorption/reduction split. They are not empirical findings about real carbons.
- The published multi-study dataset (target: 20-30 studies) has **not** been extracted
  yet. `data/cassia_master_dataset.csv` is an empty template. The seed dataset
  (Ramirez et al., 2020, doi:10.1016/j.dib.2020.105292) is registered in
  `data/study_registry.csv` but not transcribed.
- RASPA3 is not used in this release (see the manuscript, Section 4.4).
