# CASSIA systematic data-extraction protocol (PRISMA 2020-aligned)

Status: **protocol registered; extraction not yet performed.** `cassia_master_dataset.csv`
is intentionally empty. No value may be added that cannot be traced to a table, figure
or supplementary file of a cited publication.

## 1. Eligibility

Include:
- Batch experiments on Cr(VI) removal from water by activated carbon, biochar, or
  biomass-derived porous carbon (raw, activated or chemically modified, as long as
  carbon is the main phase).
- Reported: C0, residual Cr(VI) (or total Cr, with the analytical method named), pH,
  adsorbent dose, contact time and temperature.
- Peer-reviewed, any year, English.

Exclude:
- Composites in which a non-carbon phase is the declared active component
  (e.g. nZVI, layered double hydroxides, MOFs on carbon).
- Column-only studies (kept in a separate breakthrough register for later work).
- Data that exist only as fitted parameters with no underlying observations.
- Duplicate publication of the same data set (keep the most complete record).

## 2. Search

Databases: Scopus, Web of Science Core Collection, PubMed; data repositories: Zenodo,
Mendeley Data, Figshare, Data in Brief.

Core string (adapt per database syntax):

```
("hexavalent chromium" OR "Cr(VI)" OR chromate OR dichromate)
AND (adsorption OR biosorption OR removal OR reduction)
AND ("activated carbon" OR biochar OR "porous carbon" OR "carbon material" OR hydrochar)
AND (isotherm OR kinetic* OR "batch")
```

Supplementary data filter: add `AND (supplementary OR "data set" OR dataset)` for the
repository search. Record the date and hit count for every query.

## 3. Screening

Two-stage (title/abstract, then full text), with the reason for each exclusion recorded
in `study_registry.csv`. Report the counts in a PRISMA 2020 flow diagram.

## 4. Extraction rules

1. Tables and supplementary files first; digitise figures only when no table exists
   (WebPlotDigitizer, recording the figure number; digitised rows have
   `extraction_method = digitised`).
2. Units are converted to the schema units (`cassia_schema.csv`); T is always in K.
3. `pH_basis` records whether the paper reports initial, final or controlled pH. If
   unstated, leave it blank. Never infer it.
4. `Ce_total` only when total Cr was measured (AAS/ICP). `Ce_VI` only when Cr(VI) was
   measured specifically (DPC/IC-PCR). If the paper states "chromium" without the
   method, record under `method` and do not fill either column until resolved.
5. Unreported ionic strength, Eh, alkalinity and major ions stay blank. They are handled
   downstream as labelled scenario ranges (Section 2.4 of the manuscript), never
   imputed into the data file.
6. Each row carries `doi` and `source_location`.

## 5. Evidence classes (computed by `cassia.dataset.classify`)

| Class | Requirement | Permitted endpoint |
|---|---|---|
| A | Ce_total and Ce_VI measured, solid-phase Cr speciation (XPS/XANES), pH basis final or controlled, background electrolyte reported | Adsorbed Cr(VI), reduced Cr and Cr(III) partition |
| B | Ce_total **or** solid-phase speciation | Cr(VI) disappearance with a bounded partition |
| C | Residual Cr(VI) only | Apparent Cr(VI) disappearance only |

## 6. Quality appraisal

Per study, record replicates (n), SD, a blank or control without adsorbent, a mass
balance check and pH drift reporting. These are used as weights and in sensitivity
analyses, not as exclusion criteria.

## 7. Validation units

Leave-one-study-out and leave-one-material-family-out splits are defined on `study_id`
and `material_family` before any model is fitted. At least one complete study is held
back, unopened, as the external test set.

## 8. Execution records (to be completed when the search is run)

- `search_log.csv`: one row per query (date, database, exact string, filters, hits, export file, operator).
- `prisma_flow.csv`: counts for every PRISMA 2020 stage, including duplicates removed and exclusion reasons.
- `study_registry.csv`: one row per screened full text with inclusion decision and reason.
- Extraction QC: each included study is extracted by one person and checked by a second against the
  source table/figure; for digitised figures, 10% of points are re-digitised and the mean absolute
  difference is recorded in the registry. Disagreements are resolved by discussion and logged.

Status on 2026-10-30: not executed. The publisher and repository hosts were not reachable from the
computational environment in which CASSIA 0.2 was built; no records have been screened or extracted.
