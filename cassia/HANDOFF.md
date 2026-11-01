# CASSIA – project handoff (cloud session → local environment)

Handoff date: 2026-11-01. This file is the complete state of the work done in the cloud Claude Code session, so you can continue in a local project folder. Read sections 1–4 to get running. Sections 5–9 are the project memory: what was done, why, the key numbers and what is still open.

---

## 1. What this project is

**CASSIA** (*Chromium Adsorption, Speciation and Surface-reduction Identifiability Analysis*) is an open Python framework and a research article. Its current title is:

> *CASSIA: A Provenance-Aware Computational Framework and Known-Truth Benchmark for Evidence-Constrained Analysis of Cr(VI) Removal by Carbonaceous Materials*

- **Target:** an Elsevier journal. The manuscript is formatted to Elsevier review conventions: double spacing, line numbers, highlights of ≤ 85 characters, an abstract of ≤ 250 words, and CRediT and AI-use declarations.
- **Nature of the work:** purely computational. It uses:
  - PHREEQC through phreeqpython;
  - a multi-mechanism known-truth simulator with scenarios A–E;
  - identifiability and equivalence tests;
  - ML transferability (RF, XGBoost and GPR);
  - study-level conformal prediction;
  - two-endpoint treatment design with sequential verification.
- **No experiments and no empirical data.** No laboratory experiments were done, and **no published experimental data have been analysed yet**. The paper is positioned as a *methodological, known-truth benchmark*.
- **Review history.** The manuscript has been through four reviewer rounds; the current software is **CASSIA 0.3.2**.

## 2. Where everything is

**Repository:** `https://github.com/karthice1/adsorption`
- **Working branch:** `claude/exciting-hawking-ayadsw`. All work is pushed there; nothing is merged into `main` and no pull request has been opened.
- **Repository root:** the `adsorption` equation package (Amvro23, MIT), which CASSIA depends on.
- **CASSIA itself:** the `cassia/` directory, which holds everything else.

**Latest commits:**
- `8f60859`: response-letter wording (HEAD before this handoff).
- `b35d448`: revision-4 deliverables.
- `72a0b17`: **code commit cited in the manuscript**; `release/MANIFEST.sha256` was generated at this commit.

**Other repositories** attached to the cloud session (`karthice1/pyGAPS`, `karthice1/RASPA3` and `karthice1/phreeqpython`) were only audited or used as libraries. Nothing was changed in them.

### Directory map (`cassia/`)

| Path | Contents |
|---|---|
| `cassia/` (package) | `dataset.py` (37-field schema, evidence classes C/B/A1/A2), `chemistry.py` (PHREEQC adapter, 4-D speciation grid, SI bisection, electrolytes), `simulate.py` (scenarios A–E, mass/electron balance, `parameter_table()` = Table S1), `measurement.py` (XPS/XANES/inventory/phase models), `identifiability.py` (Fisher, detection tests, equivalence test `equivalence_test_v3`), `fitting.py`, `validation.py` (RF/XGB/GPR, LOSO/LOFO, conformal M0/M1/M2/M2s), `pores.py` (pyGAPS adapter) |
| `scripts/` | `run_v2.py`, `run_v3.py` (stage drivers), `typeI_validation.py`, `seqver_oc.py`, `bench_partition.py`, `endpoint_validation.py`, `m2s_family.py`, `verification_table.py`, `make_table_s1.py`, `make_supp_figs.py`, `make_release_manifest.py`, **`reproduce_all.sh`** (runs everything in order) |
| `tests/test_cassia.py` | 40 tests (mass/electron balance, SI status, interpolation, conformal, equivalence negative control, Table S1 ↔ code constants) |
| `data/` | `PROTOCOL.md` (extraction protocol, class criteria), `study_registry.csv` (seed study S001 = Ramirez et al. 2020), templates, PRISMA/search-log templates |
| `databases/` | `minteq.v4.dat`, `minteq.dat`, `llnl.dat` with SHA256SUMS and PROVENANCE.md |
| `results/v2/`, `results/v3/` | All result tables (CSV/JSON) and per-study outputs; `coverage_report.txt` |
| `release/` | `MANIFEST.sha256` (120 files), `RELEASE.md`, `test_log.txt` |
| `examples/user_data_workflow.py` | End-to-end workflow for real data in the CASSIA schema |
| `manuscript/` | Sources and deliverables (see §4) |
| `CHANGELOG.md`, `CITATION.cff`, `README.md`, `pyproject.toml`, `requirements-lock.txt` | Packaging and release metadata |

## 3. Setting up locally

### 3.1 Get the code

```bash
git clone https://github.com/karthice1/adsorption.git
```

```bash
cd adsorption
```

```bash
git checkout claude/exciting-hawking-ayadsw
```

If you would rather keep the work in a new folder or repository, copy the whole clone. CASSIA expects the `adsorption` package at the repository root, one level above `cassia/`.

As an alternative to git, unpack `cassia/manuscript/cassia-0.3.2-source-8f60859.tar.gz`. It is a `git archive` of commit `8f60859` and was sent to you in the chat; it is not stored in the repository.

### 3.2 Python environment (Python 3.10+; the cloud used 3.13.17)

```bash
python -m venv .venv
```

```bash
source .venv/bin/activate
```

On Windows, activate with `.venv\Scripts\activate` instead.

```bash
pip install -e .
```

```bash
cd cassia
```

```bash
pip install -e ".[ml,pores,dev]"
```

For the exact versions used in the cloud, install from the lock file:

```bash
pip install -r requirements-lock.txt
```

Key versions: phreeqpython 1.6.2, NumPy 2.5.3, SciPy 1.18.1, pandas 3.0.5, scikit-learn 1.9.1, XGBoost 3.4.1 and matplotlib 3.11.2. pyGAPS was a 4.3.1 development build; `pip install pygaps` (4.x) is sufficient.

Check the installation:

```bash
pytest -q tests
```

Expected result: `40 passed`, in about 40 s.

```bash
sha256sum -c --quiet release/MANIFEST.sha256
```

This should print nothing, which means every file matches the manifest.

### 3.3 Word-document build tools

- **Node.js 18+** (the cloud used v22) with the `docx` package, version 9.7.1:

```bash
cd manuscript
```

```bash
npm install
```

- **LibreOffice:** optional, only for rendering the `.docx` files to PDF for visual checks. The cloud used 24.2.

## 4. How the manuscript is produced

The `.docx` files are **generated from Markdown sources**. Edit the sources, never the `.docx` files.

| Source file (in `manuscript/`) | Section |
|---|---|
| `_v3_front_methods.md` | Title, keywords, §1 Introduction, §2 Methods (Tables 1–2) |
| `_v3_results_A.md` | §3.1 chemistry, §3.2 benchmark (Table 3–4), §3.3 identifiability/equivalence (Table 5, Fig. 3) |
| `_v3_results_ML.md` | §3.4 ML (Table 6, Figs 4–5) |
| `_v3_results_B.md` | §3.5 intervals (Tables 7, 7b), §3.6 design (Tables 8, 8b, 8c, Fig. 7), §3.7 missing chemistry (Fig. 8), §3.8 empirical status |
| `_v3_discussion.md` | §4 Discussion, §4.7 Limitations, §5 Conclusions (ML discussion is inlined) |
| `_v3_appendix.md` | Appendix A (M2s theorem/proof, LaTeX math), Appendix B (GOF/equivalence/type-I calibration) |
| `_v3_back_matter.md` | Data/code availability, funding, CRediT, declarations, references |
| `assemble_v3.py` | **Highlights and Abstract live here** (lists `HL` and `AB`); assembles everything into `CASSIA_manuscript.md` and `_build_manuscript.md` (appendix math converted by `latex2uni.py`) |
| `make_supplement_v3.py` + `_supp_base_v2.md` | Builds `CASSIA_supplementary.md` (Table S1 from code, S2–S5 frozen, S5b–S14, Fig. S1) |
| `_v5_response.md` | Source of the latest response letter (placeholders `{COMMIT}` etc. are filled into `CASSIA_response_to_reviewer.md`) |
| `_round1_…`, `_round2_…`, `_v3_response.md`, `_v4_response.md` | Earlier response letters (archive) |
| `build_docx.js` | Markdown → Word (Elsevier review format: Times 12 pt, double spacing, continuous line numbers, page numbers; tables, figures, bullets, sub/superscripts `_{}`/`^{}`, centred equations `<center>`, fenced code) |
| `postprocess_docx.py` | Stops line numbers appearing in the footer |
| `figures/` | fig1–fig8 PNG and figS1 (file names are historical: `fig8_mechanism_transfer.png` is **Fig. 5**, `fig5_conformal.png` is Fig. 6, `fig6_design.png` Fig. 7, `fig7_missing_chemistry.png` Fig. 8) |

### Rebuild all four Word files

Run these from `cassia/manuscript/`.

**1. Supplement and manuscript.** Pass the code commit you want cited:

```bash
python make_supplement_v3.py
```

```bash
python assemble_v3.py "$(git rev-parse HEAD)"
```

**2. Response letter.** Fill its placeholders from `_v5_response.md`:

```bash
python -c "s=open('_v5_response.md').read();c='$(git rev-parse HEAD)';s=s.replace('{COMMIT}',c).replace('{COMMIT_SHORT}',c[:7]).replace('{NTESTS}','40').replace('{COV}','81%');open('CASSIA_response_to_reviewer.md','w').write(s)"
```

**3. Word documents.** On Windows PowerShell, first run `$env:CASSIA_IN="_build_manuscript.md"`, then run `node build_docx.js` without the prefix:

```bash
CASSIA_IN=_build_manuscript.md node build_docx.js
```

This writes `CASSIA_manuscript.docx` and `CASSIA_highlights.docx`.

```bash
node build_docx.js plain CASSIA_supplementary.md CASSIA_supplementary.docx
```

```bash
node build_docx.js plain CASSIA_response_to_reviewer.md CASSIA_response_to_reviewer.docx
```

```bash
for f in CASSIA_manuscript CASSIA_supplementary CASSIA_response_to_reviewer CASSIA_highlights; do python postprocess_docx.py $f.docx; done
```

**4. Optional visual check:**

```bash
soffice --headless --convert-to pdf CASSIA_manuscript.docx
```

**Release workflow, in order:**
1. Commit the code and results.
2. Run `python scripts/make_release_manifest.py 0.3.2` from `cassia/`.
3. Rebuild the documents with that commit.
4. Commit again.

The manuscript cites the *code* commit, the one at which the manifest was generated.

### Re-run the analyses

From `cassia/`:

```bash
bash scripts/reproduce_all.sh
```

This takes about 9–11 h on 4 cores. The steps can also be run one at a time; see manuscript Table 2 and the comments inside the script.

**Ordering pitfall:** `run_v2.py --stage all` overwrites Figs. 3, 5, 6 and 7 with older versions. `run_v3.py --stage figures` must therefore run **after** it, and `reproduce_all.sh` already does this.

**Two-hour job limit (cloud only):** the slow stages, `equiv3` (`EQUIV_PART=1/2`) and `typeI_validation.py` (`CTRL_PART=1/2`), were split into parts to fit the cloud's 2-hour limit. Locally they can run straight through.

## 5. Review history (what each round asked and what was done)

| Round | Main requests | Outcome |
|---|---|---|
| Original | Full article from repository audit, PHREEQC, simulator, ML, UQ, design; Word in Elsevier format | v0.1 manuscript |
| 1 | Empirical validation; circular simulator (fixed biochar reductance) → multi-mechanism A–E, ≥100 replicates; PHREEQC (4-D speciation grid, SI not phase-boundary convergence, database validation, uncertainty propagation with more draws); treatment design (separate total-Cr model, joint criteria, more materials); study-level conformal; ML design; protocol; software release | v0.2: all addressed except empirical data |
| 2 | Ramirez + 20–30 studies; formal equivalence test (margins, TOST, calibrated type-I, power); classes C/B/A1/A2 with separate XPS/XANES models; sequential-verification false-acceptance; exact theorem for conformal; ML across scenarios; PHREEQC pH-threshold fix; release/DOI | v0.3.0: done except empirical data, Ball & Nordstrom validation, DOI |
| 3 | Repeat of round-1 comments (cited v0.1 pages) | v0.3.1: benchmark precision (Tables S11–S12), endpoint-model validation (Table S13), unit-of-calibration table, removed revision-history wording |
| 4 | (1) empirical still missing; (2) Table S1 out of date; (3) Table S5 lacks M2s; (4) type-I not controlled (D4n 12.5%); (5) D5 phase evidence vs test vector; (6) denominators of safety claim; (7) endpoint-specific conformal wording; (8) version/reproducibility inconsistencies | v0.3.2: all addressed; (1) by repositioning as a methodological benchmark |

## 6. Key results (all synthetic or thermodynamic)

**Chemistry.**
- The three PHREEQC databases agree on Cr(VI) speciation within 0.044, but their Cr(VI)/Cr(III) boundaries differ by 63–83 mV (minteq.v4 against llnl).
- Cr(III) stays undersaturated up to 5,000 mg L⁻¹ to pH 4.25 (minteq.v4, I ≤ 0.01 M), 4.5 (I = 0.1 M) and 3.5 (minteq).
- Changing the electrolyte (NaCl, NaNO₃ or Na₂SO₄) shifts species fractions by at most 0.024.

**Benchmark** (15,000 virtual studies; 100 corpora per scenario).
- Median reduction shares are 0.24 (B), 0.22 (C), 0.04 (D) and 0.38 (E).
- Fits with R² ≥ 0.95 occur in 41–52% of high-reduction studies (B, C, E).
- Fitted q_max exceeds the adsorption capacity in proportion to the reduction share.
- MCSE is at most 0.02 (Tables S11–S12).

**Equivalence and rejection tests** (Table 5, Table S14).
- For Class C data, adsorption-only was not rejected for 65–78% of high-reduction studies. It was shown equivalent within 10% for 43–73%.
- Power is 80–88% with total Cr and 92–98% with XPS or XANES.
- The type-I error is 4.7–5.6% for every design, calibrated with a conformal threshold on 360 controls from 18 corpora and evaluated leave-one-corpus-out (95% intervals within 2.5–8.6%). In single corpora it reaches 15–25%. The old six-corpus calibration gave 14–15% for D1.
- The false-equivalence rate is 0–4.3%.

**Model-free diagnostics.** A simple XPS threshold false-alarms in 7.1% of controls, against 1.0% for XANES. Only A2 designs constrain the precipitated share (posterior SD 0.026).

**Machine learning.**
- R² falls from 0.89 under random splits to 0.63 leave-one-study-out (LOSO) and 0.56 leave-one-family-out (LOFO). Across scenarios, LOSO R² is 0.48–0.65.
- Speciation features do not help.
- The intraclass correlation (ICC) is 0.45.
- The learning curve rises from 0.31 (12 studies) to 0.66 (48–96 studies).
- Endpoint models under LOSO reach log-R² 0.74 (C_e,VI) and 0.65 (C_e,tot).

**Prediction intervals.**
- The pooled M0 interval covers only 65–75% of studies.
- M2s is exact, covering 0.93 of studies at nominal 0.90. Its 90% width is 3.3 log₁₀ units for q_app and 3.9 for C_e,VI. For total-Cr endpoints it is unbounded at 95% and often at 90%.
- Hydrochar reaches 84% under M2s (42/50; Wilson 72–92%).

**Treatment design.**
- Point predictions were unsafe for 20.7% (S1) and 7.0% (S2) of materials.
- The simple verification rule accepted 105 (S1) and 98 (S2) of 300 materials, with 0 unsafe. The Wilson upper bound is 1.3% over all materials and 3.5–3.8% among accepted designs.
- On an independent cohort, the error-controlled rule kept the false-acceptance probability at or below 0.7% per material.

## 7. Open items (in priority order)

1. **Empirical application – the decisive gap.** The reviewer will not recommend acceptance in a leading remediation journal without it.
   - **Seed dataset:** Ramirez, A., Ocampo, R., Giraldo, S., Padilla, E., Flórez, E., Acelas, N. (2020), *Data in Brief* 29, 105292, doi:10.1016/j.dib.2020.105292. It is open access as **PMC7044665**, with raw isotherm and kinetic files as supplementary material. The companion article is *J. Environ. Chem. Eng.* 8 (2020) 103702.
   - **Why it was not done in the cloud:** the cloud network policy blocked every literature host. Locally, download the PMC supplementary files directly.
   - **Plan:** reviewer steps 1–7, listed in §3.8 of the manuscript.
     1. Import the data in the CASSIA schema (`data/` templates; `examples/user_data_workflow.py`).
     2. Record provenance and DOI.
     3. Classify as Class C.
     4. Compute apparent uptake.
     5. Fit nonlinear isotherms and kinetics.
     6. Compare with the published Langmuir q_max of 72.46 mg g⁻¹.
     7. Extract 20–30 further studies, prioritising those with total Cr, XPS/XANES, solid inventory, final pH and electrolyte data.
   - **Presentation:** report the empirical findings in a section separate from the synthetic results.
   - **Integrity rule (kept throughout):** never transcribe numbers from memory, estimate them or invent them. Digitise from figures only if it is declared as such.
2. **Independent thermodynamic validation** against Ball & Nordstrom (1998) and Rai et al. (1987), which were also blocked in the cloud. CASSIA accepts a reference-constant table for this comparison.
3. **Release:**
   - create the GitHub tag `cassia-v0.3.2` and connect Zenodo to mint a DOI;
   - update `CITATION.cff`, which still holds a placeholder DOI `10.5281/zenodo.0000000`;
   - update the Data availability text in `_v3_back_matter.md`.
4. **Author details:**
   - surname, affiliation and ORCID: placeholders on the title page are hard-coded in `build_docx.js` (`Karthik [Surname]`, `[Department, …]`, `ORCID: [to be added]`);
   - CRediT, funding and co-authors: in `_v3_back_matter.md`;
   - confirm the AI-use declaration against the target journal's policy.
5. **Optional clean-ups:**
   - rename the figure files to match their numbers;
   - fix subscript artefacts in supplement Table S3, where underscores in parameter names render as subscripts;
   - main-text Table 5 is cramped in Word.

## 8. Conventions and pitfalls learned

- **Honesty rules.**
  - Everything is labelled synthetic or thermodynamic.
  - "Not rejected" is never read as "equivalent"; there are three outcomes: rejected, equivalent within ε, and inconclusive.
  - Guarantees are stated only where they are proven: M2s (Appendix A) and the conformal threshold under exchangeability.
  - Observed frequencies are given with Wilson intervals and their denominators.
- **The corpus is the unit of exchangeability** for null distributions and calibration, because hyper-priors are redrawn per corpus. Calibrating on a few corpora does not transfer.
- **Table S1 must always be regenerated from code** (`scripts/make_table_s1.py`); a test enforces this.
- **Supplement builder** (`make_supplement_v3.py`):
  - it reads the frozen `_supp_base_v2.md` and is idempotent, so do not hand-edit `CASSIA_supplementary.md`;
  - `|` inside table cells breaks the docx table parser, and the builder replaces it.
- **Table numbering:** main-text Tables 1–8 (with 7b, 8b, 8c); supplementary S1–S14, plus S5b, S6b and Fig. S1. Adding a main table means renumbering; the round-4 renumbering was done with a regex over the `_v3_*.md` files.
- **Manuscript text must not mention the review history** (no "v0.1", "v0.2" or "the reviewer"); that belongs only in the response letters.
- **Commit trailer** used in the cloud (optional locally): `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## 9. Starter prompt for a new local Claude Code session

Paste this into Claude Code from the `adsorption/` folder. `cassia/CLAUDE.md` also loads automatically.

> I am continuing the CASSIA research article (Cr(VI) removal by carbonaceous materials). Read `cassia/HANDOFF.md` first; it describes the repository, the build process, the review history and the open items. The manuscript is generated from `cassia/manuscript/_v3_*.md` via `assemble_v3.py` and `build_docx.js` — edit the sources, not the .docx files. The next task is the empirical application (handoff §7.1): I have downloaded the Ramirez et al. (2020) supplementary files into `<path>`. Import them into the CASSIA schema, run the Class C workflow, compare with the published results, and add a separate empirical results section, without inventing or estimating any values. Then rebuild the four Word documents and update the response letter.
