# CASSIA – instructions for Claude Code

Start by reading `HANDOFF.md` (full project state, build steps, open items).

- Research integrity: never invent, estimate or transcribe-from-memory any experimental value or citation. Empirical data enter only from files the user provides or that are downloaded and recorded with provenance (`data/study_registry.csv`). Label every result as empirical, synthetic or thermodynamic.
- Manuscript: edit `manuscript/_v3_*.md`, `assemble_v3.py` (Highlights/Abstract), `make_supplement_v3.py`, `_v5_response.md`; never edit the generated `.docx` or `CASSIA_*.md` by hand. Rebuild with the commands in HANDOFF.md §4. Highlights ≤ 85 characters, Abstract ≤ 250 words.
- Manuscript text must not refer to review history ("v0.x", "the reviewer"); that goes in the response letter only.
- Table S1 is generated from code (`scripts/make_table_s1.py`); keep `tests/test_cassia.py` passing (`pytest -q tests`, 40 tests).
- Analyses: `bash scripts/reproduce_all.sh` (order matters: run_v3 figures after run_v2). Seeds are fixed; report Monte Carlo uncertainty (Wilson / cluster bootstrap) with any observed rate and state its denominator.
- After code/result changes: commit, `python scripts/make_release_manifest.py <version>`, rebuild documents citing that commit, commit again.
