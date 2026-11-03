# Response to the Reviewer (ninth round)

**Manuscript:** CASSIA: A Provenance-Aware Computational Framework and Known-Truth Benchmark for Evidence-Constrained Analysis of Cr(VI) Removal by Carbonaceous Materials (revision 9; software CASSIA 0.3.7, code commit `e165fa260c9c8c506bf1e6b3ee3ba65db104f57e`)

We thank the reviewer for verifying the archive manifest and the M2s implementation, and for finding the mismatch between the family-specific point estimates and their bootstrap intervals. This revision makes the estimator and its interval consistent, and states the reproduction that has actually been run. The items that need the author's own details, accounts or computing time are listed at the end; we have not filled them in.

**Summary of changes.**

| Comment | Change | Where |
|---|---|---|
| Major 1 Family-specific estimand | Estimand defined as the pooled rate; the cluster bootstrap now recomputes that pooled rate in each resample; Table S5c regenerated | Table S5c, §3.5, `scripts/m2s_exchangeability.py` |
| Major 2 Reproducibility | Wording states what was run and that an uninterrupted end-to-end run is still outstanding; faster stages rerun from a clean archive of this version | §2.2, §2.12, Data and code availability |
| Major 3 Empirical validation | Unchanged positioning; no values introduced | §3.8, §4.7 |
| Minor 1 Table S5b page break | Table rows can no longer split across pages; header rows repeat; S5b layout tightened | All Word tables, Table S5b |
| Minor 2 S5c caption | States the estimand | Table S5c |
| Minor 3 Archive name | Archive file name now gives both the code and the documents commit; the relation is stated | Data availability, `release/RELEASE.md` |
| Minor 4–8 | Require the author | Author actions |

---

## Major Comment 1 – Family-specific point estimates and bootstrap intervals

**Response: accepted; we use the reviewer's option 2.**
- *Estimand.* The quantity of interest is the coverage probability for a test study drawn at random from the role assignments of a given family, pooled over corpora. This matches the marginal, study-level reading of the M2s guarantee, and it is the estimand already used in Table S5b. Each corpus is therefore weighted by the number of role assignments of that family it contributes.
- *Point estimate.* Unchanged: the pooled rate.
- *Interval.* The corpus-cluster bootstrap now resamples corpora with replacement and recomputes the same pooled ratio, i.e. the sum of covered role assignments divided by the number of role assignments, in each of 4,000 resamples. Previously it averaged per-corpus rates with equal weight. `scripts/m2s_family.py` (Table S5b) already used the pooled ratio, so the two tables now use the same estimator.
- *Effect.* The point estimates are unchanged. The intervals changed by at most a few hundredths: for apparent uptake, hydrochar is now 0.85 (0.77–0.91) instead of (0.76–0.91), and manure biochar is unchanged at 0.89 (0.85–0.93). For q_tot,aq, hydrochar is 0.92 (0.81–0.99; 9 corpora). §3.5 now states the estimand, and the Table S5c caption defines it.
- *Alternative estimand.* For transparency, the unweighted mean of the per-corpus rates is written to `T_m2s_exchangeability.csv` (column `family_<family>_equal_corpus_mean`). It differs from the pooled rate by at most 0.024; the largest difference is for porous carbon on C_e,tot. This reproduces the differences the reviewer reported.
- *How the table was regenerated.* Only the summary stage was run (`m2s_exchangeability.py --summarise`), from the saved per-study results; no simulation was repeated. All other columns of the table are unchanged.

---

## Major Comment 2 – End-to-end reproducibility

**Response: accepted. The wording now matches what was run, and the full run remains an author action.**
- *§2.2.* Every number is produced by the seeded scripts of Table 2. The reported results were produced by running these scripts stage by stage, and a single uninterrupted run of `reproduce_all.sh` from the source archive has not yet been completed.
- *§2.12.* The bullet "full regeneration of every reported number from seeded scripts" has been replaced by: generation of every reported number by seeded scripts, and regeneration of the faster stages from a clean copy of the source, which reproduced the committed result tables byte for byte; a complete end-to-end run is still outstanding.
- *Data and code availability.* A new paragraph, "Reproduction status", lists what was run: the tests, the stage-by-stage runs, the clean-copy regeneration of the faster stages, and the M2s evaluation. It states that the uninterrupted run is required before publication and that the manuscript will be corrected if any output differs.
- *For this revision.* From a clean `git archive` extraction of the code commit of this version, we reran the faster stages: `typeI_corpus.py`, `seqver_oc.py`, `bench_partition.py`, `m2s_family.py`, `verification_table.py`, `make_table_s1.py`, the summary stages of `typeI_validation.py` and `m2s_exchangeability.py`, and the supplement generator. All 12 regenerated result files were byte-identical to the committed ones, except `T_m2s_roles_summary.csv`. That file differed only in the last floating-point digit (at most 4 × 10⁻¹⁶), because the summary stage recomputes it from the saved compressed per-study file. The generated supplementary text matched apart from the version label. The 42 tests passed (line coverage 82%).
- *Full run.* We agree with every step the reviewer lists: a clean extraction of the archive, the CPython 3.13.18 environment of `requirements-lock.txt`, the 42 tests in that environment, a complete execution log, and a comparison of the generated tables with the documents. The run takes 9–11 h on 4 cores. Our cloud build environment restarts long-running processes before a run of that length finishes, so it must be done on the author's own machine. It is listed below.

---

## Major Comment 3 – Experimental validation

**Response: agreed; no change of positioning.** No experimental values have been introduced, transcribed, estimated or digitised. The manuscript is positioned as a known-truth methodological benchmark for a computational-methodology journal. It states throughout that it is not an experimentally validated framework.

---

## Minor comments

| No. | Comment | Response |
|---|---|---|
| 1 | Table S5b page break | Every table row in the Word documents is now set so that it cannot break across a page, and header rows repeat on a new page. Column widths now follow the cell contents rather than the header text, and tables with ten or more columns use 8-pt type. Table S5b also merges the study and corpus counts into one column, "studies (corpora)", and shows coverage to two decimals. Its rows no longer wrap or split. |
| 2 | S5c caption | The caption now defines the estimand (Major Comment 1). |
| 3 | Archive name | The archive is a `git archive` of the commit that follows the code commit and adds only the rebuilt documents and the release record. Its code, data, result and figure files are those of the code commit, as the manifest verifies. The archive is now named `cassia-0.3.7-source-code<code commit>-docs<documents commit>.tar.gz`, and Data availability and `release/RELEASE.md` explain this. |
| 4 | DOI | Requires the author's GitHub and Zenodo accounts; `CITATION.cff` keeps the placeholder until then. |
| 5–7 | Title page, declarations, verification statement | Require the author. |
| 8 | Length | About 20,600 words including tables, Appendix A and references. The author will check this against the target journal's limit. |

---

## Items the author must complete before submission

| Item | Where |
|---|---|
| Surname, affiliation, ORCID | Title page (`build_docx.js`), `CITATION.cff`, `LICENSE` |
| Funding statement | Back matter (currently marked [Author to confirm]) |
| CRediT statement and any co-authors | Back matter |
| The author's own statement of verification of code and results | §2.12 (placeholder) |
| Generative-AI declaration, checked against the target journal's policy | Back matter |
| Word limit of the target journal | Whole manuscript |
| GitHub release `cassia-v0.3.7` and Zenodo DOI | Data availability, `CITATION.cff` (steps in `README.md`) |
| Uninterrupted run of `bash scripts/reproduce_all.sh` from a clean extraction of the archive, in the environment of `requirements-lock.txt` (CPython 3.13.18), with the 42 tests and a complete execution log; comparison with the submitted tables | `release/test_log.txt`, results, Data availability "Reproduction status" |
