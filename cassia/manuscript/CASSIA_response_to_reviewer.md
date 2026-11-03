# Response to the Reviewer (eighth round)

**Manuscript:** CASSIA: A Provenance-Aware Computational Framework and Known-Truth Benchmark for Evidence-Constrained Analysis of Cr(VI) Removal by Carbonaceous Materials (revision 8; software CASSIA 0.3.6, code commit `566e59fd2287a5144612a86f3d13a384e8db1879`)

We thank the reviewer for checking the source archive and the role-averaged M2s implementation independently, and for the clear statement of what remains. This revision makes the wording corrections requested in Section 4 of the report and shortens the manuscript for a computational-methodology journal. The items that need the author's own details, accounts or computing time are listed at the end; we have not filled them in.

**Summary of changes.**

| Comment | Change | Where |
|---|---|---|
| 2.1, 4.1–4.2 Scope of the M2s guarantee | The Abstract, Conclusions and Appendix A now say that the guarantee holds under random role assignment for a randomly held-out study of the same realised corpus. They call it a finite-corpus randomisation guarantee, not a guarantee for a study from an independent literature population. | Abstract, Highlights, §5.1 item 5, Appendix A |
| 2.4 Empirical application | Unchanged; no values introduced | §3.8, §4.7 |
| 3 End-to-end reproduction | Still requires a local run; see below | Author actions |
| 4.3 Table S5c clusters | Number of contributing corpora shown for each family interval; captions of Tables S5b and S5c say which corpora are resampled | Table S5b, Table S5c, `scripts/m2s_exchangeability.py` |
| 4.4 Python version | Response letter harmonised with the lock file (CPython 3.13.18) | This letter |
| 4.5 Width scale | Widths stated on the log₁₀(y + c) scale; the factors between bounds refer to y + c | §3.5, Tables 7, 7b, S5b, S5d |
| 4.6 Length | Shortened for a computational-methodology journal by about 1,800 words: Appendix B, the full propagation of unreported chemistry and the sequential-verification details moved to Supplementary Notes S1–S3 | §2.10, §3.6, §3.7, Supplementary Notes S1–S3, Table S15, Fig. S2 |
| 4.7–4.8 DOI and declarations | Require the author | Author actions |

---

## 2.1 and Additional corrections 1–2 – Scope of the guarantee in the Abstract and Conclusions

**Response: accepted.**
- *Abstract.* The sentence "Observed coverage of an exact study-level split-conformal interval agreed with its guarantee" has been replaced by: "Under random role assignment, an exact study-level split-conformal interval guarantees coverage only for a randomly held-out study of the same realised corpus, not for independent literature". To stay within 250 words, three details given in full in §3.3 and §3.6 were dropped from the Abstract: the number of corpora behind the negative controls, the rounded "about 5%" beside the 4.7–5.6% range, and the Wilson bound among accepted designs.
- *Conclusions (§5.1 item 5).* The item now reads that the coverage was consistent with the guarantee under random role assignment. It adds that this is a finite-corpus randomisation guarantee for a randomly selected held-out study of the same realised corpus, marginal over family and study size, and not a guarantee for a genuinely new study from an independent literature population.
- *Appendix A.* The scope paragraph after the theorem uses the same words: a finite-corpus randomisation guarantee, not a conventional distribution-free guarantee for future experimental studies.
- *Highlights.* The coverage highlight now reads "Exact coverage is marginal within one corpus, not per family; intervals are too wide" (84 characters).

---

## 2.4 – Empirical application

**Response: agreed; no change of positioning.** No experimental values have been introduced, transcribed, estimated or digitised. The Ramirez et al. (2020) source files still cannot be retrieved from our computational environment. The manuscript remains a known-truth methodological benchmark and is now aimed at a computational-methodology journal, in line with the reviewer's scope-dependent recommendation.

---

## 3 – Clean end-to-end reproduction

**Response: agreed; this remains an author action.**
- *Status.* The complete `reproduce_all.sh` run takes 9–11 h on 4 cores. Our cloud build environment restarts long-running processes before a run of this length can finish, so we could not complete it there. As in revision 7, the faster stages were rerun from a clean copy and matched the committed results byte for byte.
- *This revision.* No analysis was rerun except the summary step of the M2s evaluation (`m2s_exchangeability.py --summarise`). This step rebuilds Tables 7b and S5c from the saved per-study results and adds the corpus counts. All previously reported values in `T_m2s_exchangeability.csv` were unchanged.
- *Tests.* The test log records 42 passed (line coverage 82%). We note that the reviewer treats it as author-reported evidence, and the final run on the author's machine (below) will provide an independent log.

---

## Additional corrections

| No. | Comment | Response |
|---|---|---|
| 1–2 | Abstract and Conclusions wording | Done (section 2.1 above). |
| 3 | Table S5c cluster count | Each family-conditional rate is now followed by its bootstrap interval and the number of corpora contributing to it, e.g. "0.917 (0.79–0.99; 9 corpora)" for hydrochar on q_tot,aq. Only corpora that contain test studies of the family contribute, and these are the clusters resampled. That is 10 for every family and endpoint except hydrochar on the two total-Cr endpoints: total Cr is measured only in Class A/B studies, and one corpus has no hydrochar study with total-Cr data. The Table S5b caption now says the same. The count is written by `scripts/m2s_exchangeability.py` (column `family_<family>_ncorpora`). |
| 4 | Python version in the letter | Harmonised. The lock file records the interpreter actually used, CPython 3.13.18, with every package version unchanged from the earlier 3.13.17 environment. The author action below asks for a final run in the environment of the lock file, not in 3.13.17. |
| 5 | Width scale | Accepted. §3.5 now says that widths are on the modelling scale log₁₀(y + c) of §2.7, with c = 0.1 mg g⁻¹ for uptake, so that a width of 3.5 units corresponds to a factor of about 3,000 between the bounds of y + c, not of y itself. The same applies to the M0–M2 comparison (factor of 100–1,000), Tables 7 and 7b and the supplementary width columns. |
| 6 | Length | Accepted. The manuscript is aimed at a computational-methodology journal and has been shortened from about 22,400 to about 20,600 words including tables, Appendix A and references. The goodness-of-fit and equivalence tests (formerly Appendix B) are now Supplementary Note S1. The methods and results of the propagation of unreported chemistry are summarised in one paragraph each in §2.10 and §3.7, with the full text and the former Fig. 8 in Supplementary Note S2 (Fig. S2). The sequential-verification details of §3.6 and the former Table 8b are in Supplementary Note S3 (Table S15); §3.6 keeps the main numbers and Table 8c. The duplicate endpoint-accuracy table in §3.6 was removed in favour of Table S13. Section numbers are unchanged; all cross-references were updated. |
| 7 | Release and DOI | Require the author's GitHub and Zenodo accounts (below). |
| 8 | Title page and declarations | Require the author (below). |

---

## Items the author must complete before submission

We have not filled these in, because they require the author's own information, judgement, accounts or computing time.

| Item | Where |
|---|---|
| Surname, affiliation, ORCID | Title page (`build_docx.js`), `CITATION.cff`, `LICENSE` |
| Funding statement | Back matter (currently marked [Author to confirm]) |
| CRediT statement and any co-authors | Back matter |
| The author's own statement of verification of code and results | §2.12 (placeholder) |
| Generative-AI declaration, checked against the target journal's policy | Back matter |
| GitHub release `cassia-v0.3.6` and Zenodo DOI | Data availability, `CITATION.cff` (steps in `README.md`) |
| Clean end-to-end run of `bash scripts/reproduce_all.sh` from the supplied archive, in the environment of `requirements-lock.txt` (CPython 3.13.18), including the test suite | `release/test_log.txt`, results |
