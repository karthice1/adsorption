# Response to the Reviewer (seventh round)

**Manuscript:** CASSIA: A Provenance-Aware Computational Framework and Known-Truth Benchmark for Evidence-Constrained Analysis of Cr(VI) Removal by Carbonaceous Materials (revision 7; software CASSIA 0.3.5, code commit `{COMMIT}`)

We thank the reviewer for a precise report. The main point is that the M2s results we reported came from one fixed-seed calibration split per corpus, whereas the guarantee we proved in Appendix A rests on random role assignment. We agree that the evaluation should match the theorem. In this revision:
- the reported M2s results are now the random role assignment itself, averaged over 50 proper-training sets per corpus with every remaining study as test study, at every nominal level;
- the original leave-one-study-out run is reported only as one realisation of that procedure;
- the scope of the guarantee is stated in the reviewer's words wherever it appears.

The other comments are addressed point by point below. Items that need the author's own details or accounts are listed at the end; we have not filled them in. Section, table and figure numbers refer to the revised manuscript.

**Summary of changes.**

| Comment | Change | Where |
|---|---|---|
| 1 M2s evaluation | Reported M2s results are role-averaged (50 proper-training sets per corpus, all levels); the theoretical rate k/(J + 1) is shown beside them; the fixed-seed run is labelled a single realisation | §2.8, §3.5, Table 7b, Fig. 6, Tables S5b–S5d, Appendix A |
| 2 Scope wording | "A randomly selected held-out study from the same realised corpus, marginally over family and study size" | §2.8 (table of constructions), §3.5, §4.5, §4.7 item 7, Appendix A |
| 3 Type-I error | 15–25% single-corpus rates kept in the Abstract, Highlights, §3.3 and §4.7; recalibration on real negative controls required; extra-binomial results described as nominal and exploratory | Abstract, §3.3, §4.7 item 5 |
| 4 Empirical application | Unchanged; no values introduced | §3.8, §4.7 |
| 5 What each result depends on | New paragraph separating mathematical properties, simulator-dependent results and results that also depend on the assumed analytical error models | §4.5 |
| 6 PHREEQC | Database agreement no longer presented as evidence of accuracy; electrolyte statement restricted to the electrolytes tested; independent validation remains an open limitation | §3.1, §4.7 |
| 7 Coverage and width | Width and unbounded share reported beside every coverage value; intervals described as not practically reliable for treatment design; design rests on confirmatory testing | §3.5, Table 7b, §5.1 |
| Minor | §4.4 items per evidence class; Appendix A wording; §2.8 fixed-split clarification; Fig. 6 caption corrected; response-letter typo; length | see table below |

---

## Major Comment 1 – Align the M2s evaluation with the random-role procedure

**Response: accepted.** The reviewer is right that `StudySplitConformal(seed=k)` used the same fixed seed in every leave-one-study-out fold of a corpus, so each fold had one calibration subset. That run is one realisation of the random role assignment, not the procedure the theorem describes.

**What the reported M2s results now are (§2.8, §3.5, Table 7b, Table S5d, Fig. 6).**
- *Procedure.* In each of the 10 corpora we drew 50 proper-training sets uniformly and independently (half of the studies other than the test study). For each, every remaining study was taken in turn as the test study, calibrated on the others. This is a Monte Carlo average over uniformly random (test, proper-training, calibration) partitions, which is the setting of the proposition in Appendix A.
- *Levels and endpoints.* All five nominal levels (0.50, 0.70, 0.80, 0.90, 0.95) and all four endpoints, with 7,500 role assignments for each Cr(VI) endpoint and 4,650 for each total-Cr endpoint.
- *Results.* For apparent uptake the fraction of test studies with at least 90% of their observations covered was 0.53, 0.73, 0.80 and 0.93 at nominal 0.50, 0.70, 0.80 and 0.90, and the interval was unbounded at 0.95. The mean width of the 90% interval was 3.46 log₁₀ units (3.71 for C_e,VI, 3.68 and 2.10 for the total-Cr endpoints, bounded intervals only). For the total-Cr endpoints the 90% rate was 0.97, of which 62% of intervals were unbounded.

**Theory and realisation are distinguished explicitly.**
- *Theoretical rate.* Once every role assignment is enumerated, the observed rate equals k/(J + 1) exactly when scores are distinct. Table 7b and Table S5d show this value beside the observed rate, and §3.5 and Appendix A say plainly that the agreement is an identity: it confirms the implementation, and it is not independent evidence that the intervals are useful.
- *Single fixed-seed realisation.* The original run is retained as the last column of Table 7b (0.93, 0.95, 0.97 and 0.98 at 90%) and as Table S5b, and both are labelled as one realisation of the procedure. The methods (§2.8) say that the original run drew one calibration subset per fold with the same fixed seed in every fold of a corpus.
- *Family-conditional and new-corpus rates.* These now use the same 50 role assignments per corpus, with corpus-cluster bootstrap intervals (Table S5c). For apparent uptake: hydrochar 0.85 (95% interval 0.76–0.91), manure biochar 0.89 (0.85–0.93), other families 0.92–0.98. With model and calibration from one corpus and test studies from another, the rate averaged 0.93 and was 0.84 in the worst corpus pair.

**Code.** `scripts/m2s_exchangeability.py` is now the main M2s evaluation (Table 2) and writes `results/v3/m2s_roles_per_study.csv.gz`, `T_m2s_roles_summary.csv` and `T_m2s_exchangeability.csv`. Fig. 6 is drawn from the role-averaged results. `reproduce_all.sh` runs this script before the figures.

---

## Major Comment 2 – "New study" wording

**Response: accepted.** Every statement of the guarantee now uses the reviewer's formulation: a randomly selected held-out study from the same realised corpus, marginally over family and study size. This applies to the table of constructions in §2.8, §3.5, §4.5, §4.7 item 7 and Appendix A. Appendix A adds an explanatory paragraph: the theorem says nothing about a genuinely new study drawn from an independent literature population, nor about a given family. The "In words" sentence after the theorem now refers to the observations of the test study.

---

## Major Comment 3 – Type-I error

**Response: accepted.**
- *Visible single-corpus rates.* The 15–25% maximum single-corpus rates remain in the Abstract ("reaching 15–25% in some individual corpora: an empirical average, not a guarantee"), in §3.3 (the Highlight says "higher rates in some individual corpora"), in bold in the Table 5 footnote, and in §4.7 item 5.
- *Not a reliable 5% test.* §3.3 states that the average empirical type-I error is not evidence of reliable 5% control for real studies; §3.3 and §4.7 state that the threshold must be recalibrated on real negative controls, i.e. studies known to involve adsorption only, before it is applied to literature data.
- *Extra-binomial results.* The two results with p < 0.05 are now described as nominal, unadjusted p-values among 12 comparisons, none significant after Bonferroni correction (0.05/12 = 0.0042), and as exploratory.

---

## Major Comment 4 – Empirical application

**Response: agreed; no change of positioning.** No experimental values have been introduced, transcribed or estimated. The Ramirez et al. (2020) dataset remains Class C and its demonstration will be completed when its source files are available; the computational environment still cannot reach the hosts that serve them. The manuscript keeps the known-truth benchmark scope throughout.

---

## Major Comment 5 – Which results depend on what

**Response: accepted.** §4.5 has a new paragraph, "What each result depends on", with three groups:
- *Properties of the methods,* which follow mathematically and hold for any data meeting the stated assumptions: the M2s statement under random role assignment, the corpus-level conformal bound and the 19-corpus requirement, the minimum numbers of calibration studies, the absence of a study-level guarantee for pooled calibration, and the logic of "not rejected" versus "equivalent".
- *Results that depend on the kinetic simulator:* reduction shares, Class C equivalence, q_max inflation, family patterns, the random-versus-study R² gap and ICC, the size of the study-pooled type-I errors, family-conditional and new-corpus M2s coverage, unsafe-design rates.
- *Results that also depend on the assumed analytical error models:* rejection power with total Cr and with solid-phase data, XPS and XANES false-alarm rates, the posterior SD of the precipitated share, and sequential-verification error rates.

The paragraph ends by stating that only the first group transfers to real data without further evidence.

---

## Major Comment 6 – PHREEQC validation

**Response: partly accepted.**
- *Independent validation.* We could not obtain the Ball and Nordstrom (1998) or Rai et al. (1987) data in machine-readable form in our environment, and we have not transcribed values. Independent validation against critically evaluated Cr(VI)/Cr(III) data therefore remains an open limitation (§4.7).
- *Agreement is not accuracy.* §3.1 now says that agreement between the three databases is not evidence of thermodynamic accuracy, and that an independent comparison with critically evaluated constants is outstanding.
- *Electrolyte scope.* The electrolyte-sensitivity statement is restricted to the three sodium electrolytes tested; Ca²⁺, Mg²⁺, carbonate and dissolved organic matter were not evaluated.

---

## Major Comment 7 – Coverage, width and practical use

**Response: accepted.**
- *Width beside coverage.* Table 7b gives the mean width at 90% next to every coverage value, and the unbounded share is written inside the total-Cr coverage cells ("0.97, of which 62% unbounded"). Table S5d gives width, unbounded share and mean within-study coverage at every level.
- *Practical reliability.* §3.5 and §5.1 now say that, with literature-scale numbers of studies, the intervals are too wide, or unbounded, to serve as practically reliable predictions for treatment design, and that a design decision rests on confirmatory testing of the material under the target conditions (§3.6).

---

## Minor comments

| No. | Comment | Response |
|---|---|---|
| 1 | §4.4 items per class | New paragraph "Which items each class needs": B needs item 1; A1 needs items 1–4, including a final or controlled pH and a reported electrolyte; A2 needs A1 plus a solid Cr inventory and phase or coordination evidence; item 5 sets the analytical error. |
| 2 | Appendix A wording | Theorem "In words" and the scope paragraph revised (Major Comment 2). |
| 3 | §2.8 fixed split | The methods state that the original run drew one calibration subset per fold with the same fixed seed in every fold of a corpus, and that it is reported as a single realisation. |
| 4 | Extra-binomial results | Nominal and exploratory (Major Comment 3). |
| 5 | Unbounded proportions | Placed in the same cell as the total-Cr coverage (Major Comment 7). |
| 6 | Fig. 6 | Caption corrected to match the panels: (a) fraction of studies with ≥ 90% coverage, (b) per-study coverage, (c) mean width; M2s role-averaged. |
| 7 | DOI | Still requires the author's accounts; see below. |
| 8 | Response-letter typo | The script name was printed without its underscore because a code span inside bold text was rendered with a subscript. The correct name is `scripts/m2s_exchangeability.py`; this letter no longer places code inside bold or italic text. |
| 9 | Length | The assembled manuscript is about 22,000 words including tables and appendices. We have not shortened it further in this round, because most of this revision's changes are statements of scope that the reviewer requested; we are willing to move §3.7–§3.8 detail and Appendix B into the Supplement if the editor prefers. |

---

## Reproducibility

- *Tests.* The suite has 42 tests; `release/test_log.txt` records {NTESTS} passed, line coverage {COV}.
- *Interpreter.* We could not install CPython 3.13.17 in our build environment. Rather than claim the locked interpreter, we updated `requirements-lock.txt`: it now records that the 0.3.3–0.3.5 tests and analyses ran on CPython 3.13.18, with every package version unchanged. A final run on 3.13.17 is listed as an author action.
- *Source archive.* The source archive of code commit `{COMMIT_SHORT}` (`git archive`, repository root) is supplied with this revision; it contains `cassia/release/MANIFEST.sha256` and runs `bash scripts/reproduce_all.sh` without git metadata.
- *Regeneration.* {FULLRUN}
- *Release and DOI.* The tagged release and Zenodo DOI still require the author's accounts.

---

## Items the author must complete before submission

We have not filled these in, because they require the author's own information, judgement or accounts.

| Item | Where |
|---|---|
| Surname, affiliation, ORCID | Title page (`build_docx.js`), `CITATION.cff`, `LICENSE` |
| Funding statement | Back matter (currently marked [Author to confirm]) |
| CRediT statement and any co-authors | Back matter |
| The author's own statement of verification of code and results | §2.12 (placeholder) |
| Generative-AI declaration, checked against the target journal's policy | Back matter |
| GitHub release `cassia-v0.3.5` and Zenodo DOI | Data availability, `CITATION.cff` (steps in `README.md`) |
| Final test run in the exact locked environment (Python 3.13.17), or acceptance of the 3.13.18 lock note | `release/test_log.txt`, `requirements-lock.txt` |
| Full end-to-end rerun of `reproduce_all.sh` (9–11 h) | Local machine |
