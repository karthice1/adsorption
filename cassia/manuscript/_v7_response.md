# Response to the Reviewer (sixth round)

**Manuscript:** CASSIA: A Provenance-Aware Computational Framework and Known-Truth Benchmark for Evidence-Constrained Analysis of Cr(VI) Removal by Carbonaceous Materials (revision 6; software CASSIA 0.3.4, code commit `{COMMIT}`)

We thank the reviewer for another careful report. The main new point is that the M2s guarantee had been justified by a sentence in Appendix A ("independent draws of studies from a common population") that the hierarchical simulator does not satisfy. We agree. In this revision:
- we establish the exchangeability M2s actually uses under the exact simulation design, which is random role assignment within a corpus, and prove it;
- we restrict the guarantee to the population for which that holds;
- we report numerically what happens outside that scope.

The other comments are addressed point by point below. Items that need the author's own details or accounts are listed at the end; we have not filled them in. Section, table and figure numbers refer to the revised manuscript.

**Summary of changes.**

| Comment | Change | Where |
|---|---|---|
| 1 M2s exchangeability | Exchangeability derived from the random assignment of test, proper-training and calibration roles within a corpus, with a short proof; the reviewer's five questions answered explicitly; guarantee restricted to a study drawn at random from the same corpus; family-stratified alternative shown to be unbounded; new numerical check (Table S5c) and unit test | Appendix A, §2.8, §3.5, §4.7 item 7, Table S5c, `scripts/m2s_exchangeability.py` |
| 2 Calibration wording | "Study-pooled, empirically calibrated" and "corpus-level conformal" used consistently; "as guaranteed" replaced by the reviewer's wording | Abstract, Table 1, §2.6, §3.3, §4.1, §5.1, Appendix B |
| 3 Novelty claim | Narrowed to the methodological integration demonstrated (the reviewer's option 2); claim about existing databases removed | §1 Contribution, §4.5 |
| 4 Reproducibility | 42 tests throughout; locked package versions confirmed; source archive supplied | §2.2, Data availability, archive |
| 5 Empirical validation | Unchanged positioning; no values introduced | §3.8, §4.7 |
| Minor 1–4 | Highlight states that coverage is marginal; hydrochar Wilson interval corrected to 71–92%; Appendix A wording; "as guaranteed" | Highlights, §3.5, Table S5b, Appendix A, §3.3 |

---

## Major Comment 1 – Exchangeability for the M2s guarantee under the hierarchical simulator

**Response: accepted.** The previous Appendix A sentence was wrong for our design and has been removed. The exchangeability the theorem needs does hold in the simulation, but for a different reason than we gave: it comes from the random assignment of roles. The guarantee is now stated only for the population this argument covers.

**The design, stated exactly (Appendix A, "Exchangeability under the simulator design").**
- A corpus has 30 studies.
- Material families are assigned in a fixed rotation: study *i* has family *i* mod 6, so there are five studies per family. Families are not randomised.
- Family hyper-priors are drawn once per corpus and are shared by the corpus's studies.
- Materials and study conditions are then drawn independently for each study.
- Consequently, studies of different families are not identically distributed, and studies of one corpus are dependent through the shared hyper-priors.

**Where the exchangeability comes from.**
- *Evaluation procedure.* In the M2s evaluation each study is the test study in turn, with equal weight. The remaining studies are split at random into proper-training and calibration sets, independently of which study is the test study.
- *Consequence.* Over the test index and the split, (test study, proper-training set, calibration set) is a uniformly random partition of the corpus's studies.
- *Proposition and proof.* Given the realised corpus and the proper-training set, the test study is therefore a uniformly random member of the J + 1 non-training studies. Their scores are computed by a model that does not depend on them, so the test score is exchangeable with the J calibration scores. The theorem's assumption thus holds *conditionally on the realised corpus*, whatever the family allocation and the shared random effects. The proof takes three lines and is now in Appendix A.
- *Code check.* We checked `cassia.validation.StudySplitConformal` and `run_v3.py` (stage m2s). The model is fitted on the proper-training studies only. The calibration subset is drawn by a random generator whose draw does not depend on the test index, and every study serves as the test study once.

**Answers to the reviewer's five questions.**
1. *Family assignment.* It is fixed, in rotation, not randomised. The proof does not use it.
2. *Exchangeability given the proper-training set.* Calibration and test studies are exchangeable conditional on the proper-training set, and on the corpus, because of the random role assignment.
3. *Shared family random effects.* They change the distribution of the study scores; for example, one family's scores can be systematically larger. They do not change the scores' exchangeability under random role assignment. What they do change is family-conditional coverage.
4. *Which population.* The guarantee concerns a study drawn at random from the corpus that supplied the training and calibration studies, i.e. from the same realised hierarchical population. It is not a guarantee for a study from a new, independent corpus, nor for a new study whose family is chosen non-randomly. §2.8, §3.5 and §4.7 now say so.
5. *Held-out study from another family.* The guarantee is an average over the corpus's family mix. It does not hold conditionally on the test study's family, including a family that is rare or absent among the calibration studies.

**Numerical check (new Table S5c, `scripts/m2s_exchangeability.py`).**
- *Method.* On the same 10 corpora, model and β = 0.9, we drew 20 random proper-training sets per corpus and took every remaining study as the test study in turn.
- *Within-corpus rate.* It equals J/(J + 1) = 0.933 for the Cr(VI) endpoints in every corpus. We say plainly in Appendix A that this is an identity once all role assignments are enumerated. It confirms that the implementation realises the proposition; it is not independent evidence.
- *Outside the theorem's scope (informative).* For apparent uptake, coverage conditional on the test study's family was:
  - 0.84 for hydrochar and 0.87 for manure biochar, i.e. below 0.90;
  - 0.93–0.99 for the other families.

  With model and calibration from one corpus and test studies from another, the rate averaged 0.93 but fell to 0.85 in the worst corpus pair.

  These results show concretely why the guarantee must not be read as family-conditional or as transferring to a new corpus.

**The other two acceptable approaches.**
- *Hierarchical (Mondrian, family-stratified) conformal.* With five studies per family and corpus, it has at most four calibration studies per family. Its 90% interval is therefore unbounded, because it needs at least 9. This is stated in Appendix A.
- *Restricting the guarantee.* We therefore keep the theorem, prove its assumption for the actual design, and restrict the exact statement to the within-corpus population (your third approach). A unit test reproduces the argument with a fixed family rotation and a large shared family effect: marginal coverage holds, and coverage for the affected family fails.

---

## Major Comment 2 – Study-pooled versus corpus-level calibration

**Response: accepted.**
- *Consistent terminology.* The two procedures are named consistently throughout:
  - the *study-pooled, empirically calibrated* threshold, whose type-I error is evaluated leave-one-corpus-out;
  - the *corpus-level conformal* threshold, whose finite-sample guarantee holds only at 5.3% with 18 corpora.

  This applies in the Abstract, Table 1, §2.6, §3.3, §4.1, the Conclusions (§5.1) and Appendix B. Bare "calibrated test" no longer appears.
- *"As guaranteed".* §3.3 now reads: "the observed type-I error was 5.3–5.7%, consistent with the theoretical 5.6% bound within Monte Carlo uncertainty".
- *Retained.* The 15–25% maximum single-corpus rates, the D1 results of 14.2% and 10.8%, and the statement that both thresholds must be recalibrated, for example on real negative controls, before use on literature data are all retained. The study-pooled threshold is not recommended as a generally reliable 5% test for new literature studies (§4.7 item 5).
- *Related wording elsewhere.* For the same reason, two observed coverage results in §3.5 and §5.1 are now described as "consistent with" the M2s guarantee rather than "met" or "held".

---

## Major Comment 3 – Novelty claim

**Response: accepted; we take the reviewer's preferred option 2.**
- *Removed.* The sentence in §4.5 that existing adsorption databases and ML studies "typically pool residual-Cr(VI) observations without these distinctions" has been removed, because we have not assessed those resources systematically.
- *What §4.5 now says.* The contribution is the specific methodological integration demonstrated: a known-truth benchmark of how analytical measurement design, hierarchical study-to-study variation and uncertainty methods affect the identifiability of Cr(VI) adsorption, reduction and precipitation. §4.5 also states explicitly that we make no claim about how existing datasets, databases or software treat these distinctions, and no claim of better performance on real chromium-removal data.
- *Introduction.* The Contribution paragraph in §1 opens with the same sentence.

---

## Major Comment 4 – Source archive and test information

**A. Test count.** The suite now has 42 tests: one test was added for the corpus-level conformal threshold in revision 5, and one in this revision for the M2s role-assignment argument. The count of 42 is now used in §2.2, Data availability, the README, `reproduce_all.sh` and this letter, and it matches `release/test_log.txt` ({NTESTS} passed; line coverage {COV}). The earlier "40" in §2.2 was an oversight.

**B. Python version.**
- *Packages.* The test log and all reruns used exactly the locked versions of every scientific dependency: NumPy 2.5.3, SciPy 1.18.1, pandas 3.0.5, scikit-learn 1.9.1, XGBoost 3.4.1, phreeqpython 1.6.2 and matplotlib 3.11.2.
- *Interpreter.* Only the interpreter differs, as a patch release: Python 3.13.18 against the locked 3.13.17.
- *Why we could not use 3.13.17.* We tried to install 3.13.17 in our build environment, but no distribution of it was available there.
- *What remains.* We agree the final suite should be run in the exact locked environment, and we list it as an author action before publication.

**C. Source archive.**
- *Supplied.* The source archive of the code commit `{COMMIT_SHORT}` (`git archive`, repository root) is supplied with this revision.
- *Contents.* It contains `cassia/release/MANIFEST.sha256` for verification and runs `bash scripts/reproduce_all.sh` without git metadata.
- *Full rerun.* {FULLRUN}
- *Release and DOI.* The tagged release and the Zenodo DOI still require the author's accounts (see below).

---

## Major Comment 5 – Empirical validation

**Response: agreed; no change of positioning.**
- *No values introduced.* No experimental values have been introduced, transcribed or estimated.
- *Ramirez et al. (2020).* The dataset remains Class C, and its demonstration will be completed when the source files are available. The computational environment still cannot reach the hosts that serve them.
- *Scope.* The manuscript keeps the known-truth benchmark scope throughout and does not present simulated detection rates as experimentally established capabilities.

---

## Minor comments

| No. | Comment | Response |
|---|---|---|
| 1 | Highlight on coverage | Now "Exact study-level coverage is marginal, not per family; intervals are often too wide" (84 characters). |
| 2 | Hydrochar Wilson interval | The Wilson lower limit is 0.7149, so the correct value is 71–92%; the text had wrongly rounded it to 72%. §3.5 now reports the Wilson interval as 71–92% and the reported envelope as 71–92%, matching Table S5b. |
| 3 | Appendix A "independent draws" | Removed and replaced with the exchangeability argument above. |
| 4 | "As guaranteed" | Replaced (Major Comment 2). |
| 5 | Test count | Corrected (Major Comment 4A). |
| 6–9 | DOI, title page, funding and authorship, AI declaration and verification | Require the author; see below. |
| – | Unbounded total-Cr intervals | Agreed. Table 7b reports the unbounded share next to the total-Cr coverage (63% of test studies at 90%), and Table S5b reports it by family (53–66%). §3.5 and §4.7 now state that the total-Cr coverage must always be read together with this share. |

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
| GitHub release `cassia-v0.3.4` and Zenodo DOI | Data availability, `CITATION.cff` (steps in `README.md`) |
| Final test run in the exact locked environment (Python 3.13.17) | `release/test_log.txt` |
