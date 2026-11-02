# Response to the Reviewer (fifth round)

**Manuscript:** CASSIA: A Provenance-Aware Computational Framework and Known-Truth Benchmark for Evidence-Constrained Analysis of Cr(VI) Removal by Carbonaceous Materials (revision 5; software CASSIA 0.3.3, code commit `816b7830e1d2901917bda5f261c3433d4fe5b10e`)

We thank the reviewer for checking the response against the source code itself. The reviewer found a real gap between what we claimed and what we implemented. The rejection threshold pooled individual control studies. The finite-sample rank guarantee we cited holds only if the corpus, not the study, is the exchangeable unit. We accept this criticism in full and have done both things the reviewer proposed.
- **Option 1.** We implemented a corpus-level split-conformal calibration. It shows that a distribution-free 5% threshold is unbounded with 18 corpora.
- **Option 2.** We withdrew the finite-sample guarantee for the threshold actually used. Its leave-one-corpus-out results are now reported as empirical average performance under the simulator, with the 15–25% maximum per-corpus rate shown prominently.

We also explain the D1 results of 14.2% and 10.8%. The other points are addressed as listed below. Several submission items need the author's personal details or accounts and are listed in the last section. Section, table and figure numbers refer to the revised manuscript.

**Summary of changes.**

| Reviewer point | Change | Where |
|---|---|---|
| 3.4 / Priority 1: type-I guarantee | Corpus-level split conformal implemented (unbounded at 5% with 18 corpora; guarantee at 1/19 = 5.3% with power cost). Study-pooled threshold kept as the primary analysis but reported as empirical average cross-corpus performance; finite-sample guarantee withdrawn. D1 14.2%/10.8% explained. Between-corpus heterogeneity and calibration-set sensitivity analysed. | §2.6, §3.3, Table 5, §4.7, §5.1, Appendix B, Table S14b, `scripts/typeI_corpus.py` |
| 5 / Priority 2: Table S5b intervals | Corpus-cluster bootstrap intervals added; Wilson intervals labelled descriptive; the wider envelope is reported | Table S5b, §3.5, `scripts/m2s_family.py` |
| 3.8 / Priority 3: reproducibility | `reproduce_all.sh` no longer needs git (argument, `CASSIA_COMMIT`, git, or `release/RELEASE.md`), checked from a `git archive`; README updated to 0.3.3 | `scripts/reproduce_all.sh`, `README.md`, Data availability |
| 3.7 / Priority 5: Highlights and Abstract | Type-I highlight replaced with the reviewer's wording; Abstract states the per-corpus 15–25% and that the 5% is an empirical average; treatment-design safety labelled simulated | Highlights, Abstract |
| 7 / Priority 4: submission metadata | Cannot be completed by us without the author's details and accounts; listed below | Title page, back matter, `CITATION.cff` |

---

## 3.4 (Priority 1) – Type-I error calibration

**Response: accepted. We implemented the corpus-level calibration (Option 1) and withdrew the finite-sample guarantee for the threshold used (Option 2). The reported procedure is Option 2.**

**The reviewer is right about the implementation.** `conformal_threshold()` in `scripts/typeI_validation.py` takes the ⌈0.95(n + 1)⌉-th smallest of the pooled study-level statistics. The rank lemma then guarantees P(χ²_new > t̂) ≤ 0.05 only for a new study that is exchangeable with those n studies. Studies within a corpus share redrawn family hyper-priors, so a study from a new corpus does not satisfy this premise by construction. Leave-one-corpus-out evaluation and the cluster bootstrap measure the consequence; they do not restore the premise. The revised manuscript no longer states or implies a finite-sample 5% guarantee for this threshold (§2.6, §3.3, Appendix B, §4.7 item 5, §5.1 item 3, Abstract, Highlights).

**Option 1: corpus-level split conformal (new; Appendix B, Table S14b).**
- *Construction.* One control is drawn at random from each of the K calibration corpora. The threshold is the ⌈(1 − α)(K + 1)⌉-th smallest of these K scores (`cassia.validation.cluster_conformal_threshold`). This is the single-subsample construction of Dunn et al. (2023), already cited, with the corpus as the group.
- *Guarantee.* If corpora are i.i.d. draws from the simulator, P(score of a random study in a new corpus > threshold) ≤ α. Appendix B gives a short proof. The guarantee is marginal over corpora and does not bound the rate within any one corpus.
- *Why corpora are exchangeable here.* All 18 control corpora come from the same generator. Material properties are drawn before the scenario switches are applied, and every control is then set to scenario A. They differ only in their seeds.
- *Result, reported as the reviewer anticipated.* With K = 18 corpora the 95% threshold is **unbounded for every design and noise level**, because a finite threshold needs K ≥ 19. The smallest attainable level is 1/19 = 5.3%, where the threshold is the largest of the 18 draws.
- *Behaviour at that level.* Evaluated leave-one-corpus-out (K = 17, level 1/18 = 5.6%), its type-I error was 5.3–5.7%, consistent with the guarantee; the largest per-corpus rate was 12–19%.
- *Cost.* The threshold is unstable: for D3–D5 at 3% noise its 5th–95th percentiles over random draws span 47–9,752. Power for high-reduction studies falls from 88–98% to 75–85% under D3–D5 at 3% noise, and from 80–94% to 59–88% at 10% noise.
- *Conclusion.* With the number of corpora available, a corpus-level finite-sample guarantee exists only at 5.3%, with materially lower power. We report this limitation rather than working around it.

**Option 2: empirical cross-corpus performance (primary analysis).**
- *Primary procedure.* The study-pooled threshold from all 360 controls remains the threshold used for power and equivalence. Its leave-one-corpus-out type-I error (4.7–5.6% on average; corpus-bootstrap 95% intervals within 2.5–8.6%) is described as *empirical average performance under the specified simulator*. It is not described as a guarantee, and not as applying to a new real corpus.
- *Maximum per-corpus rate.* The 15–25% maximum per-corpus rate is now shown wherever the average is quoted: in the Abstract, a Highlight, the Table 5 footnote (in bold), §3.3, §4.7 and §5.1. One to two of the 18 corpora reached ≥ 15% in every design.
- *Context for that maximum.* With 20 controls per corpus, and even if every control had an independent 5% rejection probability, the largest of 18 per-corpus rates would reach 15% with probability 0.76, 20% with probability 0.25 and 25% with probability 0.045.
- *Extra-binomial variation.* It was detected in only 2 of 12 design × noise cells (p = 0.017 and 0.030).
- *No corpus variance component in the null statistic.* The between-corpus share of the variance of log χ² was indistinguishable from zero in all cells (intraclass correlation ≤ 0.008; F-test p ≥ 0.30).
- *How we use these analyses.* They help interpret the spread. They do not turn the average into a guarantee, and the manuscript says so.

**Correction of an earlier statement.** The previous revision asserted that "the null χ² distribution differs between corpora" and that "6 corpora are too few to estimate its 95th percentile". The heterogeneity analysis does not support the first statement as an observed fact. The corpus remains the unit that can be *assumed* exchangeable by design, but no between-corpus variance component was detected. Both sentences have been replaced.

**The D1 results of 14.2% and 10.8% (§3.3, Appendix B, Table S14b).**
- *Facts.* With the conformal order statistic of the same 120 calibration controls, D1 had a type-I error of 14.2% (Wilson 95% interval 10.3–19.1%) and 10.8% (7.5–15.4%) on the 240 independent controls. D2–D5 had 1.7–5.4% with either the empirical or the conformal threshold.
- *Why the conformal threshold does not prevent this.* Its guarantee requires the independent controls to be exchangeable with the 120 calibration controls. They come from 12 different corpora, so that is not guaranteed. In addition, the guarantee holds on average over calibration sets, not for a particular one.
- *What happened.* At 3% noise the largest D1 statistic among the 120 calibration controls was 23.6, whereas 6 of the 240 independent controls exceeded it (maximum 49.8). Refitting the two largest with 20 instead of 3 optimiser starts reproduced their χ² exactly, so they are not optimisation failures. They are genuine upper-tail values that the calibration cohort happened not to contain.
- *How unusual this is.* Over 2,000 random allocations of 6 of the 18 corpora to calibration, a D1 type-I error ≥ 14.2% occurred in 0.1% of allocations (2.5–97.5% range 0.8–12.1%), and ≥ 10.8% at 10% noise in 0.7%. Random allocation of 120 of the 360 controls, ignoring corpora, gave similar probabilities (0.1% and 2%). For at least one of the 12 design × noise cells, the probability was about 1%.
- *Interpretation in the manuscript.* The observed calibration set was an unfavourable draw for D1. It is not explained by a detectable corpus effect. We did not identify any further cause and say so. The two noise levels share the calibration studies and the standardised noise draws, so the two D1 values are not independent confirmations.

**Recommendation stated in the manuscript (§3.3, "Reporting choice"; §4.7 item 5).** We report the study-pooled threshold with its empirical leave-one-corpus-out average and maximum per-corpus rate as the primary analysis. The corpus-level threshold is reported as the procedure that carries a guarantee, at 5.3%, not 5%. Neither is claimed to give 5% type-I control for a new independent corpus. Both would need recalibration, for example on real negative controls, before use on literature data.

---

## 5 (Priority 2) – Confidence intervals in Table S5b

**Response: accepted.** `scripts/m2s_family.py` now adds a corpus-cluster percentile bootstrap: the 10 corpora are resampled with replacement, 4,000 times. Table S5b reports three intervals:
- the cluster interval;
- the Wilson interval, labelled descriptive under the study-level independence approximation;
- the wider envelope of the two, which is the interval quoted in the text.

We report the envelope because a percentile bootstrap over only 10 clusters can itself be too narrow, and it collapses to a single point when every study of a family covers.
- *Hydrochar.* The cluster interval was 76–90% and the Wilson interval 72–92%, so the reported interval is 72–92%.
- *Within-corpus dependence.* In these data it did not widen the intervals. With 10 independent corpora, the family comparison is still described as a tendency, not a demonstrated effect (§3.5).

---

## 3.8 and 4 (Priority 3) – Source reproducibility

**Problem A: `git rev-parse HEAD` in an archive. Corrected.**
- *Commit resolution.* `scripts/reproduce_all.sh` now takes the commit cited in the manuscript, in this order, from:
  1. its first argument;
  2. the environment variable `CASSIA_COMMIT`;
  3. `git rev-parse HEAD`;
  4. the commit documented in `release/RELEASE.md`.
- *Fail early.* The commit is resolved at the start, and the script stops before any computation if none is found, so it can no longer fail after hours of calculation.
- *Manifest.* `scripts/make_release_manifest.py` uses the same rule. Without git or `CASSIA_COMMIT` it records the commit as unknown instead of failing.
- *Check.* We tested this from a `git archive` of the release commit, which has no `.git`:
  - with no argument, the script found the commit in `release/RELEASE.md`;
  - with an argument and with `CASSIA_COMMIT`, it used that value;
  - with `release/RELEASE.md` removed and nothing else given, it stopped before the tests ran;
  - the final assembly and manifest commands then completed in the archive.
- *Full rerun.* The long numerical stages (about 9–11 h on 4 cores) were not rerun from the archive for this check; apart from the two scripts changed in this revision (scripts/typeI_corpus.py and scripts/m2s_family.py, both run in full), they are unchanged since the previous release.

**Problem B: outdated README. Corrected.** `README.md` now refers to CASSIA 0.3.3. It has a reproduction section with the archive usage, a list of the scripts with their tables, the manifest check and a 0.3.3 release checklist; the old v0.2.0 checklist was removed.

**Test execution in the reviewer's runtime.** We cannot change the reviewer's environment. Test collection needs phreeqpython (PHREEQC), which is a declared dependency (`pyproject.toml`, `requirements-lock.txt`), so the tests cannot run without it. With `pip install -e .` at the repository root and `pip install -e cassia[ml,pores,dev]`, the 41 tests pass, now including a test of the corpus-level conformal threshold. Line coverage is 82% (`release/test_log.txt`, `results/coverage_report.txt`). These were run in our environment with Python 3.13.18, a patch release later than the locked 3.13.17; we agree that a supplied log is not an independently reproduced result.

**Tagged release and DOI.** These need the author's GitHub and Zenodo accounts (see the last section).

---

## 3.7 and Priority 5 – Abstract and Highlights

**Response: accepted.**
- *Highlight.* The type-I highlight now uses the reviewer's wording: "Average type-I error was about 5%, with higher rates in some individual corpora". Power is stated separately: "Total-Cr and solid-phase data gave 80–98% rejection power in simulated studies". To stay within five highlights, the first two highlights were merged.
- *Abstract.* It now states that the type-I error averaged about 5% (4.7–5.6%) over held-out corpora, with 15–25% in some individual corpora, and that this is an empirical average, not a guarantee.
- *Treatment design (3.6).* The treatment-design statement now reads "no unsafe accepted design was observed in the simulated cohorts". The corresponding highlight already reads "in the simulated test cohorts".

---

## Points the reviewer marked as resolved (3.1, 3.2, 3.3, 3.5, 3.6)

No further changes, apart from the cluster intervals added to Table S5b (point 3.3) and the word "simulated" added in the Abstract (point 3.6). Table S1 is regenerated from the implementation and is identical except for the version label.

## 6 – Remaining scientific limitations

We agree with the reviewer's list and keep every item in §4.7. The type-I limitation (item 5) has been rewritten as described above. We also agree that the paper's most valuable conclusions, about what measurements can establish regarding adsorption versus reduction, are properties of the assumed simulator and measurement models, and the manuscript says so (§3.3, §4.5, §4.7).

## 7 (Priority 4) – Submission metadata: items the author must supply

The following items cannot be completed without information or actions that only the author can provide. We have not filled them in or invented them. The placeholders are kept visible so that none is overlooked.

| Item | Where | Status |
|---|---|---|
| Author surname, affiliation, ORCID | Title page (`build_docx.js`), `CITATION.cff`, `LICENSE` | Placeholder; awaiting author |
| Funding statement | Back matter | "No specific grant" text marked [Author to confirm] |
| CRediT statement, co-authors | Back matter | Awaiting author confirmation |
| Author's own verification of code and results | §2.12 | Instruction kept as a placeholder; it must be written by the author and cannot be drafted on their behalf |
| Generative-AI declaration | Back matter | Wording drafted; must be confirmed against the target journal's policy |
| GitHub release `cassia-v0.3.3`, Zenodo DOI | Data availability, `CITATION.cff` | Requires the author's GitHub and Zenodo accounts; steps listed in `README.md` |

---

## Other changes

- *Supplement tables.* Column headers containing `|` broke the Word table parser in Table S14; fixed for all supplementary tables.
- *New table.* Table S14b (two parts) is new, and Table 2 lists `scripts/typeI_corpus.py`.
- *No new simulation.* All new analyses reuse the existing 360 negative controls and 120 paired studies; no new simulations were run and no previously reported number changed except where stated above.
- *Status of all results.* All results remain synthetic or thermodynamic; no empirical data were analysed.
