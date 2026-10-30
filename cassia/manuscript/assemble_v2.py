"""Assemble CASSIA_manuscript.md (v2) from the section drafts. Run from manuscript/."""
import re

front = open("_v2_front_methods.md").read()
resA = open("_v2_results_A.md").read()
resML = open("_v2_results_ML.md").read()
resB = open("_v2_results_B.md").read()
disc = open("_v2_discussion.md").read()
tail = open("_v2_back_matter.md").read()

HL = ["Residual Cr(VI) data seldom rule out adsorption-only fits when reduction dominates",
      "Total-Cr and XPS data make reduction detectable, defining evidence classes B and A",
      "Three PHREEQC databases differ by 63–83 mV in the Cr(VI)/Cr(III) redox boundary",
      "Study-level conformal intervals cover 88–93% of new studies but span ≥2 log units",
      "Model-guided, test-confirmed design avoided all unsafe effluent recommendations"]
assert all(len(h) <= 85 for h in HL), [len(h) for h in HL]

AB = ("Loss of dissolved Cr(VI) on carbonaceous materials is routinely reported as adsorption, although Cr(VI) is also reduced to Cr(III). "
      "CASSIA is an open, provenance-aware framework that keeps apparent removal, total-aqueous-Cr removal and adsorption apart, assigns observation-level "
      "evidence classes, and asks which mechanisms a given measurement design can identify. Three PHREEQC databases agreed on Cr(VI) speciation "
      "(within 0.044 in species fraction) but differed by 63–83 mV in the Cr(VI)/Cr(III) boundary; Cr(III) remained undersaturated up to 5,000 mg L⁻¹ "
      "at pH ≤ 3.5–4.0. In known-truth simulations (five mechanism scenarios, 15,000 virtual studies), isotherm fits of residual Cr(VI) with R² ≥ 0.95 "
      "occurred in 41–52% of studies with more than 20% reduction under surface, adsorption-coupled and combined reduction, and Langmuir q_max then "
      "exceeded the adsorption capacity in proportion to the reduction share. A flexible adsorption-only model could not be rejected for 79% of "
      "high-reduction studies observed through residual Cr(VI) alone; adding total Cr and XPS reduced this to 32% and 12.5%. Across 20 replicate "
      "corpora, random cross-validation overstated transferability (random-forest R² 0.89 versus 0.63 leave-one-study-out), speciation features did "
      "not improve prediction, and 45% of residual variance was between-study. Study-level CV+ intervals covered at least 90% of observations in "
      "88–93% of unseen studies but were too wide to certify low effluent concentrations; a model-guided, test-confirmed design protocol found "
      "compliant designs for 54–62% of attainable materials with no unsafe recommendation. No empirical literature data were analysed; the "
      "extraction protocol and pipeline are released for that step.")

ML_DISC = open("_v2_ml_discussion.md").read().strip()
CONCL = {
    "4. [ML CONCLUSION]": "4. Random cross-validation overstated transferability: random-forest R² was 0.89 under random splits, 0.63 under leave-one-study-out, 0.56 under leave-one-family-out and 0.56 on external corpora. Speciation features and nested tuning did not improve prediction. Forty-five per cent of out-of-study variance was between-study, and leave-one-study-out R² rose with the number of independent studies (0.47 at 12, 0.75 at 96).",
    "5. [UQ CONCLUSION]": "5. Pooled-residual intervals reached nominal *average* coverage but covered only 65–75% of new studies at 90%. Study-level CV+ intervals with one score per study covered 88–93% of studies, at widths of 2.1–3.1 log₁₀ units. 95% intervals for total-Cr endpoints require at least 19 calibration studies with total-Cr data.",
    "6. [DESIGN CONCLUSION]": "6. Point predictions gave unsafe treatment designs for 21% (S1) and 7% (S2) of new materials, while calibrated bounds certified almost nothing. A model-guided, test-confirmed protocol found compliant designs for 54–62% of attainable materials with no unsafe recommendation, at a median dose twice the true minimum.",
}

front = front.replace("<!--HIGHLIGHTS-->", "\n".join("- " + h for h in HL)).replace("<!--ABSTRACT-->", AB)
disc = disc.replace("[ML DISCUSSION]", ML_DISC).replace("RASPA3 was audited and not used.", "RASPA3 (Ran et al., 2024) was audited and not used.")
for k, v in CONCL.items():
    disc = disc.replace(k, v)
ms = "\n\n".join([front.rstrip(), resA.strip(), resML.strip(), resB.strip(), disc.strip(), tail.strip()]) + "\n"
open("CASSIA_manuscript.md", "w").write(ms)
left = re.findall(r"\[(?:ML|UQ|DESIGN|PROFILE)[A-Z ]*\]|<!--[A-Z]+-->", ms)
print("abstract words:", len(AB.split()), "| manuscript words:", len(ms.split()), "| unfilled:", left)
