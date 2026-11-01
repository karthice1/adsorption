"""Assemble CASSIA_manuscript.md (v3) from the section drafts. Run from manuscript/."""
import re
import sys

from latex2uni import convert

COMMIT = sys.argv[1] if len(sys.argv) > 1 else "[COMMIT]"

front = open("_v3_front_methods.md").read()
parts = [open(f).read().strip() for f in ("_v3_results_A.md", "_v3_results_ML.md", "_v3_results_B.md", "_v3_discussion.md")]
appendix = open("_v3_appendix.md").read().strip()
tail = open("_v3_back_matter.md").read().strip()

HL = ["Known-truth benchmark of what Cr(VI)-removal data can show about mechanism",
      "Simulated Cr(VI)-only data were often equivalent to adsorption-only fits within 10%",
      "At about 5% type-I error, total-Cr and solid-phase data gave 80–98% rejection power",
      "Exact study-level coverage holds; intervals are too wide for some engineering uses",
      "No unsafe accepted designs were observed in the simulated test cohorts"]
assert all(len(h) <= 85 for h in HL), [len(h) for h in HL]

AB = ("This methodological study uses known-truth simulations to establish what common analyses of Cr(VI) removal by carbonaceous "
      "materials can and cannot recover. Loss of dissolved Cr(VI) is routinely reported as adsorption, although Cr(VI) is also reduced "
      "to Cr(III). CASSIA is an open, provenance-aware framework that separates apparent removal, total-aqueous-Cr removal and "
      "adsorption, assigns observation-level evidence classes (C, B, A1, A2), and tests which mechanisms a measurement design can "
      "identify. Three PHREEQC databases agreed on Cr(VI) speciation within 0.044 but differed by 63–83 mV in the Cr(VI)/Cr(III) "
      "boundary. Across five simulated mechanism scenarios, adsorption-only was not rejected for 65–78% of high-reduction studies "
      "observed through residual Cr(VI) alone, and a noncentral-χ² equivalence test showed it equivalent within a 10% margin for "
      "43–73%. Calibrated on 360 negative controls from 18 corpora, the rejection test had a type-I error of 4.7–5.6% for every "
      "measurement design on average across corpora (95% intervals within 2.5–8.6%); its power was 80–88% with total Cr and 92–98% "
      "with XPS or XANES. Random cross-validation overstated transferability (random-forest R² 0.86–0.89 versus 0.48–0.65 "
      "leave-one-study-out). An exact study-level split-conformal interval achieved its stated coverage guarantee, but the 90% interval "
      "for apparent Cr(VI) uptake spanned approximately 3.3 log₁₀ units; intervals for total-chromium endpoints were frequently "
      "unbounded with literature-scale calibration samples. With sequential confirmatory testing, no unsafe accepted design was observed "
      "(Wilson 95% upper bound 1.3% of materials, 3.5–3.8% of accepted designs). No published experimental data were analysed; the "
      "framework is benchmarked, not experimentally validated.")

front = front.replace("<!--HIGHLIGHTS-->", "\n".join("- " + h for h in HL)).replace("<!--ABSTRACT-->", AB)
parts[-1] = parts[-1].replace("RASPA3.** It was audited", "RASPA3.** RASPA3 (Ran et al., 2024) was audited")
ms = "\n\n".join([front.rstrip(), *parts, tail.split("## References")[0].strip(), appendix, "## References" + tail.split("## References")[1]]) + "\n"
ms = ms.replace("[COMMIT]", COMMIT)
open("CASSIA_manuscript.md", "w").write(ms)
open("_build_manuscript.md", "w").write(ms.replace(appendix, convert(appendix)))
open("_v3_highlights.md", "w").write("# Highlights\n\n" + "\n".join("- " + h for h in HL) + "\n")
left = re.findall(r"\[(?:ML|UQ|DESIGN|PROFILE|COMMIT|DOI)[A-Z ]*\]|<!--[A-Z]+-->", ms)
print("abstract words:", len(AB.split()), "| manuscript words:", len(ms.split()), "| unfilled:", left, "| HL len:", [len(h) for h in HL])
