"""Assemble CASSIA_manuscript.md (v3) from the section drafts. Run from manuscript/."""
import re

front = open("_v3_front_methods.md").read()
parts = [open(f).read().strip() for f in ("_v3_results_A.md", "_v3_results_ML.md", "_v3_results_B.md", "_v3_discussion.md")]
appendix = open("_v3_appendix.md").read().strip()
tail = open("_v3_back_matter.md").read().strip()

HL = ["Simulated Cr(VI)-only data were often equivalent to adsorption-only fits within 10%",
      "In the benchmark, total-Cr and solid-phase data raised rejection power to 80–96%",
      "PHREEQC databases differ by 63–83 mV in the Cr(VI)/Cr(III) redox boundary",
      "Exact study-level conformal intervals hold coverage but span over 3 log units",
      "No unsafe accepted designs were observed in the simulated test cohorts"]
assert all(len(h) <= 85 for h in HL), [len(h) for h in HL]

AB = ("Loss of dissolved Cr(VI) on carbonaceous materials is routinely reported as adsorption, although Cr(VI) is also reduced to Cr(III). "
      "CASSIA is an open, provenance-aware framework that separates apparent removal, total-aqueous-Cr removal and adsorption, assigns "
      "observation-level evidence classes (C, B, A1, A2), and tests which mechanisms a measurement design can identify. Three PHREEQC "
      "databases agreed on Cr(VI) speciation within 0.044 but differed by 63–83 mV in the Cr(VI)/Cr(III) boundary. In known-truth simulations "
      "with five mechanism scenarios, adsorption-only was not rejected by a calibrated χ² test for 63–76% of high-reduction studies observed "
      "through residual Cr(VI) alone, and a noncentral-χ² equivalence test showed it equivalent within a 10% margin for 43–65%. With total Cr, "
      "calibrated power was 80–88% (cross-validated type-I error 3–5%); with XPS or XANES, 90–96%. Random cross-validation overstated "
      "transferability (random-forest R² 0.86–0.89 versus 0.48–0.65 leave-one-study-out across scenarios), and speciation features did not "
      "improve prediction. An exact study-level split-conformal interval met its proven coverage guarantee but spanned 3.3 log₁₀ units at "
      "90%. Point-prediction treatment designs were unsafe for 7–21% of new virtual materials; with sequential confirmatory testing, no unsafe "
      "accepted design was observed (Wilson 95% upper bound 1.3%). All results are synthetic or thermodynamic: no published experimental data "
      "were analysed, and the extraction protocol and tested pipeline are released for that step.")

front = front.replace("<!--HIGHLIGHTS-->", "\n".join("- " + h for h in HL)).replace("<!--ABSTRACT-->", AB)
parts[-1] = parts[-1].replace("RASPA3.** It was audited", "RASPA3.** RASPA3 (Ran et al., 2024) was audited")
ms = "\n\n".join([front.rstrip(), *parts, tail.split("## References")[0].strip(), appendix, "## References" + tail.split("## References")[1]]) + "\n"
open("CASSIA_manuscript.md", "w").write(ms)
open("_v3_highlights.md", "w").write("# Highlights\n\n" + "\n".join("- " + h for h in HL) + "\n")
left = re.findall(r"\[(?:ML|UQ|DESIGN|PROFILE|COMMIT|DOI)[A-Z ]*\]|<!--[A-Z]+-->", ms)
print("abstract words:", len(AB.split()), "| manuscript words:", len(ms.split()), "| unfilled:", left, "| HL len:", [len(h) for h in HL])
