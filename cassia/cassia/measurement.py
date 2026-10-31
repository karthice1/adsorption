"""Solid-phase chromium measurement models (CASSIA v0.3).

Oxidation-state evidence and chromium inventory are different analytical tasks and are
modelled separately:

* XPS (surface-sensitive, top few nm): the reported Cr(III) fraction is biased by
  (i) surface enrichment of Cr(III) phases relative to the particle bulk (log-normal
  factor e) and (ii) beam/vacuum-induced reduction of Cr(VI) during analysis (fraction b of
  the Cr(VI) converted), plus fitting noise.
* XANES (bulk-sensitive in fluorescence mode; linear-combination fitting): no surface
  enrichment, smaller possible beam-induced reduction, smaller fitting noise.
* Solid Cr inventory (digestion + total-Cr analysis): total Cr per gram of solid, with
  relative error; it says nothing about oxidation state.
* Phase/coordination evidence (e.g. EXAFS or XANES pre-edge/LCF with references): fraction of
  solid Cr(III) present as precipitated Cr(OH)3 rather than surface-bound, with absolute error.

All parameter values are assumptions of the benchmark (Table S1), not instrument
specifications; they encode the qualitative behaviour reported for these techniques.
"""
from __future__ import annotations

import numpy as np

TECHNIQUES = {
    "XPS": {"enrich_logsd": 0.4, "beam_max": 0.15, "noise": 0.05, "depth": "surface (top few nm)"},
    "XANES": {"enrich_logsd": 0.0, "beam_max": 0.03, "noise": 0.03, "depth": "bulk (fluorescence mode)"},
}
INVENTORY_REL_SD = 0.05
PHASE_ABS_SD = 0.10


def oxidation_state_bias(f_true, enrich=1.0, beam=0.0):
    """Expected reported Cr(III) fraction for a technique with enrichment e and beam reduction b."""
    f = np.clip(np.asarray(f_true, float), 0, 1)
    f_s = enrich * f / (enrich * f + (1 - f)) if np.ndim(f) == 0 else np.where(f > 0, enrich * f / (enrich * f + (1 - f)), 0.0)
    return f_s + beam * (1 - f_s)


def draw_technique_bias(rng, technique):
    t = TECHNIQUES[technique]
    e = float(np.exp(rng.normal(0, t["enrich_logsd"]))) if t["enrich_logsd"] > 0 else 1.0
    b = float(rng.uniform(0, t["beam_max"]))
    return e, b


def measure_oxidation_state(rng, f_true, technique, enrich=None, beam=None):
    """Reported solid Cr(III)/(Cr(III)+Cr(VI)); enrich and beam are drawn per study if not given."""
    t = TECHNIQUES[technique]
    if enrich is None or beam is None:
        enrich, beam = draw_technique_bias(rng, technique)
    mu = oxidation_state_bias(f_true, enrich, beam)
    return np.clip(mu + rng.normal(0, t["noise"], np.shape(mu)), 0, 1)


def measure_inventory(rng, solid_mgg):
    return np.maximum(np.asarray(solid_mgg, float) * (1 + rng.normal(0, INVENTORY_REL_SD, np.shape(solid_mgg))), 0)


def measure_phase(rng, precip_frac_of_cr3):
    return np.clip(np.asarray(precip_frac_of_cr3, float) + rng.normal(0, PHASE_ABS_SD, np.shape(precip_frac_of_cr3)), 0, 1)
