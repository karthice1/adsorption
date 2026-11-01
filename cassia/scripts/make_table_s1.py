"""Supplementary Table S1 from the implementation (cassia.simulate.parameter_table()).

    python scripts/make_table_s1.py  -> results/v3/T_S1_simulator_parameters.csv
The table is generated from the same constants the simulator and measurement models use
(EVIDENCE_MIX, cassia.measurement.TECHNIQUES, INVENTORY_REL_SD, PHASE_ABS_SD); a unit test checks this.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import cassia  # noqa: E402
from cassia import simulate as sim  # noqa: E402

t = sim.parameter_table()
t.to_csv(ROOT / "results/v3/T_S1_simulator_parameters.csv", index=False)
print("Table S1:", len(t), "rows; CASSIA", cassia.__version__)
