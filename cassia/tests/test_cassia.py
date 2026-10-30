import numpy as np
import pandas as pd
import pytest

from cassia import chemistry as ch, dataset as ds, fitting as fit, simulate as sim


@pytest.mark.parametrize("db", ch.DATABASES)
def test_cr6_fractions_sum_to_one(db):
    for pH in (2.0, 5.0, 6.5, 9.0):
        f = ch.speciate_cr6(pH, 100.0, database=db).fractions
        assert sum(f.values()) == pytest.approx(1.0, abs=1e-6)


@pytest.mark.parametrize("db", ch.DATABASES)
def test_bichromate_chromate_crossover_near_pKa(db):
    lo = ch.speciate_cr6(4.0, 10.0, database=db).fractions
    hi = ch.speciate_cr6(8.5, 10.0, database=db).fractions
    assert lo["HCrO4-"] > 0.95 and hi["CrO4-2"] > 0.95


def test_dimer_grows_with_concentration():
    a = ch.speciate_cr6(3.0, 10.0).fractions.get("Cr2O7-2", 0)
    b = ch.speciate_cr6(3.0, 500.0).fractions.get("Cr2O7-2", 0)
    assert b > 5 * a


def test_high_concentration_high_pH_converges():
    f = ch.speciate_cr6(9.0, 500.0).fractions
    assert f["CrO4-2"] > 0.9


def test_redox_scenario_monotone_in_pe():
    f = [ch.redox_scenario(4.0, 50.0, pe)["frac_CrVI"] for pe in (8, 12, 14, 16)]
    assert all(np.diff(f) >= -1e-12) and f[0] < 1e-3 and f[-1] > 0.99


def test_cr3_solubility_decreases_from_pH5_to_8():
    s = [ch.cr3_solubility_mgL(p) for p in (5.0, 6.0, 7.0, 8.0)]
    assert all(np.diff(s) < 0)


def test_unvalidated_database_rejected():
    with pytest.raises(ValueError):
        ch.engine("phreeqc.dat")


def test_langmuir_recovery():
    rng = np.random.default_rng(0)
    Ce = np.geomspace(0.5, 200, 10)
    q = 50 * 0.1 * Ce / (1 + 0.1 * Ce) * (1 + rng.normal(0, 0.01, Ce.size))
    fits = fit.fit_isotherms(Ce, q)
    k, qm = fits["Langmuir"].params
    assert qm == pytest.approx(50, rel=0.05) and k == pytest.approx(0.1, rel=0.15)
    tab = fit.selection_table(fits)
    assert tab["Akaike_w"].sum() == pytest.approx(1.0)


def test_endpoints_and_classes():
    df = pd.DataFrame({"C0": [100, 100], "Ce_VI": [20, 20], "Ce_total": [30, np.nan], "D": [2, 2],
                       "solid_phase_cr_speciation": ["XPS", "none"], "pH_basis": ["final", ""],
                       "I": [0.01, np.nan]})
    e = ds.endpoints(df)
    assert e["q_app_VI"].tolist() == [40, 40]
    assert e.loc[0, "q_total"] == 35 and e.loc[0, "reduction_index"] == pytest.approx(10 / 80)
    c = ds.classify(df)
    assert c["evidence_class"].tolist() == ["A", "C"]
    assert c.loc[1, "endpoint_label"] == "apparent_CrVI_disappearance"


def test_validate_flags_celsius_and_inconsistency():
    df = ds.template()
    df.loc[0, ["C0", "Ce_VI", "Ce_total", "D", "t", "T"]] = [50, 10, 5, 1, 60, 25]
    df = df.astype({c: float for c in ["C0", "Ce_VI", "Ce_total", "D", "t", "T"]})
    probs = ds.validate(df)
    assert any("Celsius" in p for p in probs) and any("Ce_VI > Ce_total" in p for p in probs)


def test_simulator_mass_balance():
    rng = np.random.default_rng(1)
    m = sim.random_material(rng, "biochar", 0)
    tr = sim.run_batch(m, 100.0, 2.5, 2.0, 1440.0)
    removed = 100.0 - tr["Ce_VI_true"]
    accounted = tr["q_adsVI_true"] * 2.0 + tr["reduced_mgL"]
    assert removed == pytest.approx(accounted, rel=1e-4)
    assert tr["Cr3_aq_true"] + tr["Cr3_bound_mgL"] + tr["Cr3_precip_mgL"] == pytest.approx(tr["reduced_mgL"], rel=1e-6)


def test_pygaps_bet_adapter_recovers_analytic_area():
    from cassia.pores import bet_from_published_isotherm
    p = np.linspace(0.01, 0.3, 25)
    n = 10 * 200 * p / ((1 - p) * (1 - p + 200 * p))  # n_m = 10 mmol/g, C = 200
    out = bet_from_published_isotherm(p, n)
    assert out["SBET"] == pytest.approx(10e-3 * 6.02214076e23 * 0.162e-18, rel=1e-3)
