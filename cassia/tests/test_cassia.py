import numpy as np
import pandas as pd
import pytest

from cassia import chemistry as ch, dataset as ds, fitting as fit, identifiability as idf, simulate as sim
from cassia import validation as val


# ----------------------------------------------------------------- chemistry
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
    assert ch.speciate_cr6(9.0, 500.0).fractions["CrO4-2"] > 0.9


def test_redox_scenario_monotone_in_pe():
    f = [ch.redox_scenario(4.0, 50.0, pe)["frac_CrVI"] for pe in (8, 12, 14, 16)]
    assert all(np.diff(f) >= -1e-12) and f[0] < 1e-3 and f[-1] > 0.99


def test_cr3_solubility_by_direct_si():
    r3 = ch.cr3_solubility_si(3.0)
    assert r3["status"] == "undersaturated_to_cmax" and r3["SI_at_cmax"] < 0
    r6 = ch.cr3_solubility_si(6.0)
    assert r6["status"] == "bracketed" and abs(r6["SI_at_S"]) < 0.02
    s = [ch.cr3_solubility_si(p)["S_mgL"] for p in (5.0, 6.0, 7.0)]
    assert s[0] > s[1] > s[2]


def test_llnl_has_no_cr3_hydroxide():
    v, ok = ch.si_cr3(6.0, 1.0, database="llnl.dat")
    assert not ok


def test_speciation_grid_matches_phreeqc_at_node_and_interior():
    g = ch.grid()
    direct = ch.speciate_cr6(3.0, 100.0, ionic_strength=0.01, temperature_C=25.0).fractions
    assert g(3.0, 100.0, 0.01, 298.15)[0][0] == pytest.approx(direct["HCrO4-"], abs=1e-6)
    assert g.validate(n=20, seed=3)["max_abs_err"] < 0.03


def test_database_provenance_checksums():
    p = ch.database_provenance()
    assert set(p.database) == set(ch.DATABASES) and p.sha256.str.len().eq(64).all()


def test_unvalidated_database_rejected():
    with pytest.raises(ValueError):
        ch.engine("phreeqc.dat")


# ----------------------------------------------------------------- fitting
def test_langmuir_recovery():
    rng = np.random.default_rng(0)
    Ce = np.geomspace(0.5, 200, 10)
    q = 50 * 0.1 * Ce / (1 + 0.1 * Ce) * (1 + rng.normal(0, 0.01, Ce.size))
    fits = fit.fit_isotherms(Ce, q)
    k, qm = fits["Langmuir"].params
    assert qm == pytest.approx(50, rel=0.05) and k == pytest.approx(0.1, rel=0.15)
    assert fit.selection_table(fits)["Akaike_w"].sum() == pytest.approx(1.0)


# ----------------------------------------------------------------- dataset
def test_endpoints_and_observation_level_classes():
    df = pd.DataFrame({"study_id": ["s1", "s1"], "C0": [100, 100], "Ce_VI": [20, 20], "Ce_total": [30, np.nan],
                       "D": [2, 2], "xps_cr3_frac": [0.4, np.nan], "pH_final": [3.1, np.nan], "I": [0.01, 0.01]})
    e = ds.endpoints(df)
    assert e["q_app_VI"].tolist() == [40, 40] and e.loc[0, "q_tot_aq"] == 35
    c = ds.classify(df)
    assert c["evidence_class"].tolist() == ["A1", "C"]
    assert c["study_best_class"].tolist() == ["A1", "A1"]
    assert c.loc[1, "endpoint_label"] == "apparent_CrVI_disappearance"


def test_validate_flags_celsius_and_inconsistency():
    df = ds.template()
    df.loc[0, ["C0", "Ce_VI", "Ce_total", "D", "t", "T"]] = [50, 10, 5, 1, 60, 25]
    df = df.astype({c: float for c in ["C0", "Ce_VI", "Ce_total", "D", "t", "T"]})
    probs = ds.validate(df)
    assert any("Celsius" in p for p in probs) and any("Ce_VI > Ce_total" in p for p in probs)


# ----------------------------------------------------------------- simulator
@pytest.mark.parametrize("scen", sim.SCENARIOS)
def test_mass_and_electron_balance(scen):
    rng = np.random.default_rng(1)
    hyper = sim.corpus_hyperpriors(rng)
    m = sim.random_material(rng, "biochar_wood", 0, hyper, scen)
    D, C0 = 2.0, 100.0
    tr = sim.run_batch(m, C0, 2.5, D, [30.0, 1440.0], T=303.15, I=0.01)
    total = tr["Ce_VI_true"] + D * tr["q_adsVI_true"] + tr["Cr3_surface_true"] + tr["Cr3_aq_true"] + tr["Cr3_precip_true"]
    assert np.allclose(total, C0, rtol=1e-5)
    electrons = D * (m.R0 - tr["Rs"]) + (tr["Lrel"] - tr["L"])
    assert np.allclose(3 * tr["reduced_mgL"] / ch.CR_MW, electrons, atol=1e-6)
    if scen == "A":
        assert np.all(tr["reduced_mgL"] < 1e-9)


def test_family_reductance_not_fixed_by_construction():
    tops = set()
    for seed in range(12):
        h = sim.corpus_hyperpriors(np.random.default_rng(seed))
        tops.add(max(h, key=lambda f: h[f]["mu_R"]))
    assert len(tops) >= 3


def test_parameter_table_complete():
    t = sim.parameter_table()
    assert set(sim.PARAMETERS) <= set(t.parameter)


# ----------------------------------------------------------------- identifiability
def _one_study(scen="E", seed=4):
    _, mats = sim.virtual_corpus(2, seed=seed, scenario=scen)
    meta = mats.iloc[0].to_dict()
    return sim.Material(**{k: meta[k] for k in sim.Material.__dataclass_fields__}), meta


def test_posterior_sd_not_larger_than_prior_and_more_data_helps():
    m, meta = _one_study()
    fs = idf.fisher_study(m, meta)
    r1, r4 = idf.posterior_share_sd(fs, "D1"), idf.posterior_share_sd(fs, "D4")
    assert r1["post_sd_share"] <= r1["prior_sd_share"] + 1e-9
    assert r4["post_sd_share"] <= r1["post_sd_share"] + 1e-9


def test_equivalence_negative_control_passes():
    m, meta = _one_study()
    r = idf.equivalence_test(sim.apply_scenario(m, "A"), meta, "D3", sd_rel=0.10, seed=1, n_starts=2)
    assert r["share_red"] < 1e-9 and r["p_adsorption_only"] > 0.01


# ----------------------------------------------------------------- validation / conformal
def _toy(n_st=12, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for s in range(n_st):
        off = rng.normal(0, 0.2)
        for _ in range(10):
            x = rng.uniform(0, 1)
            rows.append({"study_id": f"s{s}", "material_family": f"f{s % 3}", "x": x,
                         "q": 10 ** (1 + x + off + rng.normal(0, 0.05)) - 0.1})
    return pd.DataFrame(rows)


def test_study_cvplus_intervals_and_coverage():
    d = _toy()
    cp = val.StudyCVPlus("RF", ["x"], K=5).fit(d[d.study_id != "s0"], "q")
    test = d[d.study_id == "s0"]
    for meth in ("M0", "M1", "M2"):
        lo, hi, c = cp.interval(test, 0.2, meth)
        assert np.all(lo <= c + 1e-9) and np.all(c <= hi + 1e-9)


def test_m2_is_infinite_with_too_few_studies():
    d = _toy(n_st=5)
    cp = val.StudyCVPlus("RF", ["x"], K=4).fit(d[d.study_id != "s0"], "q")
    lo, hi, _ = cp.interval(d[d.study_id == "s0"], 0.05, "M2")
    assert np.isinf(hi).all()


def test_target_censoring_and_smearing():
    assert val.to_target(np.array([-5.0]))[0] == pytest.approx(np.log10(0.1))
    assert val.duan_smear(np.zeros(10)) == pytest.approx(1.0)


# ----------------------------------------------------------------- pores
def test_pygaps_bet_adapter_recovers_analytic_area():
    from cassia.pores import bet_from_published_isotherm
    p = np.linspace(0.01, 0.3, 25)
    n = 10 * 200 * p / ((1 - p) * (1 - p + 200 * p))  # n_m = 10 mmol/g, C = 200
    out = bet_from_published_isotherm(p, n)
    assert out["SBET"] == pytest.approx(10e-3 * 6.02214076e23 * 0.162e-18, rel=1e-3)


def test_evaluate_conformal_and_variance_on_toy():
    d = _toy(n_st=10)
    d["q_app_VI"] = d["q"]
    res, pl, dd, qp = val.evaluate(d, "RF", ["x"], "LOSO")
    assert res["n_studies"] == 10 and res["R2_lo"] <= res["R2"] + 1e-9
    assert len(val.per_study(dd, dd.q_app_VI.values, qp)) == 10
    icc = val.between_study_variance(dd, pl)
    assert 0 <= icc["ICC"] <= 1
    cf = val.conformal_loso(d, target="q_app_VI", features=["x"], K=4, alphas=(0.2,))
    assert set(cf.method) == {"M0", "M1", "M2", "M2s"} and cf.coverage.between(0, 1).all()


def test_chemistry_tables_and_features_flag_missing_I():
    t = ch.speciation_table([3.0, 7.0], cr_mgL=10.0, databases=("minteq.v4.dat",))
    assert len(t) == 2
    assert {"HCrO4-", "Cr2O7-2"} <= set(ch.logk_table().index)
    df = pd.DataFrame({"pH": [3.0, 3.0], "C0": [50.0, 50.0], "I": [0.01, np.nan], "T": [298.15, 298.15]})
    f = ch.speciation_features(df)
    assert f.I_imputed.tolist() == [False, True] and f.f_HCrO4.iloc[0] == pytest.approx(f.f_HCrO4.iloc[1])


def test_kinetics_and_bootstrap():
    t = np.geomspace(1, 600, 9)
    q = 20 * (1 - np.exp(-0.02 * t))
    tab = fit.selection_table(fit.fit_kinetics(t, q))
    assert tab.loc[0, "model"] == "PFO"
    Ce = np.geomspace(0.5, 200, 8)
    fl = fit.fit_isotherms(Ce, 30 * 0.2 * Ce / (1 + 0.2 * Ce) * (1 + np.random.default_rng(1).normal(0, 0.02, 8)))["Langmuir"]
    ci, draws = fit.bootstrap_ci(fl, Ce, B=50)
    assert ci[0, 1] <= fl.params[1] <= ci[1, 1]


def test_material_descriptors_provenance():
    from cassia.pores import material_descriptors
    p = np.linspace(0.01, 0.3, 25)
    n = 10 * 200 * p / ((1 - p) * (1 - p + 200 * p))
    out = material_descriptors({"SBET": 1000.0}, isotherm=(p, n))
    assert out["SBET_source"] == "pyGAPS_from_isotherm" and abs(out["SBET_rel_diff"]) < 0.05


def test_adapted_halfwidth_shrinks_with_offset_removal():
    d = _toy(n_st=12)
    cp = val.StudyCVPlus("RF", ["x"], K=6).fit(d, "q")
    lo, hi, _ = cp.interval(d.iloc[:5], 0.1, "M1")
    hw_m1 = float(np.median((hi - lo) / 2))
    hw_k4 = cp.adapted_halfwidth(4, 0.1)
    assert np.isfinite(hw_k4) and hw_k4 < hw_m1



def test_a2_requires_inventory_and_phase_and_solid_without_total_is_c():
    base = {"C0": 100, "Ce_VI": 20, "D": 2, "Ce_total": 30, "xanes_cr3_frac": 0.4, "pH_final": 3.0, "I": 0.01}
    assert ds.evidence_class(pd.Series(base)) == "A1"
    assert ds.evidence_class(pd.Series(base | {"solid_cr_mgg": 35.0, "cr3_precip_frac": 0.1})) == "A2"
    assert ds.evidence_class(pd.Series({k: v for k, v in base.items() if k != "Ce_total"})) == "C"


def test_measurement_models_xps_biased_xanes_less():
    from cassia import measurement as ms
    rng = np.random.default_rng(0)
    xps0 = np.mean([ms.measure_oxidation_state(rng, 0.0, "XPS") for _ in range(400)])
    xan0 = np.mean([ms.measure_oxidation_state(rng, 0.0, "XANES") for _ in range(400)])
    assert xps0 > xan0 > 0          # beam reduction creates apparent Cr(III) even without reduction
    assert ms.oxidation_state_bias(0.5, enrich=2.0, beam=0.0) == pytest.approx(2 / 3)


def test_split_conformal_guarantee_on_toy():
    hits, n = 0, 0
    for seed in range(6):
        d = _toy(n_st=24, seed=seed)
        test = d[d.study_id == "s0"]; train = d[d.study_id != "s0"]
        sp = val.StudySplitConformal("RF", ["x"], frac_cal=0.5, beta=0.9, seed=seed).fit(train, "q")
        lo, hi, _ = sp.interval(test, 0.2)
        y = val.to_target(test.q.values)
        hits += ((y >= lo) & (y <= hi)).mean() >= 0.9
        n += 1
    assert hits / n >= 0.5 and np.isfinite(hi).all()


def test_equivalence_v3_control_and_reduction():
    m, meta = _one_study()
    ctrl = idf.equivalence_test_v3(sim.apply_scenario(m, "A"), meta, "D3", sd_rel=0.10, seed=1, n_starts=2)
    assert ctrl["lambda_true"] < ctrl["lambda0_0.1"]          # control is truly within the margin
    assert set(k for k in ctrl if k.startswith("equiv_")) == {"equiv_0.05", "equiv_0.1", "equiv_0.2"}
