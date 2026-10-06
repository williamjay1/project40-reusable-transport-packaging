"""Continuous same-policy joint transport budgets under a declared design prior.

No new network optimization or empirical probability is claimed. Raw data and
all frozen financial policies stay unchanged. Independent verifier uses slice
quadrature and a separate root implementation, not this generator.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from calendar_pool_network_v23 import source_polytopes, dump_csv
from integrate_joint_prior_v30 import integrate

TARGET = .95
CONTEXT_FIELDS = ["case", "parent_scenario", "annual_discount_rate",
                  "terminal_net_eur_per_asset", "missing_fate_virgin_equivalent_rho",
                  "electricity_price_context_eur_per_kwh",
                  "declared_effective_heat_price_eur_per_kwh"]


def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def fingerprint(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def context_key(row):
    return tuple(row[k] for k in CONTEXT_FIELDS)


def upper_root(probability):
    """Supremum of a nonincreasing >=TARGET event on gamma>=0."""
    p0 = probability(0.)
    if p0 < TARGET:
        return None
    lo, hi = 0., 1.
    while probability(hi) >= TARGET:
        lo, hi = hi, hi * 2
        if hi > 2**30:
            raise ValueError("No finite root bracket; inspect model before extrapolating")
    for _ in range(70):
        mid = (lo + hi) / 2
        if probability(mid) >= TARGET:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--financial-dir", type=Path, required=True)
    ap.add_argument("--input-dir", type=Path, required=True)
    ap.add_argument("--stage-facts", type=Path, required=True)
    ap.add_argument("--integration-dir", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    if (out / "joint_travel_thresholds.csv").exists():
        raise ValueError("Use a fresh output directory: do not overwrite frozen results")
    network = args.financial_dir / "calendar_network_results.json"
    integration = args.integration_dir / "integrated_joint_decisions.csv"
    tracked = [network, integration, args.input_dir / "tua_reported_grid.csv",
               args.input_dir / "tua_process_contributions.csv", args.stage_facts]
    before = {str(p): fingerprint(p) for p in tracked}
    source = source_polytopes(args.input_dir, args.stage_facts)
    results = json.loads(network.read_text())["outputs"]
    index = {(r["financial_case_id"], r["mode"]): r for r in results}
    assert len(index) == len(results) == 56
    old = read_csv(integration)
    assert len(old) == 432
    contexts = {}
    baselines = {}
    for row in old:
        key = context_key(row)
        contexts.setdefault(key, {k: row[k] for k in CONTEXT_FIELDS})
        baselines[key, float(row["travel_factor"])] = row
    assert len(contexts) == 216 and len(baselines) == 432
    # Analytical root control: P(u<1-gamma)=1-gamma^2 for 0<=gamma<=1.
    analytic = upper_root(lambda g: integrate([(1., 0., 1. - g)]))
    assert abs(analytic - math.sqrt(.05)) < 1e-12
    assert upper_root(lambda g: 0.) is None
    assert abs(upper_root(lambda g: integrate([(0., 0., 1. - g)])) - 1.) < 1e-14
    cache = {}
    emitted_coefficients = []
    rows, baseline_errors = [], []
    for key, context in contexts.items():
        case = context["case"]
        rho = float(context["missing_fate_virgin_equivalent_rho"])
        pe = float(context["electricity_price_context_eur_per_kwh"])
        ph = float(context["declared_effective_heat_price_eur_per_kwh"])
        P = index[case, "pooled"]
        D = index[case, "dedicated_adaptive"]
        counts, scenario = P["counts"], P["scenario"]
        N, W = counts["delivered_services"], counts["wash_inputs"]
        M = counts["empty_dirty_crate_km"] + scenario["clean_leg_factor"] * counts["empty_clean_crate_km"]
        assert M > 0 and N > 0 and W > 0
        environment_key = case, rho
        if environment_key not in cache:
            vector = np.array([(counts["virgin_assets"] + rho * counts["missing_assets"]) / 100,
                               counts["terminal_recovered_assets"] / 100,
                               counts["closed_loop_rejections"] / .55,
                               counts["closed_loop_rejections"] / .55, W / 100])
            records = []
            for indicator, poly in source.items():
                opt = linprog(-vector, A_ub=poly["A"], b_ub=poly["b"],
                              bounds=[(None, None)] * 5, method="highs")
                assert opt.success
                den = N / 100 * poly["u_lower_ratio"]
                transport_low, transport_high = poly["process"]["return_transport"]
                assert 0 < transport_low <= transport_high
                record = dict(indicator=indicator, source_support_upper=float(vector @ opt.x),
                              denominator=den, electricity_slope=W / 100 * poly["process"]["electricity"][1] / den,
                              heat_slope=W / 100 * poly["process"]["water_heating"][1] / den,
                              transport_low=transport_low, transport_high=transport_high,
                              total_empty_crate_km=M, wash_inputs=W)
                assert record["electricity_slope"] >= 0 and record["heat_slope"] >= 0
                records.append(record)
                emitted_coefficients.append(dict(case=case, missing_fate_virgin_equivalent_rho=rho, **record))
            cache[environment_key] = records
        records = cache[environment_key]

        def environmental_planes(gamma):
            t = (gamma * M / 139 - W) / 100
            return [(r["electricity_slope"], r["heat_slope"],
                     1 - (r["source_support_upper"] + max(t * r["transport_low"], t * r["transport_high"])) / r["denominator"])
                    for r in records]

        Npv = P["discounted_service_equivalents"]
        delta_cost = D["included_present_cost_eur"] - P["included_present_cost_eur"]
        delta_transport = D["present_components_eur"]["empty_transport"] - P["present_components_eur"]["empty_transport"]
        delta_wash = P["discounted_wash_equivalents"] - D["discounted_wash_equivalents"]
        neutral = (abs(delta_cost) < 1e-6 and abs(delta_transport) < 1e-6
                   and abs(delta_wash) < 1e-6)
        if case.startswith("flat_"):
            assert neutral, (case, delta_cost, delta_transport, delta_wash)
        else:
            assert not neutral and delta_transport < 0, (case, delta_transport)
        if neutral:
            delta_cost = delta_transport = delta_wash = 0.
        e = .0248 * pe * delta_wash / Npv
        f = (.322 / 3.6) * ph * delta_wash / Npv

        def environment(gamma):
            return integrate(environmental_planes(gamma))

        def joint(gamma, kappa):
            g = (delta_cost + (gamma - 1) * delta_transport) / Npv
            return integrate([*environmental_planes(gamma), (e, f, g - kappa)])

        probes = [0., .25, 1., 3., 10., 30.]
        env_prob = [environment(g) for g in probes]
        assert all(a >= b - 1e-12 for a, b in zip(env_prob, env_prob[1:]))
        envroot = upper_root(environment)
        for gamma in (1., 3.):
            row = baselines[key, gamma]
            err_env = abs(environment(gamma) - float(row["integrated_environmental_prior_coverage"]))
            err_joint = abs(joint(gamma, 0.) - float(row["integrated_joint_prior_coverage_at_zero_extra_cost"]))
            assert max(err_env, err_joint) < 1e-10, (key, gamma, err_env, err_joint)
            baseline_errors.extend([err_env, err_joint])
        for kappa in (0., .005):
            probabilities = [joint(g, kappa) for g in probes]
            assert all(a >= b - 1e-12 for a, b in zip(probabilities, probabilities[1:]))
            root = upper_root(lambda gamma: joint(gamma, kappa))
            if root is not None:
                assert envroot is not None and root <= envroot + 1e-9
                delta = max(1e-7, abs(root) * 1e-6)
                below = joint(max(0., root - delta), kappa)
                above = joint(root + delta, kappa)
                assert below >= TARGET - 1e-9 and above <= TARGET + 1e-9
                status = "joint_margin_above_reference" if root > 1 else "joint_threshold_below_reference"
            else:
                below = above = None
                status = "strict_cost_neutrality" if neutral else "no_nonnegative_joint_budget"
            rows.append(dict(**context, extra_net_cost_eur_per_discounted_service=kappa,
                             environmental_gamma95=envroot, joint_gamma95=root, joint_status=status,
                             empty_crate_km_per_delivered_service=M / N,
                             joint_intensity_ceiling_empty_crate_km_per_service=None if root is None else root * M / N,
                             environmental_intensity_ceiling_empty_crate_km_per_service=None if envroot is None else envroot * M / N,
                             environment_coverage_at_zero_travel=env_prob[0],
                             joint_coverage_at_zero_travel=probabilities[0],
                             joint_coverage_at_gamma_one=probabilities[2],
                             joint_coverage_at_gamma_three=probabilities[3],
                             joint_coverage_at_root=None if root is None else joint(root, kappa),
                             joint_coverage_below_root=below, joint_coverage_above_root=above,
                             cost_saving_at_gamma_zero_eur_per_discounted_service=(delta_cost - delta_transport) / Npv,
                             differential_transport_cost_slope_eur_per_discounted_service=delta_transport / Npv,
                             cost_energy_slope_e=e, cost_energy_slope_h=f,
                             numerical_cost_neutrality=neutral))
    assert len(rows) == 432 and len(cache) == 36
    dump_csv(out / "joint_travel_thresholds.csv", rows)
    dump_csv(out / "continuous_environment_coefficients.csv", emitted_coefficients)
    assert all(fingerprint(Path(p)) == h for p, h in before.items())
    classifications = {status: sum(r["joint_status"] == status for r in rows)
                       for status in sorted({r["joint_status"] for r in rows})}
    report = dict(contexts=len(contexts), charge_coordinates=[0., .005], joint_threshold_rows=len(rows),
                  unique_environmental_vectors=len(cache), analytical_root=analytic,
                  analytical_expected=math.sqrt(.05), classifications=classifications,
                  maximum_gamma_one_three_probability_difference=max(baseline_errors),
                  protected_inputs_unchanged=len(before), input_sha256=before,
                  scope="Fixed terminal-aware policies; declared prior; no new optimization, empirical EU probability or universal distance",
                  independent_verification_pending=True)
    (out / "joint_travel_summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
