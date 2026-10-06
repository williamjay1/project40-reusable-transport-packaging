"""Independent continuous fixed-policy joint transport-budget verification.

Only independent_joint_travel_verification.json is written in --results-dir.
No transport/integration producer or polygon moments are imported. Three
already-independent v30 routines supply raw source constraints, cash-ledger
reaggregation, and vertical-slice scipy.integrate.quad. Thresholds use a
separate normalized Brent inversion, with explicit strict constant suprema.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, linprog

from verify_joint_integration_v30 import (
    rebuild_financial_components, slice_probability, source_constraints,
)


TARGET = .95
PROBABILITY_TOL = 1e-9
COEFFICIENT_TOL = 2e-10
GAMMA_TOL = 2e-9
NEUTRAL_TOL = 1e-6
CHARGES = (0., .005)


def rows(path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def fingerprint(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


class Audit:
    def __init__(self):
        self.checks = defaultdict(int)
        self.maximum_errors = defaultdict(float)
        self.failures = []
        self.integrals = self.intervals = 0
        self.quad_error = 0.

    def check(self, condition, group, label):
        self.checks[group] += 1
        if not condition:
            self.failures.append(dict(group=group, location=label))

    def close(self, actual, expected, group, label, atol, rtol=0.):
        error = abs(actual - expected)
        self.maximum_errors[group] = max(self.maximum_errors[group], error)
        self.check(math.isfinite(actual) and math.isfinite(expected)
                   and error <= atol + rtol * max(abs(actual), abs(expected)),
                   group, f"{label}: actual={actual:.17g}, expected={expected:.17g}")

    def probability(self, planes):
        result = slice_probability(planes)
        self.integrals += 1
        self.intervals += result.intervals
        self.quad_error = max(self.quad_error, result.error_bound)
        self.check(result.error_bound < PROBABILITY_TOL / 10,
                   "quadrature_error", "summed QUADPACK error estimate")
        return result.value


@dataclass(frozen=True)
class Inverse:
    gamma: float | None
    kind: str
    probability_at_zero: float


def continuous_supremum(probability):
    """Supremum of a nonincreasing continuous coverage, gamma >= 0.

    A failed zero test is absent, not zero. Brent operates on gamma/hi,
    independent of the producer's bisection. Exact .95 plateaus require their
    upper endpoint rather than Brent's arbitrary zero of the residual.
    Strict constant cost jumps are handled separately by joint_supremum.
    """
    p0 = probability(0.)
    if p0 < TARGET:
        return Inverse(None, "absent_at_zero", p0)
    hi, previous = 1., p0
    for _ in range(40):
        current = probability(hi)
        require(current <= previous + 1e-12, "Increasing probability in doubling bracket")
        if current < TARGET:
            break
        previous, hi = current, hi * 2
    else:
        raise ValueError("Doubling bracket exhausted; an infinite budget is NOT inferred")

    normalized = brentq(lambda z: probability(z * hi) - TARGET,
                         0., 1., xtol=5e-16, rtol=4e-15, maxiter=200)
    root = normalized * hi
    probe = min(hi, root + max(1e-7, abs(root) * 1e-6))
    if probe > root and probability(probe) >= TARGET:
        # The residual has a plateau at zero: explicitly retain its supremum.
        lo = root
        for _ in range(80):
            mid = lo + (hi - lo) / 2
            if mid == lo or mid == hi:
                break
            if probability(mid) >= TARGET:
                lo = mid
            else:
                hi = mid
        root = lo
        kind = "continuous_plateau_supremum"
    else:
        kind = "continuous_95_percent_root"
    return Inverse(root, kind, p0)


def joint_supremum(probability, environment_inverse, cost_at_zero, slope, e, h):
    """A constant strict saving plane may fail AT its finite supremum."""
    require(slope <= 0., "Positive differential transport slope cannot use monotone inversion")
    p0 = probability(0.)
    if p0 < TARGET:
        return Inverse(None, "absent_at_zero", p0)
    require(environment_inverse.gamma is not None,
            "Joint zero coverage passed while environmental zero coverage failed")
    if e == h == 0. and slope < 0.:
        cost_cut = -cost_at_zero / slope
        if cost_cut <= environment_inverse.gamma:
            return Inverse(cost_cut, "strict_constant_cost_supremum", p0)
    return continuous_supremum(probability)


def controls(audit):
    known = {
        "square": ([], 1.),
        "empty": ([(0., 0., -1.)], 0.),
        "constant_strict_tie": ([(0., 0., 0.)], 0.),
        "constant_pass": ([(0., 0., 1.)], 1.),
        "u_below_half": ([(1., 0., .5)], .75),
        "diagonal": ([(1., 1., 1.)], 5 / 6),
        "diagonal_complement": ([(-1., -1., -1.)], 1 / 6),
        "horizontal_lower": ([(0., -1., -.5)], .25),
        "half_square": ([(1., 0., .5), (0., 1., .5)], .75**2),
        "inconsistent_vertical": ([(1., 0., .5), (-1., 0., -.5)], 0.),
        "tiny_scaled_half": ([(1e-20, 0., .5e-20)], .75),
        "crossing_envelope": ([(1., 1., 1.), (-1., 1., 0.)], 5 / 12),
    }
    values = {}
    for name, (planes, expected) in known.items():
        value = audit.probability(planes)
        audit.close(value, expected, "analytical_probability_controls", name, 1e-12)
        values[name] = value
    for scale in (1., 1e-20):
        probability = lambda g: audit.probability([(scale, 0., scale * (1 - g))])
        inverse = continuous_supremum(probability)
        audit.close(inverse.gamma, math.sqrt(.05), "analytical_root_controls",
                    f"u<1-gamma, scale={scale}", 1e-13)
        audit.check(probability(inverse.gamma - 1e-6) > TARGET
                    and probability(inverse.gamma + 1e-6) < TARGET,
                    "analytical_root_controls", "continuous root direction")
    absent = continuous_supremum(lambda g: audit.probability([(1., 0., -.1 - g)]))
    audit.check(absent.gamma is None, "analytical_root_controls", "absence not clipped to zero")
    genuine_zero = continuous_supremum(lambda g: max(0., TARGET - .1 * g))
    audit.check(genuine_zero.gamma == 0., "analytical_root_controls",
                "genuine boundary at zero is retained, unlike an absent budget")
    plateau = continuous_supremum(lambda g: .95 if g <= 2 else max(0., .95 - .1 * (g - 2)))
    audit.close(plateau.gamma, 2., "analytical_root_controls", "upper end of exact .95 plateau", 1e-13)
    constant = lambda g: audit.probability([(0., 0., 1 - g)])
    inverse = joint_supremum(constant, Inverse(10., "control", 1.), 1., -1., 0., 0.)
    audit.close(inverse.gamma, 1., "analytical_root_controls", "constant strict cost supremum", 0.)
    audit.close(constant(inverse.gamma), 0., "analytical_root_controls", "strict tie AT supremum", 0.)
    audit.check(constant(inverse.gamma - 1e-6) == 1. and constant(inverse.gamma + 1e-6) == 0.,
                "analytical_root_controls", "constant-jump direction")
    for charge in CHARGES:
        flat = lambda g: audit.probability([(0., 0., -charge)])
        audit.check(joint_supremum(flat, Inverse(10., "control", 1.),
                                  -charge, 0., 0., 0.).gamma is None,
                    "analytical_root_controls", f"flat strict tie with charge={charge}")
    try:
        joint_supremum(lambda g: 1., Inverse(10., "control", 1.), 1., 1., 0., 0.)
    except ValueError:
        rejected = True
    else:
        rejected = False
    audit.check(rejected, "analytical_root_controls", "positive financial slope rejected")
    try:
        continuous_supremum(lambda g: min(1., .96 + .01 * g))
    except ValueError as error:
        rejected = "Increasing probability" in str(error)
    else:
        rejected = False
    audit.check(rejected, "analytical_root_controls", "increasing coverage rejected before inversion")
    for t, expected in [(-2., -2.), (0., 0.), (2., 6.)]:
        audit.close(max(t * 1., t * 3.), expected, "analytical_transport_controls",
                    f"signed endpoint, t={t}", 0.)
    values.update(analytical_gamma95=math.sqrt(.05), strict_constant_supremum=1.,
                  plateau_supremum=plateau.gamma)
    return values


def context_key(row):
    return (row["case"], float(row["missing_fate_virgin_equivalent_rho"]),
            float(row["electricity_price_context_eur_per_kwh"]),
            float(row["declared_effective_heat_price_eur_per_kwh"]))


def optional_float(row, name):
    value = row[name]
    return float(value) if value != "" else None


def flow_reaggregation(args, outputs, cash, audit):
    """Use independently supplied geometry and actual positive-flow events."""
    meta = read_json(args.input_dir / "geometry" / "geometry_provenance.json")
    geometry = rows(args.input_dir / "geometry" / "geometry_numeric_facts.csv")
    audit.check(meta == read_json(args.financial_dir / "geometry_provenance.json"),
                "geometry", "frozen financial geometry provenance")
    audit.check(geometry == rows(args.financial_dir / "geometry_numeric_facts.csv"),
                "geometry", "financial geometry equals release numeric inputs")
    lookup = {(r["origin"], r["destination"]): float(r["distance_km"]) for r in geometry}
    require(len(geometry) == len(lookup) == 490, "Expected 490 unique geometry distances")
    distance = {(j, i): lookup[origin, destination]
                for j, origin in enumerate(meta["sorters"])
                for i, destination in enumerate(meta["washers"])}
    require(all(math.isfinite(d) and d >= 0. for d in distance.values()), "Invalid geometry")
    accum = defaultdict(lambda: defaultdict(list))
    unique = set()
    flow_rows = rows(args.financial_dir / "calendar_positive_flows.csv")
    for row in flow_rows:
        key = row["scenario"], row["mode"]
        require(key in outputs, f"Unknown positive-flow policy: {key}")
        tag, coord = row["variable"], tuple(map(int, row["coordinates"].split(":")))
        identifier = key, tag, coord
        require(identifier not in unique, f"Duplicate flow event: {identifier}")
        unique.add(identifier)
        value = float(row["value"])
        require(math.isfinite(value) and value >= 0., f"Invalid flow: {identifier}")
        if tag not in ("wash", "dispatch"):
            continue
        result = outputs[key]
        s, rate = result["scenario"], result["annual_discount_rate"]
        km = value * distance[coord[1], coord[2]]
        week = coord[0] + 1
        factor = math.exp(-math.log1p(rate) * week / 52)
        name = "dirty" if tag == "wash" else "clean"
        accum[key][name].append(km)
        accum[key][tag].append(value)
        clean_factor = 1. if tag == "wash" else s["clean_leg_factor"]
        accum[key]["transport_pv"].append(
            km * clean_factor * s["cost_per_truck_km"] / s["declared_empty_crates_per_truck"] * factor)
    require(set(accum) == set(outputs), "Missing policy transport-flow records")
    computed = {}
    for key, result in outputs.items():
        sums = {name: math.fsum(items) for name, items in accum[key].items()}
        counts, s = result["counts"], result["scenario"]
        for field, name in [("empty_dirty_crate_km", "dirty"),
                            ("empty_clean_crate_km", "clean"),
                            ("wash_inputs", "wash"), ("delivered_services", "dispatch")]:
            audit.close(sums[name], counts[field], "flow_counters", f"{key} {field}", 2e-6, 2e-12)
        audit.close(sums["transport_pv"], cash[key]["transport"], "flow_transport_cash_pv",
                    str(key), 1e-7, 2e-12)
        audit.close(sums["dispatch"], math.fsum(math.fsum(d) for d in result["calendar_demands"]),
                    "flow_services", str(key), 2e-7, 2e-12)
        require(s["clean_leg_factor"] == s["transport_burden_factor"] == 1.,
                f"Physical crate-km interpretation requires the frozen unit factors: {key}")
        computed[key] = dict(M=sums["dirty"] + s["clean_leg_factor"] * sums["clean"],
                             W=sums["wash"], N=sums["dispatch"])
    return computed, len(flow_rows)


@dataclass
class Environment:
    records: dict
    M: float
    W: float
    N: float
    audit: Audit

    def planes(self, gamma):
        require(math.isfinite(gamma) and gamma >= 0., "Negative/nonfinite physical gamma")
        t = (gamma * self.M / 139 - self.W) / 100
        return [(r["electricity_slope"], r["heat_slope"],
                 1 - (r["source_support_upper"]
                      + max(t * r["transport_low"], t * r["transport_high"])) / r["denominator"])
                for r in self.records.values()]

    def probability(self, gamma):
        return self.audit.probability(self.planes(gamma))


def build_environment(result, rho, flow, polys, audit):
    counts = result["counts"]
    vector = np.array([(counts["virgin_assets"] + rho * counts["missing_assets"]) / 100,
                       counts["terminal_recovered_assets"] / 100,
                       counts["closed_loop_rejections"] / .55,
                       counts["closed_loop_rejections"] / .55, flow["W"] / 100])
    require(flow["M"] > 0. and flow["W"] > 0. and flow["N"] > 0., "Invalid frozen flow denominators")
    records = {}
    for indicator, (A, b, su_lower, process) in polys.items():
        support = linprog(-vector, A_ub=A, b_ub=b, bounds=[(None, None)] * 5, method="highs")
        require(support.success, f"Source support LP failed for {indicator}: {support.message}")
        tl, th = process["return_transport"]
        require(0. < tl <= th, f"Positive transport-endpoint prerequisite failed: {indicator}")
        audit.check(0. < tl <= th, "monotonicity_prerequisites", indicator)
        denominator = flow["N"] / 100 * su_lower
        electric = max(flow["W"] / 100 * x for x in process["electricity"]) / denominator
        heat = max(flow["W"] / 100 * x for x in process["water_heating"]) / denominator
        require(electric >= 0. and heat >= 0., "Unexpected negative source energy slopes")
        records[indicator] = dict(source_support_upper=-float(support.fun), denominator=denominator,
                                  electricity_slope=electric, heat_slope=heat,
                                  transport_low=tl, transport_high=th,
                                  total_empty_crate_km=flow["M"], wash_inputs=flow["W"])
    return Environment(records, flow["M"], flow["W"], flow["N"], audit)


def verify(args, audit):
    fdir, idir, root = args.financial_dir, args.integration_dir, args.results_dir
    paths = [Path(__file__), Path(__file__).with_name("verify_joint_integration_v30.py"),
             root / "analysis_protocol.md", fdir / "calendar_network_results.json",
             fdir / "calendar_positive_flows.csv", fdir / "discounted_cashflow_ledger.csv",
             fdir / "independent_financial_verification.json",
             fdir / "independent_calendar_verification.json",
             fdir / "geometry_numeric_facts.csv", fdir / "geometry_provenance.json",
             args.input_dir / "geometry" / "geometry_numeric_facts.csv",
             args.input_dir / "geometry" / "geometry_provenance.json",
             args.input_dir / "tua_reported_grid.csv", args.input_dir / "tua_process_contributions.csv",
             args.stage_facts, idir / "integrated_joint_decisions.csv",
             idir / "environmental_halfplanes.csv", idir / "integration_summary.json",
             idir / "independent_integration_verification.json"]
    producer_paths = [root / name for name in ("joint_travel_thresholds.csv",
                      "continuous_environment_coefficients.csv", "joint_travel_summary.json")]
    missing = [str(p) for p in producer_paths if not p.exists()]
    paths += [p for p in producer_paths if p.exists()]
    before = {str(p): fingerprint(p) for p in paths}
    protocol = (root / "analysis_protocol.md").read_text(encoding="utf-8-sig")
    require("Continuous joint transport-budget protocol" in protocol, "Unexpected analysis protocol")
    results = read_json(fdir / "calendar_network_results.json")["outputs"]
    outputs = {(r["financial_case_id"], r["mode"]): r for r in results}
    require(len(results) == len(outputs) == 56, "Expected 56 unique frozen financial policies")
    require(all(r["financial_case_id"] == r["scenario"]["name"] for r in results),
            "Policy/flow identifier mismatch")
    physical = read_json(fdir / "independent_calendar_verification.json")
    audit.check(physical["successful_network_schedules"] == 56
                and all(c["flow_and_asset_balance"] and c["cost_recalculation"]
                        and c["terminal_closure"] and c["impact_envelope_recalculation"]
                        for c in physical["checks"]),
                "frozen_policy_verification", "prior independent physical/source checks")
    audit.check(read_json(idir / "independent_integration_verification.json")["all_passed"] is True,
                "frozen_integration_verification", "v30 independent report")
    network_hash = fingerprint(fdir / "calendar_network_results.json")
    audit.check(read_json(idir / "integration_summary.json")["source_network_sha256"] == network_hash,
                "provenance", "v30 frozen network hash")
    cash = rebuild_financial_components(fdir, outputs, audit)
    flows, flow_count = flow_reaggregation(args, outputs, cash, audit)
    polys = source_constraints(args.input_dir, args.stage_facts)

    decisions = rows(idir / "integrated_joint_decisions.csv")
    require(len(decisions) == 432, "Expected all 432 frozen v30 decisions")
    baseline, contexts = {}, {}
    for row in decisions:
        key, gamma = context_key(row), float(row["travel_factor"])
        require(gamma in (1., 3.) and (key, gamma) not in baseline, "Duplicate/unexpected v30 gamma")
        baseline[key, gamma] = row
        metadata = {name: row[name] for name in ("case", "parent_scenario", "annual_discount_rate",
                    "terminal_net_eur_per_asset", "missing_fate_virgin_equivalent_rho",
                    "electricity_price_context_eur_per_kwh", "declared_effective_heat_price_eur_per_kwh")}
        if key in contexts:
            audit.check(contexts[key] == metadata, "context_alignment", f"{key} gamma-free metadata")
        else:
            contexts[key] = metadata
    require(len(contexts) == 216 and set(baseline) == {(k, g) for k in contexts for g in (1., 3.)},
            "Expected exactly 216 paired gamma-free contexts")
    emitted_v30 = rows(idir / "environmental_halfplanes.csv")
    planes_v30 = {(r["case"], float(r["travel_factor"]),
                   float(r["missing_fate_virgin_equivalent_rho"]), r["indicator"]): r for r in emitted_v30}
    require(len(planes_v30) == len(emitted_v30) == 1008, "Expected 1008 v30 coefficient rows")

    supplied, producer_coefficients, producer_summary = {}, {}, None
    if not missing:
        produced = rows(producer_paths[0])
        for row in produced:
            key = context_key(row), float(row["extra_net_cost_eur_per_discounted_service"])
            require(key not in supplied, f"Duplicate producer threshold: {key}")
            supplied[key] = row
        require(len(produced) == 432 and set(supplied) == {(k, c) for k in contexts for c in CHARGES},
                "Missing/extra joint threshold contexts or charge coordinates")
        coefficient_rows = rows(producer_paths[1])
        producer_coefficients = {(r["case"], float(r["missing_fate_virgin_equivalent_rho"]),
                                  r["indicator"]): r for r in coefficient_rows}
        require(len(producer_coefficients) == len(coefficient_rows) == 504,
                "Expected 504 unique continuous environmental coefficient rows")
        producer_summary = read_json(producer_paths[2])
        audit.check(producer_summary["contexts"] == 216
                    and producer_summary["joint_threshold_rows"] == 432
                    and producer_summary["unique_environmental_vectors"] == 36
                    and producer_summary["charge_coordinates"] == list(CHARGES),
                    "producer_summary", "complete matrix")
        for path, digest in producer_summary["input_sha256"].items():
            audit.check(fingerprint(Path(path)) == digest, "provenance", f"producer frozen input: {path}")

    env_cache, pair_cache = {}, {}
    verified_rows, classifications = [], Counter()
    seen_coefficient_keys = set()
    root_residuals = []
    baseline_differences = []
    neutral_residuals = []
    for key in sorted(contexts):
        case, rho, pe, ph = key
        label = str(key)
        P, D = outputs[case, "pooled"], outputs[case, "dedicated_adaptive"]
        require(pe >= 0. and ph >= 0. and rho >= 0., "Unexpected negative context coordinate")
        for field in ("parent_scenario", "annual_discount_rate", "terminal_net_eur_per_asset"):
            actual = contexts[key][field]
            for policy in (P, D):
                audit.check(actual == str(policy[field]), "context_alignment", f"{label} {field}")
        if case not in pair_cache:
            pc, dc = cash[case, "pooled"], cash[case, "dedicated_adaptive"]
            audit.close(pc["services"], dc["services"], "cash_service_denominator", case, 1e-7, 2e-12)
            delta_cost = dc["cost"] - pc["cost"]
            delta_transport = dc["transport"] - pc["transport"]
            Wgap = pc["washes"] - dc["washes"]
            neutral = all(abs(x) < NEUTRAL_TOL for x in (delta_cost, delta_transport, Wgap))
            if P["parent_scenario"] == "flat_tau1":
                require(neutral, f"Flat components not roundoff-neutral: {case}")
                for name in P["present_components_eur"]:
                    audit.close(D["present_components_eur"][name], P["present_components_eur"][name],
                                "flat_component_neutrality", f"{case} {name}", NEUTRAL_TOL)
                neutral_residuals.append(dict(case=case, delta_cost_eur=delta_cost,
                                               delta_transport_eur=delta_transport, Wgap=Wgap))
            else:
                require(not neutral and delta_transport < 0.,
                        f"Nonflat monotonic inversion prerequisites failed: {case}")
            audit.check(delta_transport <= 0. or neutral, "monotonicity_prerequisites", case)
            if neutral:
                delta_cost = delta_transport = Wgap = 0.
            pair_cache[case] = dict(delta_cost=delta_cost, delta_transport=delta_transport,
                                    Wgap=Wgap, Npv=pc["services"], neutral=neutral)
        pair = pair_cache[case]
        g0 = (pair["delta_cost"] - pair["delta_transport"]) / pair["Npv"]
        slope = pair["delta_transport"] / pair["Npv"]
        e = .0248 * pe * pair["Wgap"] / pair["Npv"]
        h = (.322 / 3.6) * ph * pair["Wgap"] / pair["Npv"]
        env_key = case, rho
        if env_key not in env_cache:
            env = build_environment(P, rho, flows[case, "pooled"], polys, audit)
            for indicator, record in env.records.items():
                emitted_key = case, rho, indicator
                seen_coefficient_keys.add(emitted_key)
                if supplied:
                    require(emitted_key in producer_coefficients, f"Missing coefficient: {emitted_key}")
                    for field, expected in record.items():
                        actual = float(producer_coefficients[emitted_key][field])
                        normalized = field in ("electricity_slope", "heat_slope", "transport_low", "transport_high")
                        audit.close(actual, expected,
                                    "environment_coefficients_normalized" if normalized else "environment_coefficients_raw",
                                    f"{emitted_key} {field}", COEFFICIENT_TOL if normalized else 2e-6,
                                    0. if normalized else 2e-12)
            for gamma in (1., 3.):
                planes = env.planes(gamma)
                for (indicator, record), (a, b, c) in zip(env.records.items(), planes):
                    supplied_plane = planes_v30[case, gamma, rho, indicator]
                    for field, expected in [("a", a), ("b", b), ("c", c), ("base_ratio", 1 - c)]:
                        audit.close(float(supplied_plane[field]), expected, "v30_environment_coefficients",
                                    f"{env_key} gamma={gamma} {indicator} {field}", COEFFICIENT_TOL)
            zero_transport = 139 * env.W / env.M
            probes = sorted({0., .125, .25, .5, 1., 3., 10., 30.,
                             max(0., zero_transport - 1e-7), zero_transport, zero_transport + 1e-7})
            probabilities = [env.probability(g) for g in probes]
            audit.check(all(a >= b - 1e-12 for a, b in zip(probabilities, probabilities[1:])),
                        "arbitrary_gamma_environment_direction", str(env_key))
            at_kink = env.planes(zero_transport)
            for indicator, plane in zip(env.records, at_kink):
                expected = 1 - env.records[indicator]["source_support_upper"] / env.records[indicator]["denominator"]
                audit.close(plane[2], expected, "transport_sign_switch", f"{env_key} {indicator}", 1e-12)
            inverse = continuous_supremum(env.probability)
            if inverse.gamma is not None:
                probability = env.probability(inverse.gamma)
                audit.close(probability, TARGET, "environment_root_probability", str(env_key), PROBABILITY_TOL)
                delta = max(1e-7, inverse.gamma * 1e-6)
                audit.check(env.probability(max(0., inverse.gamma - delta)) >= TARGET - PROBABILITY_TOL
                            and env.probability(inverse.gamma + delta) <= TARGET + PROBABILITY_TOL,
                            "environment_root_direction", str(env_key))
            env_cache[env_key] = env, inverse
        env, env_inverse = env_cache[env_key]
        intensity = env.M / env.N

        def joint(gamma, charge):
            return audit.probability([*env.planes(gamma), (e, h, g0 + slope * gamma - charge)])

        for gamma in (1., 3.):
            old = baseline[key, gamma]
            for computed, field in [(env.probability(gamma), "integrated_environmental_prior_coverage"),
                                    (joint(gamma, 0.), "integrated_joint_prior_coverage_at_zero_extra_cost")]:
                expected = float(old[field])
                audit.close(computed, expected, "v30_gamma_one_three_probabilities",
                            f"{label} gamma={gamma} {field}", PROBABILITY_TOL)
                baseline_differences.append(abs(computed - expected))
            for field, expected in [("cost_base_allowance_eur_per_discounted_service", g0 + slope * gamma),
                                    ("cost_energy_slope_e", e), ("cost_energy_slope_h", h)]:
                audit.close(float(old[field]), expected, "v30_cost_coefficients",
                            f"{label} gamma={gamma} {field}", COEFFICIENT_TOL)

        for charge in CHARGES:
            row_key = key, charge
            row_label = str(row_key)
            fun = lambda gamma: joint(gamma, charge)
            inverse = joint_supremum(fun, env_inverse, g0 - charge, slope, e, h)
            gamma_root = inverse.gamma
            if gamma_root is None:
                status = "strict_cost_neutrality" if pair["neutral"] else "no_nonnegative_joint_budget"
                audit.check(fun(0.) < TARGET, "absent_budget_zero_test", row_label)
                root_probability = below = above = None
            else:
                status = "joint_margin_above_reference" if gamma_root > 1 else "joint_threshold_below_reference"
                require(env_inverse.gamma is not None, "Joint root without environmental root")
                audit.check(gamma_root <= env_inverse.gamma + GAMMA_TOL, "joint_environment_order", row_label)
                delta = max(1e-7, gamma_root * 1e-6)
                below, above = fun(max(0., gamma_root - delta)), fun(gamma_root + delta)
                root_probability = (0. if inverse.kind == "strict_constant_cost_supremum" else fun(gamma_root))
                if inverse.kind != "strict_constant_cost_supremum":
                    root_residuals.append(abs(root_probability - TARGET))
                    audit.close(root_probability, TARGET, "joint_root_probability", row_label, PROBABILITY_TOL)
                audit.check(below >= TARGET - PROBABILITY_TOL and above <= TARGET + PROBABILITY_TOL,
                            "independent_joint_root_direction", row_label)
            classifications[status] += 1
            probe_gammas = sorted({0., .125, .25, .5, 1., 3., 10., 30.,
                                  139 * env.W / env.M, *([] if gamma_root is None else [gamma_root])})
            probe_probabilities = [fun(g) for g in probe_gammas]
            audit.check(all(a >= b - 1e-12 for a, b in zip(probe_probabilities, probe_probabilities[1:])),
                        "arbitrary_gamma_joint_direction", row_label)
            audit.check(all(p <= env.probability(g) + 1e-12 for g, p in zip(probe_gammas, probe_probabilities)),
                        "joint_is_environment_subset", row_label)
            if pair["neutral"]:
                audit.check(g0 == slope == e == h == 0. and gamma_root is None
                            and all(p == 0. for p in probe_probabilities),
                            "flat_strict_global_controls", row_label)
            expected_values = {
                "environmental_gamma95": env_inverse.gamma, "joint_gamma95": gamma_root,
                "empty_crate_km_per_delivered_service": intensity,
                "joint_intensity_ceiling_empty_crate_km_per_service": None if gamma_root is None else gamma_root * intensity,
                "environmental_intensity_ceiling_empty_crate_km_per_service": None if env_inverse.gamma is None else env_inverse.gamma * intensity,
                "environment_coverage_at_zero_travel": env_inverse.probability_at_zero,
                "joint_coverage_at_zero_travel": inverse.probability_at_zero,
                "joint_coverage_at_gamma_one": fun(1.), "joint_coverage_at_gamma_three": fun(3.),
                "cost_saving_at_gamma_zero_eur_per_discounted_service": g0,
                "differential_transport_cost_slope_eur_per_discounted_service": slope,
                "cost_energy_slope_e": e, "cost_energy_slope_h": h,
            }
            supplied_root = None
            if supplied:
                produced = supplied[row_key]
                for field, metadata in contexts[key].items():
                    audit.check(produced[field] == metadata, "producer_context_alignment", f"{row_label} {field}")
                for field, expected in expected_values.items():
                    actual = optional_float(produced, field)
                    audit.check((actual is None) == (expected is None), "producer_numeric_existence", f"{row_label} {field}")
                    if expected is None or actual is None:
                        continue
                    if field.endswith("gamma95"):
                        group, tolerance, rtol = "producer_gamma_roots", GAMMA_TOL, 2e-11
                    elif "coverage" in field:
                        group, tolerance, rtol = "producer_probabilities", PROBABILITY_TOL, 0.
                    elif "crate_km" in field:
                        group, tolerance, rtol = "producer_intensities", 2e-6, 2e-12
                    else:
                        group, tolerance, rtol = "producer_cost_coefficients", COEFFICIENT_TOL, 0.
                    audit.close(actual, expected, group, f"{row_label} {field}", tolerance, rtol)
                audit.check(produced["joint_status"] == status, "producer_classifications", row_label)
                audit.check(produced["numerical_cost_neutrality"] == str(pair["neutral"]),
                            "producer_neutrality", row_label)
                supplied_root = optional_float(produced, "joint_gamma95")
                root_fields = ("joint_coverage_at_root", "joint_coverage_below_root", "joint_coverage_above_root")
                if supplied_root is None:
                    audit.check(all(produced[field] == "" for field in root_fields),
                                "absent_root_fields", row_label)
                    audit.check(produced["joint_intensity_ceiling_empty_crate_km_per_service"] == "",
                                "absent_root_fields", f"{row_label} no clipped intensity")
                else:
                    # These abscissae are the published producer probe convention,
                    # not its root algorithm. Re-integrate at the SUPPLIED root.
                    delta = max(1e-7, abs(supplied_root) * 1e-6)
                    supplied_root_probability = (0. if inverse.kind == "strict_constant_cost_supremum"
                                                  else fun(supplied_root))
                    prob_values = (supplied_root_probability,
                                   fun(max(0., supplied_root - delta)), fun(supplied_root + delta))
                    for field, expected in zip(root_fields, prob_values):
                        audit.close(float(produced[field]), expected, "producer_root_probe_probabilities",
                                    f"{row_label} {field}", PROBABILITY_TOL)
                    if inverse.kind != "strict_constant_cost_supremum":
                        audit.close(supplied_root_probability, TARGET, "supplied_root_probability",
                                    row_label, PROBABILITY_TOL)
                    audit.check(prob_values[1] >= TARGET - PROBABILITY_TOL
                                and prob_values[2] <= TARGET + PROBABILITY_TOL,
                                "supplied_root_direction", row_label)
            verified_rows.append(dict(context=list(key), charge=charge, status=status,
                                      root_kind=inverse.kind, **expected_values,
                                      independent_coverage_at_root=root_probability,
                                      independent_coverage_below_root=below,
                                      independent_coverage_above_root=above,
                                      arbitrary_probe_gammas=probe_gammas,
                                      arbitrary_probe_probabilities=probe_probabilities))

    audit.check(len(env_cache) == 36 and len(pair_cache) == 28 and len(verified_rows) == 432,
                "matrix_completeness", "36 environmental vectors, 28 pairs, 432 joint rows")
    audit.check(sum(r["status"] == "strict_cost_neutrality" for r in verified_rows) == 48,
                "matrix_completeness", "all 48 strict flat controls retained")
    if supplied:
        audit.check(seen_coefficient_keys == set(producer_coefficients), "matrix_completeness", "all 504 coefficients used")
        audit.check(dict(classifications) == producer_summary["classifications"],
                    "producer_summary", "independent classifications")
        audit.close(producer_summary["analytical_root"], math.sqrt(.05),
                    "producer_summary", "analytic monotone root", 1e-12)
    for key in contexts:
        zero, extra = [next(r for r in verified_rows if r["context"] == list(key) and r["charge"] == charge)
                       for charge in CHARGES]
        zero_probes = dict(zip(zero["arbitrary_probe_gammas"], zero["arbitrary_probe_probabilities"]))
        extra_probes = dict(zip(extra["arbitrary_probe_gammas"], extra["arbitrary_probe_probabilities"]))
        shared_gammas = set(zero_probes) & set(extra_probes)
        audit.check(all(extra_probes[g] <= zero_probes[g] + PROBABILITY_TOL for g in shared_gammas),
                    "extra_charge_direction", str(key))
        if extra["joint_gamma95"] is not None:
            audit.check(zero["joint_gamma95"] is not None
                        and extra["joint_gamma95"] <= zero["joint_gamma95"] + GAMMA_TOL,
                        "extra_charge_root_direction", str(key))
    for path in paths:
        audit.check(fingerprint(path) == before[str(path)], "provenance", f"unchanged: {path}")
    return dict(producer_comparison_pending=bool(missing), missing_producer_files=missing,
                financial_policies=56, financial_pairs=28, gamma_free_contexts=216,
                joint_threshold_evaluations=432, continuous_environment_vectors=len(env_cache),
                continuous_environment_coefficients=len(seen_coefficient_keys),
                original_v30_probability_comparisons=len(baseline_differences),
                maximum_v30_gamma_one_three_probability_difference=max(baseline_differences),
                positive_flow_records_reaggregated=flow_count,
                classifications=dict(classifications), flat_component_residuals=neutral_residuals,
                maximum_independent_joint_root_probability_residual=max(root_residuals, default=0.),
                maximum_quadpack_error_estimate=audit.quad_error,
                quadrature_evaluations=audit.integrals, quadrature_intervals=audit.intervals,
                independent_threshold_records=verified_rows, source_sha256=before)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("financial-dir", "input-dir", "stage-facts", "integration-dir", "results-dir"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    audit = Audit()
    report = dict(
        method="Independent vertical-slice scipy.integrate.quad; raw-source LP; normalized Brent roots with strict constant suprema",
        density="4(1-u)(1-h) on [0,1]^2; independent triangular(1,1,2) input multipliers",
        scope="Frozen v29 policies/cash ledgers; declared prior and empty-crate exposure, not empirical EU distance or a new network optimization",
        independence="Only independent v30 source_constraints, rebuild_financial_components and slice_probability imported; no producers or polygon moments",
        probability_absolute_tolerance=PROBABILITY_TOL,
        gamma_absolute_tolerance=GAMMA_TOL,
        coefficient_absolute_tolerance=COEFFICIENT_TOL,
        neutrality_threshold=NEUTRAL_TOL,
        tolerance_justification="All vertical-envelope intersections split cubic marginal integrands; quadrature uses 2e-13 targets. Cost planes use discounted-service normalization. Brent roots use gamma/bracket normalization and independently checked coverage residuals and both directions.",
        limitations=["The independent .95 values are numerical design-prior integrations, not empirical confidence limits.",
                     "Exposure scaling freezes fleet, flows, return deadlines and source vehicle/loading; it is not a feasible geographic expansion certificate.",
                     "The raw proprietary background LCIA and real operator cash costs are not reconstructed.",
                     "Fresh-directory producer cold reproduction belongs to the main agent and is not duplicated here."],
    )
    try:
        report["analytical_controls"] = controls(audit)
        report.update(verify(args, audit))
    except Exception as error:
        audit.failures.append(dict(group="fatal", location=f"{type(error).__name__}: {error}"))
    pending = report.get("producer_comparison_pending", False)
    report.update(all_passed=not audit.failures and not pending,
                  independent_checks_passed=not audit.failures,
                  checks=dict(audit.checks), maximum_absolute_errors=dict(audit.maximum_errors),
                  failures=audit.failures)
    destination = args.results_dir / "independent_joint_travel_verification.json"
    destination.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    compact = {k: v for k, v in report.items() if k not in ("independent_threshold_records", "source_sha256", "analytical_controls")}
    print(json.dumps(compact, indent=2, allow_nan=False))
    return 1 if audit.failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
