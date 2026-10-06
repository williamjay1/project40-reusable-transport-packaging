"""Independent rank-dependence and strict support-budget verification.

Negative ranks: analytic angular intersections, not producer z crossings.
Positive ranks: intersection of scalar bounds and the triangular marginal CDF.
Full support: signed, affine corner ceilings for BOTH transport endpoints.
Only independent_joint_dependence_verification.json is written in results-dir.
No producer, producer rank routine, manuscript or package is imported/modified.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

from verify_joint_integration_v30 import rebuild_financial_components, source_constraints
from verify_joint_travel_v31 import (
    Audit, Environment, Inverse, build_environment, context_key,
    continuous_supremum, fingerprint, flow_reaggregation, read_json, require, rows,
)


TARGET = .95
PROBABILITY_TOL = 1e-10
ROOT_TOL = 2e-8
COEFFICIENT_TOL = 2e-10
NEUTRAL_TOL = 1e-6
CHARGES = (0., .005)
MODES = ("independent", "positive_rank", "negative_rank")


@dataclass(frozen=True)
class RankProbability:
    value: float
    intervals: int


def positive_rank_probability(planes):
    """u=h=x; F(x)=2x-x^2. Cancelling rank slopes are strict constants."""
    left, right = 0., 1.
    for a, b, c in planes:
        require(all(math.isfinite(v) for v in (a, b, c)), "Nonfinite rank half-plane")
        slope = math.fsum((a, b))
        if slope == 0.:
            if c <= 0.:
                return RankProbability(0., 0)
            continue
        bound = c / slope
        if slope > 0.:
            right = min(right, bound)
        else:
            left = max(left, bound)
    if right <= left:
        return RankProbability(0., 0)
    # F(right)-F(left), factored to preserve small interval mass.
    return RankProbability((right - left) * (2 - left - right), 1)


def negative_rank_probability(planes):
    """u=1-cos(theta), h=1-sin(theta); density 2 sin(theta) cos(theta).

    a*u+b*h<c is a*cos(theta)+b*sin(theta)>a+b-c. Solve ALL
    theta=atan2(b,a) +/- acos((a+b-c)/hypot(a,b)) + 2*pi*k crossings.
    Strict nonconstant endpoints have zero mass. A midpoint determines each
    interval after the exhaustive cuts; no monotonicity in z is assumed.
    """
    half_pi = math.pi / 2
    cuts = {0., half_pi}
    normalized = []
    for a, b, c in planes:
        require(all(math.isfinite(v) for v in (a, b, c)), "Nonfinite rank half-plane")
        scale = max(abs(a), abs(b))
        if scale == 0.:
            if c <= 0.:
                return RankProbability(0., 0)
            continue
        a, b, c = a / scale, b / scale, c / scale
        d = math.fsum((a, b, -c))
        radius = math.hypot(a, b)
        normalized.append((a, b, d))
        if -radius <= d <= radius:
            angle = math.acos(max(-1., min(1., d / radius)))
            phase = math.atan2(b, a)
            for sign in (-1., 1.):
                for turn in (-1, 0, 1):
                    crossing = phase + sign * angle + turn * 2 * math.pi
                    if 0. < crossing < half_pi:
                        cuts.add(crossing)
    accepted = []
    for left, right in zip(sorted(cuts), sorted(cuts)[1:]):
        midpoint = left + (right - left) / 2
        cosine, sine = math.cos(midpoint), math.sin(midpoint)
        if all(math.fsum((a * cosine, b * sine, -d)) > 0.
               for a, b, d in normalized):
            # sin^2(right)-sin^2(left), avoiding subtraction cancellation.
            accepted.append(math.sin(right + left) * math.sin(right - left))
    probability = math.fsum(accepted)
    require(-1e-13 <= probability <= 1 + 1e-13, "Angular probability outside [0,1]")
    return RankProbability(min(1., max(0., probability)), len(accepted))


def rank_supremum(environment, g0, slope, e, h, charge, mode):
    """Handle strict constants before applying continuous Brent inversion."""
    require(slope <= 0., "Positive cash slope invalidates monotone inversion")
    method = positive_rank_probability if mode == "positive_rank" else negative_rank_probability
    probability = lambda gamma: method([*environment(gamma), (e, h, g0 + slope * gamma - charge)]).value
    p0 = probability(0.)
    if p0 < TARGET:
        return Inverse(None, "absent_at_zero", p0)
    constant = math.fsum((e, h)) == 0. if mode == "positive_rank" else e == h == 0.
    if constant and slope < 0.:
        cut = -(g0 - charge) / slope
        environment_at_cut = method(environment(cut)).value
        if environment_at_cut >= TARGET:
            return Inverse(cut, "strict_constant_supremum", p0)
    return continuous_supremum(probability)


def minimum_corner_margin(planes):
    return min(math.fsum((c, -a * u, -b * h))
               for a, b, c in planes for u in (0., 1.) for h in (0., 1.))


@dataclass(frozen=True)
class Support:
    gamma: float | None
    minimum_margin_at_zero: float
    binding: dict | None


def corner_supremum(env, g0, slope, e, h, charge):
    """Analytic strict guarantee, including the gamma=139*W/M sign switch.

    max(t*low,t*high)<L is equivalent to BOTH affine inequalities.
    For each corner, their ceilings are 139/M*(W+100*L/endpoint).
    A positive L binds the high endpoint; a negative L binds the low one.
    The cash corner remains signed: g0-charge-e*u-h*v+slope*gamma>0.
    """
    require(slope <= 0., "Positive cash slope invalidates support ceiling inversion")
    at_zero = [*env.planes(0.), (e, h, g0 - charge)]
    minimum_zero = minimum_corner_margin(at_zero)
    if minimum_zero <= 0.:
        return Support(None, minimum_zero, None)
    candidates = []
    kink = 139 * env.W / env.M
    for u in (0., 1.):
        for heat in (0., 1.):
            for indicator, record in env.records.items():
                allowance = record["denominator"] * (
                    1 - record["electricity_slope"] * u - record["heat_slope"] * heat
                ) - record["source_support_upper"]
                local = []
                for endpoint_name in ("transport_low", "transport_high"):
                    endpoint = record[endpoint_name]
                    require(endpoint > 0., "Nonpositive source transport endpoint")
                    ceiling = (139 / env.M) * (env.W + 100 * allowance / endpoint)
                    certificate = dict(indicator=indicator, corner=[u, heat],
                                       endpoint=endpoint_name, ceiling=ceiling,
                                       transport_sign_switch_gamma=kink,
                                       transport_allowance=allowance)
                    candidates.append((ceiling, certificate))
                    local.append(ceiling)
                # This is an independent branch identity, not a sign assumption.
                selected = record["transport_high"] if allowance >= 0. else record["transport_low"]
                signed_ceiling = (139 / env.M) * (env.W + 100 * allowance / selected)
                require(abs(min(local) - signed_ceiling) <= 1e-12 * max(1., abs(signed_ceiling)),
                        "Sign-aware transport ceiling identity failed")
                require((signed_ceiling >= kink) == (allowance >= 0.) or abs(signed_ceiling - kink) < 1e-12,
                        "Transport ceiling on incorrect side of the source sign switch")
            cash_zero = math.fsum((g0, -charge, -e * u, -h * heat))
            require(cash_zero > 0., "Strict cash corner failed despite positive minimum at zero")
            if slope < 0.:
                ceiling = cash_zero / (-slope)
                candidates.append((ceiling, dict(indicator="strict_cash_saving", corner=[u, heat],
                                                 endpoint=None, ceiling=ceiling)))
    require(candidates, "No finite corner ceiling; an unbounded guarantee is not silently inferred")
    gamma, binding = min(candidates, key=lambda pair: pair[0])
    require(gamma > 0., "Positive strict zero margins yielded a nonpositive ceiling")
    return Support(gamma, minimum_zero, binding)


def analytic_controls(audit):
    disconnected = 1 - 1.2 * math.sqrt(2 - 1.2**2)
    known = {
        "square": ([], 1., 1.),
        "empty": ([(0., 0., -1.)], 0., 0.),
        "strict_constant_tie": ([(0., 0., 0.)], 0., 0.),
        "positive_constant": ([(0., 0., 1.)], 1., 1.),
        "u_below_half": ([(1., 0., .5)], .75, .75),
        "h_below_half": ([(0., 1., .5)], .75, .75),
        "diagonal": ([(1., 1., 1.)], .75, 1.),
        "mixed_sign": ([(1., -1., 0.)], 0., .5),
        "mixed_sign_reverse": ([(-1., 1., 0.)], 0., .5),
        "disconnected_negative_event": ([(-1., -1., -.8)], .36, disconnected),
        "central_negative_event": ([(1., 1., .8)], .64, 1 - disconnected),
        "intersected_half_square": ([(1., 0., .5), (0., 1., .5)], .75, .5),
        "inconsistent_vertical": ([(1., 0., .5), (-1., 0., -.5)], 0., 0.),
        "positive_cancelled_pass": ([(1., -1., 2.)], 1., 1.),
        "negative_complement_diagonal": ([(-1., -1., -1.)], .25, 0.),
    }
    values = {}
    for scale in (1., 1e-20):
        for name, (planes, expected_positive, expected_negative) in known.items():
            scaled = [(a * scale, b * scale, c * scale) for a, b, c in planes]
            positive = positive_rank_probability(scaled)
            negative = negative_rank_probability(scaled)
            audit.close(positive.value, expected_positive, "analytic_positive_controls", f"{name} scale={scale}", 1e-12)
            audit.close(negative.value, expected_negative, "analytic_negative_controls", f"{name} scale={scale}", 1e-12)
            if name == "disconnected_negative_event":
                audit.check(negative.intervals == 2, "analytic_disconnected_controls", f"scale={scale}")
            values[f"{name}_scale_{scale}"] = dict(positive=positive.value, negative=negative.value,
                                                   negative_intervals=negative.intervals)
        for mode in ("positive_rank", "negative_rank"):
            inverse = rank_supremum(lambda gamma: [], scale, -scale, scale, 0., 0., mode)
            audit.close(inverse.gamma, math.sqrt(.05), "analytic_rank_root_controls", f"{mode} scale={scale}", 1e-13)
            constant = rank_supremum(lambda gamma: [], scale, -scale, 0., 0., 0., mode)
            audit.close(constant.gamma, 1., "analytic_rank_root_controls", f"{mode} strict supremum", 0.)
            method = positive_rank_probability if mode == "positive_rank" else negative_rank_probability
            audit.check(method([(0., 0., scale * (1 - constant.gamma))]).value == 0.,
                        "analytic_rank_root_controls", f"{mode} strict equality fails AT supremum")
            audit.check(rank_supremum(lambda gamma: [], 0., -scale, 0., 0., 0., mode).gamma is None,
                        "analytic_rank_root_controls", f"{mode} absence never clipped to zero")
    cancellation = rank_supremum(lambda gamma: [], 1., -1., 1., -1., 0., "positive_rank")
    audit.check(cancellation.gamma == 1. and cancellation.kind == "strict_constant_supremum",
                "analytic_rank_root_controls", "rank-specific cancelling cash slopes")

    def synthetic(support, electric=0., heat=0.):
        record = dict(source_support_upper=support, denominator=1., electricity_slope=electric,
                      heat_slope=heat, transport_low=1., transport_high=2.)
        return Environment({"synthetic": record}, 13900., 100., 100., audit)

    for support, expected in [(0., 1.5), (1.5, .5)]:
        env = synthetic(support)
        result = corner_supremum(env, 10., 0., 0., 0., 0.)
        audit.close(result.gamma, expected, "analytic_corner_controls", f"transport support={support}", 1e-14)
        audit.close(minimum_corner_margin([*env.planes(result.gamma), (0., 0., 10.)]), 0.,
                    "analytic_corner_controls", "boundary is a strict supremum, not attained", 1e-14)
    audit.check(corner_supremum(synthetic(2.), 10., 0., 0., 0., 0.).gamma is None,
                "analytic_corner_controls", "failed zero margin not clipped")
    audit.close(corner_supremum(synthetic(0., .2, .3), 10., 0., 0., 0., 0.).gamma, 1.25,
                "analytic_corner_controls", "all four environmental corners", 1e-14)
    empty_environment = Environment({}, 13900., 100., 100., audit)
    for scale in (1., 1e-20):
        signed = corner_supremum(empty_environment, 3 * scale, -scale, scale, -2 * scale, 0.)
        audit.close(signed.gamma, 2., "analytic_corner_controls", f"signed cash scale={scale}", 1e-14)
        audit.check(signed.binding["corner"] == [1., 0.], "analytic_corner_controls", "signed cash binding corner")
    audit.check(corner_supremum(empty_environment, 1., -1., 1., 0., 0.).gamma is None,
                "analytic_corner_controls", "one exact strict corner tie invalidates guarantee")
    constant = corner_supremum(empty_environment, 1., -1., 0., 0., 0.)
    audit.close(constant.gamma, 1., "analytic_corner_controls", "strict constant cash supremum", 0.)
    return values


def nullable_number(row, name):
    return float(row[name]) if row[name] != "" else None


def row_key(row):
    return context_key(row), float(row["extra_net_cost_eur_per_discounted_service"])


def status_for(gamma):
    return "absent_at_zero" if gamma is None else ("above_reference" if gamma > 1. else "below_reference")


def compare_optional(audit, row, field, expected, label, group, tolerance, rtol=0.):
    actual = nullable_number(row, field)
    audit.check((actual is None) == (expected is None), "numeric_existence", f"{label} {field}")
    if actual is not None and expected is not None:
        audit.close(actual, expected, group, f"{label} {field}", tolerance, rtol)


def verify(args, audit):
    root, travel, financial = args.results_dir, args.travel_dir, args.financial_dir
    tracked = [Path(__file__), Path(__file__).with_name("verify_joint_travel_v31.py"),
               Path(__file__).with_name("verify_joint_integration_v30.py"), root / "analysis_protocol.md",
               financial / "calendar_network_results.json", financial / "calendar_positive_flows.csv",
               financial / "discounted_cashflow_ledger.csv", financial / "independent_financial_verification.json",
               financial / "independent_calendar_verification.json", financial / "geometry_numeric_facts.csv",
               financial / "geometry_provenance.json", args.input_dir / "geometry" / "geometry_numeric_facts.csv",
               args.input_dir / "geometry" / "geometry_provenance.json", args.input_dir / "tua_reported_grid.csv",
               args.input_dir / "tua_process_contributions.csv", args.stage_facts,
               travel / "joint_travel_thresholds.csv", travel / "continuous_environment_coefficients.csv",
               travel / "independent_joint_travel_verification.json"]
    producer_paths = [root / "joint_dependence_thresholds.csv", root / "joint_dependence_summary.json"]
    missing = [str(path) for path in producer_paths if not path.exists()]
    tracked.extend(path for path in producer_paths if path.exists())
    before = {str(path): fingerprint(path) for path in tracked}
    protocol = (root / "analysis_protocol.md").read_text(encoding="utf-8-sig")
    require("Post-result joint-dependence robustness protocol" in protocol, "Unexpected analysis protocol")
    previous = read_json(travel / "independent_joint_travel_verification.json")
    audit.check(previous.get("all_passed") is True, "prior_verification", "v31 independent verification")
    # Both unchanged live inputs and newly packaged input copies must match the
    # independently frozen v31 fingerprints; no package is edited.
    old_hashes = previous["source_sha256"]
    for path in tracked:
        if str(path) in old_hashes:
            audit.check(before[str(path)] == old_hashes[str(path)], "frozen_v31_hashes", str(path))
    for path in [args.input_dir / "tua_reported_grid.csv", args.input_dir / "tua_process_contributions.csv",
                 args.stage_facts, args.input_dir / "geometry" / "geometry_numeric_facts.csv",
                 args.input_dir / "geometry" / "geometry_provenance.json"]:
        matching = [digest for old_path, digest in old_hashes.items()
                    if Path(old_path).name == path.name]
        require(matching, f"No independently frozen fingerprint for {path.name}")
        audit.check(all(before[str(path)] == digest for digest in matching),
                    "frozen_v31_hashes", f"release copy {path}")

    original_results = read_json(financial / "calendar_network_results.json")["outputs"]
    outputs = {(r["financial_case_id"], r["mode"]): r for r in original_results}
    require(len(outputs) == len(original_results) == 56, "Expected all 56 frozen policies")
    cash = rebuild_financial_components(financial, outputs, audit)
    flows, flow_count = flow_reaggregation(args, outputs, cash, audit)
    polys = source_constraints(args.input_dir, args.stage_facts)
    base_rows = rows(travel / "joint_travel_thresholds.csv")
    baseline = {row_key(r): r for r in base_rows}
    require(len(base_rows) == len(baseline) == 432, "Expected 432 unique v31 context/charge records")
    contexts = {key for key, charge in baseline}
    require(len(contexts) == 216 and set(baseline) == {(key, c) for key in contexts for c in CHARGES},
            "Expected all 216 contexts at both charges")
    coefficient_rows = rows(travel / "continuous_environment_coefficients.csv")
    frozen_coefficients = {(r["case"], float(r["missing_fate_virgin_equivalent_rho"]), r["indicator"]): r
                           for r in coefficient_rows}
    require(len(coefficient_rows) == len(frozen_coefficients) == 504, "Expected 504 frozen environment records")
    produced, summary = {}, None
    if not missing:
        producer_rows = rows(producer_paths[0])
        produced = {row_key(r): r for r in producer_rows}
        require(len(produced) == len(producer_rows) == 432 and set(produced) == set(baseline),
                "Producer missing/extra/duplicate context-charge rows")
        summary = read_json(producer_paths[1])
        audit.check(summary["rows"] == 432, "producer_summary", "row count")
        for path, digest in summary["input_sha256"].items():
            audit.check(fingerprint(Path(path)) == digest, "producer_hashes", path)
        audit.check(summary["new_optimization_solves"] == 0, "producer_summary", "no new policy optimization")

    environment_cache, pair_cache = {}, {}
    independent_records = []
    classifications = {mode: Counter() for mode in (*MODES, "full_support")}
    shifts = {mode: [] for mode in ("positive_rank", "negative_rank", "full_support")}
    root_residuals = {mode: [] for mode in MODES}
    used_coefficients = set()
    angular_evaluations = angular_intervals = 0
    charge_records = defaultdict(dict)
    for key in sorted(contexts):
        case, rho, pe, ph = key
        P, D = outputs[case, "pooled"], outputs[case, "dedicated_adaptive"]
        require(rho >= 0. and pe >= 0. and ph >= 0., "Unexpected negative context coordinate")
        if case not in pair_cache:
            pc, dc = cash[case, "pooled"], cash[case, "dedicated_adaptive"]
            audit.close(pc["services"], dc["services"], "cash_common_services", case, 1e-7, 2e-12)
            delta_cost = dc["cost"] - pc["cost"]
            delta_transport = dc["transport"] - pc["transport"]
            Wgap = pc["washes"] - dc["washes"]
            neutral = all(abs(value) < NEUTRAL_TOL for value in (delta_cost, delta_transport, Wgap))
            if P["parent_scenario"] == "flat_tau1":
                require(neutral, f"Flat case not independently cost-neutral: {case}")
                for component in P["present_components_eur"]:
                    audit.close(D["present_components_eur"][component], P["present_components_eur"][component],
                                "flat_component_neutrality", f"{case} {component}", NEUTRAL_TOL)
            else:
                require(not neutral and delta_transport < 0., f"Nonflat cash monotonicity failed: {case}")
            if neutral:
                delta_cost = delta_transport = Wgap = 0.
            pair_cache[case] = dict(g0=(delta_cost - delta_transport) / pc["services"],
                                   slope=delta_transport / pc["services"], Wgap=Wgap,
                                   Npv=pc["services"], neutral=neutral)
        pair = pair_cache[case]
        g0, slope = pair["g0"], pair["slope"]
        e = .0248 * pe * pair["Wgap"] / pair["Npv"]
        h = (.322 / 3.6) * ph * pair["Wgap"] / pair["Npv"]
        env_key = case, rho
        if env_key not in environment_cache:
            env = build_environment(P, rho, flows[case, "pooled"], polys, audit)
            for indicator, record in env.records.items():
                coefficient_key = case, rho, indicator
                used_coefficients.add(coefficient_key)
                require(coefficient_key in frozen_coefficients, f"Missing frozen coefficient: {coefficient_key}")
                frozen = frozen_coefficients[coefficient_key]
                for field, expected in record.items():
                    normalized = field in ("electricity_slope", "heat_slope", "transport_low", "transport_high")
                    audit.close(float(frozen[field]), expected,
                                "source_coefficients_normalized" if normalized else "source_coefficients_raw",
                                f"{coefficient_key} {field}", COEFFICIENT_TOL if normalized else 2e-6,
                                0. if normalized else 2e-12)
            environment_cache[env_key] = env
        env = environment_cache[env_key]
        intensity = env.M / env.N

        for charge in CHARGES:
            label = str((key, charge))
            base = baseline[key, charge]
            for field in ("parent_scenario", "annual_discount_rate", "terminal_net_eur_per_asset"):
                for policy in (P, D):
                    audit.check(base[field] == str(policy[field]), "context_metadata", f"{label} {field}")
            for field, expected in [("cost_saving_at_gamma_zero_eur_per_discounted_service", g0),
                                    ("differential_transport_cost_slope_eur_per_discounted_service", slope),
                                    ("cost_energy_slope_e", e), ("cost_energy_slope_h", h)]:
                audit.close(float(base[field]), expected, "independent_cash_coefficients", f"{label} {field}", COEFFICIENT_TOL)
            audit.check(base["numerical_cost_neutrality"] == str(pair["neutral"]), "cost_neutrality", label)

            def planes(gamma):
                return [*env.planes(gamma), (e, h, g0 + slope * gamma - charge)]

            def probability(mode, gamma):
                nonlocal angular_evaluations, angular_intervals
                if mode == "independent":
                    return audit.probability(planes(gamma))
                if mode == "positive_rank":
                    return positive_rank_probability(planes(gamma)).value
                answer = negative_rank_probability(planes(gamma))
                angular_evaluations += 1
                angular_intervals += answer.intervals
                return answer.value

            inverses = {}
            inverses["independent"] = continuous_supremum(lambda gamma: probability("independent", gamma))
            for mode in ("positive_rank", "negative_rank"):
                inverses[mode] = rank_supremum(env.planes, g0, slope, e, h, charge, mode)
            support = corner_supremum(env, g0, slope, e, h, charge)
            supplied = produced.get((key, charge))
            audit.close(float(base["empty_crate_km_per_delivered_service"]), intensity,
                        "intensity_base", label, 1e-9, 2e-12)
            base_root = nullable_number(base, "joint_gamma95")
            compare_optional(audit, base, "joint_gamma95", inverses["independent"].gamma,
                             label, "v31_independent_roots", ROOT_TOL)
            base_status = ("strict_cost_neutrality" if pair["neutral"] else "no_nonnegative_joint_budget")
            if inverses["independent"].gamma is not None:
                base_status = "joint_margin_above_reference" if inverses["independent"].gamma > 1 else "joint_threshold_below_reference"
            audit.check(base_status == base["joint_status"], "v31_independent_status", label)
            for gamma, field in [(0., "joint_coverage_at_zero_travel"), (1., "joint_coverage_at_gamma_one"),
                                 (3., "joint_coverage_at_gamma_three")]:
                audit.close(probability("independent", gamma), float(base[field]),
                            "v31_independent_probabilities", f"{label} {field}", PROBABILITY_TOL)
            if supplied:
                metadata_fields = ("case", "parent_scenario", "annual_discount_rate", "terminal_net_eur_per_asset",
                                   "missing_fate_virgin_equivalent_rho", "electricity_price_context_eur_per_kwh",
                                   "declared_effective_heat_price_eur_per_kwh", "extra_net_cost_eur_per_discounted_service")
                for field in metadata_fields:
                    audit.check(supplied[field] == base[field], "producer_contexts", f"{label} {field}")
                compare_optional(audit, supplied, "independent_joint_gamma95", inverses["independent"].gamma,
                                 label, "producer_independent_roots", ROOT_TOL)
                audit.check(supplied["independent_joint_status"] == base_status, "producer_independent_status", label)
                audit.close(float(supplied["empty_crate_km_per_delivered_service"]), intensity,
                            "producer_intensities", label, 1e-9, 2e-12)

            record = dict(context=list(key), extra_cost=charge, intensity=intensity, rank_results={})
            for mode in MODES:
                inverse = inverses[mode]
                gamma_root = inverse.gamma
                status = status_for(gamma_root)
                classifications[mode][status] += 1
                if gamma_root is None:
                    audit.check(inverse.probability_at_zero < TARGET, "absent_rank_zero_test", f"{label} {mode}")
                    at_root = below = above = None
                else:
                    at_root = probability(mode, gamma_root)
                    if inverse.kind != "strict_constant_supremum":
                        root_residuals[mode].append(abs(at_root - TARGET))
                        audit.close(at_root, TARGET, "rank_root_probability", f"{label} {mode}", PROBABILITY_TOL)
                    else:
                        audit.close(at_root, 0., "rank_strict_constant_ties", f"{label} {mode}", PROBABILITY_TOL)
                    delta = max(1e-7, gamma_root * 1e-6)
                    below, above = probability(mode, max(0., gamma_root - delta)), probability(mode, gamma_root + delta)
                    audit.check(below >= TARGET - PROBABILITY_TOL and above <= TARGET + PROBABILITY_TOL,
                                "rank_root_direction", f"{label} {mode}")
                probes = sorted({0., .125, .25, .5, 1., 3., 10., 30., 139 * env.W / env.M,
                                 *([] if gamma_root is None else [gamma_root])})
                probe_values = [probability(mode, gamma) for gamma in probes]
                audit.check(all(a >= b - PROBABILITY_TOL for a, b in zip(probe_values, probe_values[1:])),
                            "gamma_monotonicity", f"{label} {mode}")
                reported_values = {"coverage_at_gamma_zero": probability(mode, 0.),
                                   "coverage_at_gamma_one": probability(mode, 1.),
                                   "coverage_at_gamma_three": probability(mode, 3.)}
                if supplied and mode != "independent":
                    compare_optional(audit, supplied, mode + "_joint_gamma95", gamma_root,
                                     label, "producer_rank_roots", ROOT_TOL)
                    audit.check(supplied[mode + "_status"] == status, "producer_rank_status", f"{label} {mode}")
                    compare_optional(audit, supplied, mode + "_intensity_ceiling",
                                     None if gamma_root is None else gamma_root * intensity,
                                     label, "producer_intensities", ROOT_TOL * max(1., intensity))
                    for suffix, expected in reported_values.items():
                        audit.close(float(supplied[mode + "_" + suffix]), expected,
                                    "producer_rank_probabilities", f"{label} {mode} {suffix}", PROBABILITY_TOL)
                    producer_root = nullable_number(supplied, mode + "_joint_gamma95")
                    if producer_root is None:
                        audit.check(supplied[mode + "_coverage_at_supremum"] == "", "absent_supremum_fields", f"{label} {mode}")
                    else:
                        expected = probability(mode, producer_root)
                        audit.close(float(supplied[mode + "_coverage_at_supremum"]), expected,
                                    "producer_supremum_probabilities", f"{label} {mode}", PROBABILITY_TOL)
                        if inverse.kind != "strict_constant_supremum":
                            audit.close(expected, TARGET, "producer_supremum_95_probability", f"{label} {mode}", PROBABILITY_TOL)
                        delta = max(1e-7, producer_root * 1e-6)
                        audit.check(probability(mode, max(0., producer_root - delta)) >= TARGET - PROBABILITY_TOL
                                    and probability(mode, producer_root + delta) <= TARGET + PROBABILITY_TOL,
                                    "producer_rank_root_direction", f"{label} {mode}")
                record["rank_results"][mode] = dict(gamma=gamma_root, status=status, root_kind=inverse.kind,
                                                    intensity_ceiling=None if gamma_root is None else gamma_root * intensity,
                                                    **reported_values, coverage_at_supremum=at_root,
                                                    coverage_below_supremum=below, coverage_above_supremum=above,
                                                    arbitrary_probe_gammas=probes, arbitrary_probe_coverages=probe_values)
                if mode != "independent" and gamma_root is not None and base_root is not None:
                    shifts[mode].append(dict(context=list(key), extra_cost=charge, difference=gamma_root - base_root))
            support_status = status_for(support.gamma)
            classifications["full_support"][support_status] += 1
            margin_one = minimum_corner_margin(planes(1.))
            if support.gamma is not None:
                delta = max(1e-7, support.gamma * 1e-6)
                below = minimum_corner_margin(planes(max(0., support.gamma - delta)))
                above = minimum_corner_margin(planes(support.gamma + delta))
                audit.check(below > 0. and above < 0., "support_root_direction", label)
                audit.close(minimum_corner_margin(planes(support.gamma)), 0., "support_boundary_margin", label, PROBABILITY_TOL)
                for mode, inverse in inverses.items():
                    audit.check(inverse.gamma is not None and support.gamma <= inverse.gamma + ROOT_TOL,
                                "support_threshold_order", f"{label} {mode}")
                if base_root is not None:
                    shifts["full_support"].append(dict(context=list(key), extra_cost=charge, difference=support.gamma - base_root))
            else:
                audit.check(support.minimum_margin_at_zero <= 0., "absent_support_zero_test", label)
                below = above = None
            if supplied:
                compare_optional(audit, supplied, "full_support_joint_gamma_supremum", support.gamma,
                                 label, "producer_support_roots", ROOT_TOL)
                audit.check(supplied["full_support_status"] == support_status, "producer_support_status", label)
                compare_optional(audit, supplied, "full_support_intensity_ceiling",
                                 None if support.gamma is None else support.gamma * intensity,
                                 label, "producer_intensities", ROOT_TOL * max(1., intensity))
                audit.close(float(supplied["full_support_minimum_margin_at_gamma_zero"]), support.minimum_margin_at_zero,
                            "producer_corner_margins", label + " gamma=0", PROBABILITY_TOL)
                audit.close(float(supplied["full_support_minimum_margin_at_gamma_one"]), margin_one,
                            "producer_corner_margins", label + " gamma=1", PROBABILITY_TOL)
                supplied_support = nullable_number(supplied, "full_support_joint_gamma_supremum")
                if supplied_support is not None:
                    delta = max(1e-7, supplied_support * 1e-6)
                    audit.check(minimum_corner_margin(planes(max(0., supplied_support - delta))) > 0.
                                and minimum_corner_margin(planes(supplied_support + delta)) < 0.,
                                "producer_support_root_direction", label)
            if pair["neutral"]:
                audit.check(all(inverse.gamma is None for inverse in inverses.values()) and support.gamma is None,
                            "all_flat_absent", label)
            record["full_support"] = dict(gamma=support.gamma, status=support_status,
                                           existence="absent" if support.gamma is None else "finite_strict_supremum",
                                           attained=False, binding=support.binding,
                                           intensity_ceiling=None if support.gamma is None else support.gamma * intensity,
                                           minimum_margin_at_zero=support.minimum_margin_at_zero,
                                           minimum_margin_at_one=margin_one,
                                           minimum_margin_below_supremum=below, minimum_margin_above_supremum=above)
            independent_records.append(record)
            charge_records[key][charge] = record

    for key, charges in charge_records.items():
        zero, extra = charges[0.], charges[.005]
        for mode in MODES:
            a, b = zero["rank_results"][mode], extra["rank_results"][mode]
            if b["gamma"] is not None:
                audit.check(a["gamma"] is not None and b["gamma"] <= a["gamma"] + ROOT_TOL,
                            "extra_charge_threshold_order", f"{key} {mode}")
            for field in ("coverage_at_gamma_zero", "coverage_at_gamma_one", "coverage_at_gamma_three"):
                audit.check(b[field] <= a[field] + PROBABILITY_TOL, "extra_charge_coverage_order", f"{key} {mode} {field}")
        a, b = zero["full_support"]["gamma"], extra["full_support"]["gamma"]
        if b is not None:
            audit.check(a is not None and b <= a + ROOT_TOL, "extra_charge_support_order", str(key))
    audit.check(len(environment_cache) == 36 and len(pair_cache) == 28 and len(independent_records) == 432,
                "matrix_completeness", "36 environments, 28 policy pairs, 432 rows")
    audit.check(used_coefficients == set(frozen_coefficients), "matrix_completeness", "all source coefficient records used")
    if summary:
        for mode in ("positive_rank", "negative_rank", "full_support"):
            audit.check(dict(classifications[mode]) == summary["classifications"][mode], "producer_summary", mode)
        audit.close(summary["maximum_independent_baseline_difference"], 0., "producer_summary", "frozen independent baseline copied", 0.)
    for path in tracked:
        audit.check(fingerprint(path) == before[str(path)], "input_unchanged_during_verification", str(path))
    shift_summary = {}
    for mode, values in shifts.items():
        differences = [r["difference"] for r in values]
        shift_summary[mode] = dict(minimum_gamma_shift=min(differences, default=0.),
                                   maximum_gamma_shift=max(differences, default=0.),
                                   maximum_absolute_gamma_shift=max(map(abs, differences), default=0.),
                                   maximum_absolute_shift_context=max(values, key=lambda r: abs(r["difference"])) if values else None)
    return dict(producer_comparison_pending=bool(missing), missing_producer_files=missing,
                contexts=216, joint_context_charge_rows=432, financial_policies=56, environmental_vectors=36,
                independent_source_support_lp_solves=36 * 14, new_policy_optimization_solves=0,
                positive_flow_records_reaggregated=flow_count,
                classifications={mode: dict(counts) for mode, counts in classifications.items()},
                support_existence=dict(Counter("absent" if r["full_support"]["gamma"] is None else "finite" for r in independent_records)),
                maximum_independent_root_probability_residuals={mode: max(values, default=0.) for mode, values in root_residuals.items()},
                dependence_shifts_relative_to_independent=shift_summary,
                angular_probe_evaluations=angular_evaluations, angular_accepted_intervals=angular_intervals,
                independent_quad_evaluations=audit.integrals, maximum_quadpack_error_estimate=audit.quad_error,
                source_sha256=before, independent_context_records=independent_records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("financial-dir", "input-dir", "stage-facts", "travel-dir"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--results-dir", "--output-dir", dest="results_dir", type=Path, required=True)
    args = parser.parse_args()
    audit = Audit()
    report = dict(
        method="Raw source LP/cash ledger; negative-rank analytic angular crossings and interval masses; positive-rank triangular CDF; exact signed affine four-corner support ceilings",
        independence="No producer imports or producer rank/crossing routines. Only independent v30/v31 verification helpers reused.",
        angular_parameterization="u=1-cos(theta), h=1-sin(theta), theta in [0,pi/2]; density 2*sin(theta)*cos(theta)",
        probability_absolute_tolerance=PROBABILITY_TOL, gamma_absolute_tolerance=ROOT_TOL,
        coefficient_absolute_tolerance=COEFFICIENT_TOL,
        scope="Frozen v29 terminal-aware policies and declared triangular-marginal rank stresses; whole-support strict guarantee. No empirical dependence, universal copula bounds, EU radius or new policy optimization.",
        strict_support_interpretation="All planes must have strictly positive margin at every energy corner. A finite numerical boundary is an unattained supremum, not a feasible equality point.",
        runtime=dict(python=__import__("sys").executable, TEMP=os.environ.get("TEMP"), TMP=os.environ.get("TMP")),
        changed_files=[str(Path(__file__)), str(args.results_dir / "independent_joint_dependence_verification.json")],
        limitations=["Perfect positive/negative ranks are specified stresses, not extremal copula bounds.",
                     "Whole-support guarantees apply only to the declared energy square and fixed source accounting.",
                     "Distance scaling freezes flows, stock, time, loading and technology; operational geographic feasibility is not established.",
                     "Minimum raw corner margin mixes normalized environmental and EUR/service cash inequalities; only its sign is a guarantee certificate.",
                     "Background LCIA, operator cash costs and empirical correlations are not newly estimated."],
        tolerance_justification="Angular crossing formulas normalize tiny planes and use factored sine interval masses; positive CDF differences are factored. Roots use independent normalized Brent plus constant-supremum handling. Support ceilings are solved algebraically, not by producer bisection.",
    )
    try:
        report["analytical_controls"] = analytic_controls(audit)
        report.update(verify(args, audit))
    except Exception as error:
        audit.failures.append(dict(group="fatal", location=f"{type(error).__name__}: {error}"))
    pending = report.get("producer_comparison_pending", False)
    report.update(all_passed=not audit.failures and not pending, independent_checks_passed=not audit.failures,
                  failures=audit.failures, checks=dict(audit.checks), maximum_absolute_errors=dict(audit.maximum_errors))
    destination = args.results_dir / "independent_joint_dependence_verification.json"
    destination.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    compact = {key: value for key, value in report.items()
               if key not in ("analytical_controls", "source_sha256", "independent_context_records")}
    print(json.dumps(compact, indent=2, allow_nan=False))
    return 1 if audit.failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
