"""Independent vertical-slice quadrature of the frozen v30 design-prior decisions.

No integration generator, polygon clipping/moments, or financial optimizer is
imported. Source constraints are rebuilt from the printed Tua inputs. Financial
components are reaggregated from the independently verified v29 cash ledger.
Only independent_integration_verification.json is written, in --results-dir.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, linprog


PROBABILITY_TOL = 1e-9
COEFFICIENT_TOL = 2e-10
BUDGET_TOL = 1e-10  # EUR per discounted service, not campaign EUR.
TARGET = 0.95
QUAD_ABS_TOL = 2e-13
QUAD_REL_TOL = 2e-13


def read_rows(path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


class Audit:
    def __init__(self):
        self.checks = defaultdict(int)
        self.maximum_errors = defaultdict(float)
        self.failures = []

    def check(self, condition, group, label):
        self.checks[group] += 1
        if not condition:
            self.failures.append({"group": group, "location": label})

    def close(self, actual, expected, group, label, atol, rtol=0.0):
        error = abs(actual - expected)
        self.maximum_errors[group] = max(self.maximum_errors[group], error)
        self.check(
            math.isfinite(actual) and math.isfinite(expected)
            and error <= atol + rtol * max(abs(actual), abs(expected)),
            group, f"{label}: actual={actual:.17g}, expected={expected:.17g}",
        )


@dataclass(frozen=True)
class Probability:
    value: float
    error_bound: float
    intervals: int


def slice_probability(planes):
    """Integrate a*u+b*h<c on [0,1]^2 using independent 1D quadrature.

    All nonvertical boundaries, including h=0 and h=1, participate in every
    pairwise intersection. Thus the upper/lower envelope and feasibility cannot
    change inside a quadrature interval. The marginal integrand is cubic there.
    Strict nonconstant boundaries have zero prior mass; constant ties do not.
    """
    left, right = 0.0, 1.0
    lower = [(0.0, 0.0)]
    upper = [(0.0, 1.0)]
    for a, b, c in planes:
        require(all(math.isfinite(x) for x in (a, b, c)), "Nonfinite half-plane")
        scale = max(abs(a), abs(b))
        if scale == 0:
            if c <= 0:
                return Probability(0.0, 0.0, 0)
            continue
        a, b, c = a / scale, b / scale, c / scale
        if b == 0:
            intercept = c / a
            if a > 0:
                right = min(right, intercept)
            else:
                left = max(left, intercept)
        else:
            line = (-a / b, c / b)
            (upper if b > 0 else lower).append(line)
    if right <= left:
        return Probability(0.0, 0.0, 0)

    boundaries = lower + upper
    cuts = {left, right}
    for i, (slope, intercept) in enumerate(boundaries):
        for other_slope, other_intercept in boundaries[i + 1:]:
            if slope != other_slope:
                crossing = (other_intercept - intercept) / (slope - other_slope)
                if left < crossing < right:
                    cuts.add(crossing)
    cuts = sorted(cuts)

    def marginal(u):
        lo = max(slope * u + intercept for slope, intercept in lower)
        hi = min(slope * u + intercept for slope, intercept in upper)
        if hi <= lo:
            return 0.0
        # Stable integral of 4(1-u)(1-h) from lo to hi, without moment subtraction.
        return 4 * (1 - u) * (hi - lo) * (1 - (lo + hi) / 2)

    pieces, errors = [], []
    for start, stop in zip(cuts, cuts[1:]):
        value, error = quad(
            marginal, start, stop, epsabs=QUAD_ABS_TOL / (len(cuts) - 1),
            epsrel=QUAD_REL_TOL, limit=100,
        )
        pieces.append(value)
        errors.append(error)
    value, error = math.fsum(pieces), math.fsum(errors)
    require(-PROBABILITY_TOL <= value <= 1 + PROBABILITY_TOL,
            f"Quadrature probability outside [0,1]: {value}")
    return Probability(min(1.0, max(0.0, value)), error, len(pieces))


def independent_budget(planes, g, e, h, environmental_probability):
    """Find a strict-event supremum with cost slopes normalized before rooting."""
    if environmental_probability < TARGET:
        return None
    scale = max(abs(e), abs(h))
    if scale == 0:
        return g  # At this supremum the constant strict event fails, by design.
    E, H = e / scale, h / scale
    offsets = [-E * u - H * v for u in (0.0, 1.0) for v in (0.0, 1.0)]
    lo, hi = min(offsets) - 1.0, max(offsets) + 1.0

    def residual(offset):
        return slice_probability([*planes, (E, H, -offset)]).value - TARGET

    if abs(environmental_probability - TARGET) <= 5e-14:
        # Handle a possible .95 plateau by locating its upper endpoint.
        for _ in range(60):
            midpoint = (lo + hi) / 2
            if residual(midpoint) >= -5e-14:
                lo = midpoint
            else:
                hi = midpoint
        offset = (lo + hi) / 2
    else:
        offset = brentq(residual, lo, hi, xtol=2e-13, rtol=2e-14)
    return g + scale * offset


def analytic_controls(audit):
    controls = {
        "square": ([], 1.0),
        "empty": ([(0.0, 0.0, -1.0)], 0.0),
        "constant_strict_tie": ([(0.0, 0.0, 0.0)], 0.0),
        "constant_pass": ([(0.0, 0.0, 1.0)], 1.0),
        "half_axis": ([(1.0, 0.0, 0.5)], 0.75),
        "diagonal": ([(1.0, 1.0, 1.0)], 5 / 6),
        "diagonal_complement": ([(-1.0, -1.0, -1.0)], 1 / 6),
        "lower_horizontal": ([(0.0, -1.0, -0.5)], 0.25),
        "half_square": ([(1.0, 0.0, 0.5), (0.0, 1.0, 0.5)], 0.75**2),
        "inconsistent_vertical": ([(1.0, 0.0, 0.5), (-1.0, 0.0, -0.5)], 0.0),
        "small_scaled_half_axis": ([(1e-20, 0.0, 0.5e-20)], 0.75),
        "crossing_upper_envelope": ([(1.0, 1.0, 1.0), (-1.0, 1.0, 0.0)], 5 / 12),
    }
    values = {}
    for name, (planes, expected) in controls.items():
        computed = slice_probability(planes)
        audit.close(computed.value, expected, "analytical_controls", name, 1e-12)
        values[name] = computed.value
    for g, scale in [(1.0, 1.0), (0.01, 1e-8)]:
        root = independent_budget([], g, scale, 0.0, 1.0)
        expected = g - scale * (1 - math.sqrt(0.05))
        audit.close(root, expected, "analytical_root_controls", str((g, scale)), 2e-14)
    audit.check(independent_budget([], 0.0, 0.0, 0.0, 1.0) == 0.0,
                "analytical_root_controls", "constant zero supremum")
    audit.check(independent_budget([], 0.0, 1.0, 0.0, 0.5) is None,
                "analytical_root_controls", "insufficient environmental probability")
    return values


def source_constraints(input_dir, stage_facts):
    grid = read_rows(input_dir / "tua_reported_grid.csv")
    process = read_rows(input_dir / "tua_process_contributions.csv")
    stages = json.loads(stage_facts.read_text(encoding="utf-8"))["facts"]
    fixed = dict(zip([
        "Climate change", "Photochemical ozone formation", "Acidification",
        "Terrestrial eutrophication", "Freshwater ecotoxicity",
        "Cumulative energy demand", "Water depletion",
    ], [.01, .01, .01, .01, 1, 1, .01]))
    comparison = dict(zip([
        "Climate change", "Ozone depletion", "Human toxicity, non-cancer effects",
        "Human toxicity, cancer effects", "Particulate matter",
        "Photochemical ozone formation", "Acidification", "Terrestrial eutrophication",
        "Freshwater eutrophication", "Marine eutrophication", "Freshwater ecotoxicity",
        "Mineral, fossil & renewable resource depletion", "Cumulative energy demand",
        "Water depletion",
    ], [.01, 1e-7, 1e-7, 1e-7, 1e-4, 1e-3, 1e-3, 1e-3, 1e-4, 1e-4, .01, 1e-5, 1., .01]))
    polys = {}
    for stage in stages:
        indicator = stage["indicator"]
        selected = [r for r in grid if r["source_scenario"] == "W1_RE1"
                    and r["indicator"] == indicator]
        require(len(selected) == 12, f"Expected 12 source totals for {indicator}")
        reference = float(selected[0]["single_use_impact_per_service_100_crates"])
        A, b = [], []

        def constrain(vector, value, halfwidth):
            A.extend([vector, [-x for x in vector]])
            b.extend([(value + halfwidth) / reference, (-value + halfwidth) / reference])

        for row in selected:
            n = int(row["reported_rotations"])
            value = float(row["reusable_cumulative_impact_100_crates"])
            precision = fixed.get(indicator, 10**(math.floor(math.log10(value)) - 2))
            constrain([1, 1, n - 1, n - 1, n - 1], value, precision / 2)
        for n in (20, 125):
            for tag, vector in [
                ("P", [1, 0, n - 1, 0, 0]), ("E", [0, 1, 0, n - 1, 0]),
                ("W", [0, 0, 0, 0, n - 1]),
            ]:
                cell = stage["cells"][f"{tag}{n}"]
                constrain(vector, cell["value"], (cell["upper"] - cell["lower"]) / 2)
        contributions = {}
        for row in process:
            if row["indicator"] != indicator or row["source_scenario"] != "W1_RE1":
                continue
            mantissa, _, exponent = row["displayed_value"].upper().partition("E")
            exponent = int(exponent) if exponent else 0
            value = Fraction(mantissa) * Fraction(10)**exponent
            halfwidth = Fraction(1, 2) * Fraction(10)**(
                exponent - len(mantissa.partition(".")[2]))
            contributions[row["process"]] = (
                float(value - halfwidth) / reference, float(value + halfwidth) / reference)
        require(set(contributions) == {"electricity", "return_transport", "water_heating"},
                f"Incomplete source processes for {indicator}")
        polys[indicator] = (np.array(A), np.array(b),
                            1 - comparison[indicator] / (2 * reference), contributions)
    require(len(polys) == 14, "Expected 14 source indicators")
    return polys


def rebuild_financial_components(financial_dir, outputs, audit):
    verified = json.loads((financial_dir / "independent_financial_verification.json").read_text())
    audit.check(verified.get("all_passed") is True, "cash_reconstruction", "v29 independent checks")
    audit.check(verified.get("financial_schedules") == len(outputs),
                "cash_reconstruction", "v29 verified schedule count")
    parts = defaultdict(lambda: defaultdict(list))
    for row in read_rows(financial_dir / "discounted_cashflow_ledger.csv"):
        key = (row["scenario"], row["mode"])
        require(key in outputs, f"Unknown cash-ledger policy {key}")
        rate = outputs[key]["annual_discount_rate"]
        pv = float(row["cost_eur"]) * math.exp(-math.log1p(rate) * int(row["week"]) / 52)
        audit.close(pv, float(row["present_value_eur"]), "cash_invoice_pv", str(key),
                    1e-8, 1e-12)
        parts[key][row["component"]].append(pv)
    washes = defaultdict(list)
    for row in read_rows(financial_dir / "calendar_positive_flows.csv"):
        if row["variable"] == "wash":
            key = (row["scenario"], row["mode"])
            require(key in outputs, f"Unknown wash-flow policy {key}")
            week = int(row["coordinates"].split(":")[0]) + 1
            rate = outputs[key]["annual_discount_rate"]
            washes[key].append(float(row["value"]) * math.exp(-math.log1p(rate) * week / 52))
    computed = {}
    for key, result in outputs.items():
        components = {name: math.fsum(values) for name, values in parts[key].items()}
        audit.check(set(components) == set(result["present_components_eur"]),
                    "cash_reconstruction", f"{key} component names")
        for name, value in components.items():
            audit.close(value, result["present_components_eur"][name], "cash_reconstruction",
                        f"{key} {name}", 1e-7, 1e-12)
        cost = math.fsum(components.values())
        rate = result["annual_discount_rate"]
        services = math.fsum(math.fsum(demand) * math.exp(-math.log1p(rate) * (t + 1) / 52)
                             for t, demand in enumerate(result["calendar_demands"]))
        Wpv = math.fsum(washes[key])
        audit.close(cost, result["included_present_cost_eur"], "cash_reconstruction",
                    f"{key} included cost", 1e-7, 1e-12)
        audit.close(services, result["discounted_service_equivalents"], "cash_reconstruction",
                    f"{key} discounted services", 1e-7, 1e-12)
        audit.close(Wpv, result["discounted_wash_equivalents"], "cash_reconstruction",
                    f"{key} discounted washing", 1e-7, 1e-12)
        computed[key] = dict(cost=cost, transport=components["empty_transport"],
                             washes=Wpv, services=services)
    return computed


def environmental_planes(result, polys, gamma, rho):
    counts, scenario = result["counts"], result["scenario"]
    wash = counts["wash_inputs"] / 100
    vector = np.array([
        (counts["virgin_assets"] + rho * counts["missing_assets"]) / 100,
        counts["terminal_recovered_assets"] / 100,
        counts["closed_loop_rejections"] / .55,
        counts["closed_loop_rejections"] / .55, wash,
    ])
    excess_transport = (gamma * (counts["empty_dirty_crate_km"]
                        + scenario["clean_leg_factor"] * counts["empty_clean_crate_km"])
                        / 139 - counts["wash_inputs"]) / 100
    records = {}
    for indicator, (A, b, single_use_low, process) in polys.items():
        support = linprog(-vector, A_ub=A, b_ub=b,
                          bounds=[(None, None)] * 5, method="highs")
        require(support.success, f"Source support failed: {indicator}: {support.message}")
        denominator = counts["delivered_services"] / 100 * single_use_low
        baseline = (-support.fun + max(excess_transport * p
                                      for p in process["return_transport"])) / denominator
        electric = max(wash * p for p in process["electricity"]) / denominator
        heat = max(wash * p for p in process["water_heating"]) / denominator
        records[indicator] = dict(a=electric, b=heat, c=1 - baseline, base_ratio=baseline)
    return records


def decision_key(row):
    return (row["case"], float(row["travel_factor"]),
            float(row["missing_fate_virgin_equivalent_rho"]),
            float(row["electricity_price_context_eur_per_kwh"]),
            float(row["declared_effective_heat_price_eur_per_kwh"]))


def verify(args, audit, controls):
    root, financial = args.results_dir, args.financial_dir
    network_path = financial / "calendar_network_results.json"
    tracked = [network_path, financial / "discounted_cashflow_ledger.csv",
               financial / "calendar_positive_flows.csv",
               financial / "independent_financial_verification.json",
               financial / "joint_environment_financial_budgets.csv",
               args.input_dir / "tua_reported_grid.csv",
               args.input_dir / "tua_process_contributions.csv", args.stage_facts,
               root / "integrated_joint_decisions.csv", root / "environmental_halfplanes.csv",
               root / "integration_summary.json"]
    fingerprints = {str(path): sha256(path) for path in tracked}
    summary = json.loads((root / "integration_summary.json").read_text())
    audit.check(summary["source_network_sha256"] == fingerprints[str(network_path)],
                "provenance", "producer frozen-policy hash")
    results = json.loads(network_path.read_text())["outputs"]
    outputs = {(r["financial_case_id"], r["mode"]): r for r in results}
    require(len(outputs) == len(results) == 56, "Expected 56 unique v29 policies")
    financial_components = rebuild_financial_components(financial, outputs, audit)
    polys = source_constraints(args.input_dir, args.stage_facts)
    old_rows = read_rows(financial / "joint_environment_financial_budgets.csv")
    decisions = read_rows(root / "integrated_joint_decisions.csv")
    require(len(decisions) == len(old_rows) == 432, "Expected all 432 probability evaluations")
    old = {decision_key(row): row for row in old_rows}
    require(len(old) == 432 and len({decision_key(row) for row in decisions}) == 432,
            "Duplicate probability evaluation keys")
    require(set(old) == {decision_key(row) for row in decisions}, "Missing/extra decision contexts")
    emitted = read_rows(root / "environmental_halfplanes.csv")
    producer_planes = {(r["case"], float(r["travel_factor"]),
                        float(r["missing_fate_virgin_equivalent_rho"]), r["indicator"]): r
                       for r in emitted}
    require(len(producer_planes) == len(emitted) == 72 * 14,
            "Expected 1008 unique environmental coefficient records")
    env_cache, used_planes = {}, set()
    max_quad_error, total_intervals = 0.0, 0
    root_rows = constant_rows = missing_rows = flat_rows = 0
    max_root_residual = 0.0
    mc_probability_differences, mc_budget_differences, changed = [], [], []

    def integrate(planes):
        nonlocal max_quad_error, total_intervals
        probability = slice_probability(planes)
        max_quad_error = max(max_quad_error, probability.error_bound)
        total_intervals += probability.intervals
        audit.check(probability.error_bound < PROBABILITY_TOL / 10,
                    "quadrature_error", "sum of per-slice QUADPACK error estimates")
        return probability.value

    for row in decisions:
        key = decision_key(row)
        case, gamma, rho, pe, ph = key
        label = str(key)
        reference = old[key]
        P, D = outputs[case, "pooled"], outputs[case, "dedicated_adaptive"]
        for field in ("parent_scenario", "annual_discount_rate", "terminal_net_eur_per_asset"):
            audit.check(row[field] == reference[field], "context_alignment", f"{label} {field}")
        env_key = (case, gamma, rho)
        if env_key not in env_cache:
            records = environmental_planes(P, polys, gamma, rho)
            for indicator, record in records.items():
                coefficient_key = (*env_key, indicator)
                require(coefficient_key in producer_planes, f"Missing coefficient {coefficient_key}")
                used_planes.add(coefficient_key)
                for name, value in record.items():
                    audit.close(value, float(producer_planes[coefficient_key][name]),
                                "environmental_coefficients", f"{coefficient_key} {name}",
                                COEFFICIENT_TOL)
            planes = [(r["a"], r["b"], r["c"]) for r in records.values()]
            env_cache[env_key] = (planes, integrate(planes))
        planes, environment = env_cache[env_key]
        pc, dc = financial_components[case, "pooled"], financial_components[case, "dedicated_adaptive"]
        audit.close(pc["services"], dc["services"], "cash_reconstruction",
                    f"{label} common service denominator", 1e-7, 1e-12)
        Npv = pc["services"]
        wash_difference = pc["washes"] - dc["washes"]
        g = (dc["cost"] - pc["cost"] + (gamma - 1) * (dc["transport"] - pc["transport"])) / Npv
        e = .0248 * pe * wash_difference / Npv
        h = (.322 / 3.6) * ph * wash_difference / Npv
        if max(abs(g - e * u - h * v) * Npv for u in (0.0, 1.0) for v in (0.0, 1.0)) < 1e-6:
            g = e = h = 0.0
        for field, value in [("cost_base_allowance_eur_per_discounted_service", g),
                             ("cost_energy_slope_e", e), ("cost_energy_slope_h", h)]:
            audit.close(value, float(row[field]), "cost_coefficients", f"{label} {field}",
                        COEFFICIENT_TOL)
        joint = integrate([*planes, (e, h, g)])
        audit.close(environment, float(row["integrated_environmental_prior_coverage"]),
                    "environmental_probabilities", label, PROBABILITY_TOL)
        audit.close(joint, float(row["integrated_joint_prior_coverage_at_zero_extra_cost"]),
                    "joint_probabilities", label, PROBABILITY_TOL)
        expected_budget = independent_budget(planes, g, e, h, environment)
        supplied_budget = (float(row["integrated_joint95_budget_eur_per_discounted_service"])
                           if row["integrated_joint95_budget_eur_per_discounted_service"] else None)
        audit.check((supplied_budget is None) == (expected_budget is None),
                    "budget_existence", label)
        if expected_budget is None:
            missing_rows += 1
            audit.check(environment < TARGET, "budget_existence", f"{label} environmental failure")
        elif supplied_budget is not None:
            audit.close(supplied_budget, expected_budget, "normalized_budget_roots", label, BUDGET_TOL)
            scale = max(abs(e), abs(h))
            if scale == 0:
                constant_rows += 1
                audit.close(integrate([*planes, (0.0, 0.0, g - supplied_budget)]), 0.0,
                            "constant_ties", f"{label} strict equality", PROBABILITY_TOL)
                audit.close(integrate([*planes, (0.0, 0.0, 1.0)]), environment,
                            "constant_ties", f"{label} immediately below supremum", PROBABILITY_TOL)
                audit.close(integrate([*planes, (0.0, 0.0, -1.0)]), 0.0,
                            "constant_ties", f"{label} above supremum", PROBABILITY_TOL)
            else:
                root_rows += 1
                normalized = (g - supplied_budget) / scale
                root_probability = integrate([*planes, (e / scale, h / scale, normalized)])
                max_root_residual = max(max_root_residual, abs(root_probability - TARGET))
                audit.close(root_probability, TARGET, "root_probabilities", label, PROBABILITY_TOL)
                # Perturb the normalized intercept, avoiding campaign-scale root tolerances.
                below = integrate([*planes, (e / scale, h / scale, normalized + 1e-6)])
                above = integrate([*planes, (e / scale, h / scale, normalized - 1e-6)])
                audit.check(below >= TARGET - PROBABILITY_TOL and above <= TARGET + PROBABILITY_TOL,
                            "root_direction", label)
        if case.startswith("flat_"):
            flat_rows += 1
            audit.check(g == e == h == joint == 0.0 and supplied_budget == 0.0,
                        "flat_controls", label)

        for producer_field, old_field in [
            ("MC_environmental_coverage", "environmental_prior_coverage"),
            ("MC_joint_coverage", "joint_env_and_pool_cost_coverage_at_zero_extra_cost"),
        ]:
            audit.close(float(row[producer_field]), float(reference[old_field]),
                        "MC_alignment", f"{label} {producer_field}", 1e-14)
        for probability, old_field, difference_field in [
            (environment, "environmental_prior_coverage", "integrated_minus_MC_environmental_probability"),
            (joint, "joint_env_and_pool_cost_coverage_at_zero_extra_cost", "integrated_minus_MC_joint_probability"),
        ]:
            difference = probability - float(reference[old_field])
            mc_probability_differences.append(abs(difference))
            audit.close(difference, float(row[difference_field]), "difference_alignment",
                        f"{label} {difference_field}", PROBABILITY_TOL)
        prior_budget = (float(reference["joint95_budget_eur_per_discounted_service"])
                        if reference["joint95_budget_eur_per_discounted_service"] else None)
        if expected_budget is not None and prior_budget is not None:
            difference = expected_budget - prior_budget
            mc_budget_differences.append(abs(difference))
            audit.close(difference, float(row["integrated_minus_MC_budget_eur_per_discounted_service"]),
                        "difference_alignment", f"{label} budget difference", BUDGET_TOL)
        else:
            audit.check(row["integrated_minus_MC_budget_eur_per_discounted_service"] == "",
                        "difference_alignment", f"{label} absent budget difference")
        classification_changed = ((expected_budget is None) != (prior_budget is None)
                                  or (expected_budget is not None and prior_budget is not None
                                      and (expected_budget > 1e-9) != (prior_budget > 1e-9)))
        audit.check((row["classification_changed"] == "True") == classification_changed,
                    "classification", label)
        if classification_changed:
            changed.append(dict(case=case, gamma=gamma, rho=rho, pe=pe, ph=ph,
                                MC_budget=prior_budget, integrated_budget=expected_budget))
    audit.check(used_planes == set(producer_planes), "matrix_completeness", "all emitted coefficients used")
    audit.check(len(env_cache) == summary["unique_environmental_vectors"] == 72,
                "matrix_completeness", "72 environmental vectors")
    audit.check(summary["decisions"] == len(decisions), "matrix_completeness", "432 decisions")
    audit.check(flat_rows == 48, "matrix_completeness", "48 flat strict-equality controls")
    for name in summary["controls"]:
        audit.close(controls[name], summary["controls"][name], "producer_controls", name, 1e-12)
    audit.close(max(mc_probability_differences), summary["maximum_absolute_probability_difference"],
                "summary_alignment", "maximum MC probability difference", PROBABILITY_TOL)
    audit.close(max(mc_budget_differences), summary["maximum_absolute_budget_difference_eur_per_discounted_service"],
                "summary_alignment", "maximum MC budget difference", BUDGET_TOL)
    audit.check(len(changed) == len(summary["classification_changes"]),
                "summary_alignment", "classification change count")
    for path in tracked:
        audit.check(sha256(path) == fingerprints[str(path)], "provenance", f"unchanged during verification: {path}")
    return dict(
        probability_evaluations=len(decisions), financial_policies=len(outputs),
        environmental_vectors=len(env_cache), environmental_halfplanes=len(used_planes),
        nonconstant_95_percent_budget_roots=root_rows, constant_budget_suprema=constant_rows,
        absent_budgets=missing_rows, flat_control_evaluations=flat_rows,
        maximum_root_probability_residual=max_root_residual,
        maximum_quadpack_absolute_error_estimate=max_quad_error,
        reported_quadrature_intervals=total_intervals,
        classification_changes=changed, source_sha256=fingerprints,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--financial-dir", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--stage-facts", type=Path, required=True)
    args = parser.parse_args()
    audit = Audit()
    controls = analytic_controls(audit)
    producer_outputs = [args.results_dir / name for name in (
        "integrated_joint_decisions.csv", "environmental_halfplanes.csv", "integration_summary.json")]
    missing = [str(path) for path in producer_outputs if not path.exists()]
    if missing:
        print(json.dumps(dict(status="ready_waiting_for_producer", missing_files=missing,
                              analytical_controls=controls, failures=audit.failures), indent=2))
        return 1 if audit.failures else 0
    report = dict(
        method="Independent scipy.integrate.quad over vertical slices with all line-intersection breakpoints",
        density="4(1-u)(1-h) on [0,1]^2; u=lambda_E-1, h=lambda_H-1",
        scope="Frozen-policy declared design-prior integration; not empirical probability or policy optimization",
        independence="No producer or polygon-moment imports; raw Tua constraints and verified v29 cash ledger rebuilt",
        probability_absolute_tolerance=PROBABILITY_TOL,
        coefficient_absolute_tolerance=COEFFICIENT_TOL,
        budget_absolute_tolerance_eur_per_discounted_service=BUDGET_TOL,
        tolerance_justification="All envelope-switch and feasibility intersections split the cubic marginal integrand; quadrature uses 2e-13 absolute/relative targets. Budget roots use normalized slopes and independent probability residual checks.",
        analytical_controls=controls,
    )
    try:
        report.update(verify(args, audit, controls))
    except Exception as error:
        audit.failures.append(dict(group="fatal", location=f"{type(error).__name__}: {error}"))
    report.update(all_passed=not audit.failures, checks=dict(audit.checks),
                  maximum_absolute_errors=dict(audit.maximum_errors), failures=audit.failures)
    destination = args.results_dir / "independent_integration_verification.json"
    destination.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in report.items() if k != "source_sha256"}, indent=2, allow_nan=False))
    return 0 if report["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
