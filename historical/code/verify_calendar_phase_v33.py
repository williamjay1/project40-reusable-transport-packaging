"""Independent arithmetic and event-ledger audit of the phase/size extension.

Only independent_phase_verification.json is written. No producer is imported.
Demand uses scalar math.sin, bounds use explicit non-circular loops, and source
support is rebuilt from attributed printed facts. The frozen independent event
verifier runs in process with its report captured and its calibration read from
the frozen release if no new calibration was supplied. This is an arithmetic
audit of declared scenarios, not empirical validation or a full LCA/LCC audit.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
import runpy
import sys
from collections import Counter, defaultdict
from contextlib import redirect_stdout
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from unittest.mock import patch

import numpy as np
from scipy.optimize import linprog


WEIGHT_NAMES = ("uniform", "ramp", "reverse_ramp")
ORIGINS = (0., math.pi / 2, 3 * math.pi / 2)
DELAYS = (1, 3, 5)
NETWORK_COMPRESSIONS = (0., .5, .9, 1.)
MODES = ("dedicated_adaptive", "pooled")
INDICATORS = (
    "Climate change", "Ozone depletion", "Human toxicity, non-cancer effects",
    "Human toxicity, cancer effects", "Particulate matter",
    "Photochemical ozone formation", "Acidification", "Terrestrial eutrophication",
    "Freshwater eutrophication", "Marine eutrophication", "Freshwater ecotoxicity",
    "Mineral, fossil & renewable resource depletion", "Cumulative energy demand",
    "Water depletion",
)
GRID_FIXED = dict(zip(
    ("Climate change", "Photochemical ozone formation", "Acidification",
     "Terrestrial eutrophication", "Freshwater ecotoxicity",
     "Cumulative energy demand", "Water depletion"),
    (.01, .01, .01, .01, 1., 1., .01),
))
SU_PRECISION = dict(zip(INDICATORS, (
    .01, 1e-7, 1e-7, 1e-7, 1e-4, 1e-3, 1e-3, 1e-3,
    1e-4, 1e-4, .01, 1e-5, 1., .01,
)))
SCREEN_FIELDS = {
    "weights", "phase_origin_radians", "compression", "return_delay_weeks",
    "phase_coherence", "delivered_services", "pooled_stock_lower_bound",
    "dedicated_stock_lower_bound", "shareable_reserve_lower_bound",
    "relative_stock_gap", "pooled_climate_floor", "dedicated_climate_floor",
    "pooled_all_routing_failure", "dedicated_all_routing_failure",
}
PAIR_FIELDS = {
    "case", "weights", "compression", "return_delay_weeks", "phase_coherence",
    "delivered_services", "dedicated_success", "pooled_success",
    "dedicated_stock_lower_bound", "pooled_stock_lower_bound",
    "dedicated_initial_stock", "pooled_initial_stock",
    "reference_partial_cost_reduction_fraction", "dedicated_maximum_upper_ratio",
    "pooled_maximum_upper_ratio", "dedicated_selected_policy_outcome",
    "pooled_selected_policy_outcome",
}


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def rows(path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def fingerprint(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


class Audit:
    def __init__(self):
        self.counts = Counter()
        self.max_absolute_error = defaultdict(float)
        self.issues = []

    def check(self, group, condition, message):
        self.counts[group] += 1
        if not condition:
            self.issues.append({"severity": "HIGH", "group": group, "message": message})

    def near(self, group, actual, expected, label, atol=1e-7, rtol=2e-11):
        actual, expected = float(actual), float(expected)
        error = abs(actual - expected)
        finite = math.isfinite(actual) and math.isfinite(expected)
        if finite:
            self.max_absolute_error[group] = max(self.max_absolute_error[group], error)
        self.check(group, finite and error <= atol + rtol * abs(expected),
                   f"{label}: observed={actual!r}, independent={expected!r}")


def parse_bool(value):
    if value is True or value == "True":
        return True
    if value is False or value == "False":
        return False
    raise ValueError(f"Invalid Boolean {value!r}")


def printed_interval(printed):
    """Base-10 fractions preserve the precision of the printed mantissa."""
    mantissa, _, exponent = printed.strip().upper().partition("E")
    power = int(exponent) if exponent else 0
    value = Fraction(mantissa) * Fraction(10) ** power
    half = Fraction(1, 2) * Fraction(10) ** (power - len(mantissa.partition(".")[2]))
    return float(value), float(half)


def independent_source(inputs, audit):
    grid = rows(inputs / "tua_reported_grid.csv")
    facts_document = read_json(inputs / "replacement_stage_facts_v21.json")
    facts = facts_document["facts"]
    processes = rows(inputs / "tua_process_contributions.csv")
    require(len(facts) == 14 and {r["indicator"] for r in facts} == set(INDICATORS),
            "Stage facts must contain each of 14 indicators exactly once")
    require(facts_document["source_table"] == "S12" and
            facts_document["source_doi"] == "10.3390/resources8020110",
            "Unexpected stage attribution")
    support, minima, lp_solves = {}, {}, 0
    for fact in facts:
        indicator = fact["indicator"]
        selected = [r for r in grid if r["source_scenario"] == "W1_RE1" and
                    r["indicator"] == indicator]
        require(len(selected) == 12 and len({r["reported_rotations"] for r in selected}) == 12,
                f"Incomplete source grid: {indicator}")
        unit = float(selected[0]["single_use_impact_per_service_100_crates"])
        require(unit > 0., f"Nonpositive single-use source denominator: {indicator}")
        matrix, rhs = [], []

        def add_interval(vector, centre, half):
            matrix.append(list(vector))
            rhs.append((centre + half) / unit)
            matrix.append([-v for v in vector])
            rhs.append((half - centre) / unit)

        for row in selected:
            audit.near("source_grid_denominators", row["single_use_impact_per_service_100_crates"],
                       unit, indicator, atol=0., rtol=0.)
            n = int(row["reported_rotations"])
            centre = float(row["reusable_cumulative_impact_100_crates"])
            half = GRID_FIXED[indicator] / 2 if indicator in GRID_FIXED else \
                10. ** (math.floor(math.log10(centre)) - 2) / 2
            add_interval((1, 1, n - 1, n - 1, n - 1), centre, half)
        for n in (20, 125):
            for tag, vector in (
                ("P", (1, 0, n - 1, 0, 0)),
                ("E", (0, 1, 0, n - 1, 0)),
                ("W", (0, 0, 0, 0, n - 1)),
            ):
                cell = fact["cells"][f"{tag}{n}"]
                centre, half = printed_interval(cell["printed"])
                for field, expected in (("value", centre), ("lower", centre - half),
                                        ("upper", centre + half)):
                    audit.near("attributed_stage_values", cell[field], expected,
                               f"{indicator}/{tag}{n}/{field}", atol=0., rtol=3e-14)
                add_interval(vector, centre, half)
        process = {}
        selected_processes = [r for r in processes if r["indicator"] == indicator and
                              r["source_scenario"] == "W1_RE1"]
        require(len(selected_processes) == 3, f"Incomplete process facts: {indicator}")
        for row in selected_processes:
            centre, half = printed_interval(row["displayed_value"])
            require(row["source_doi"] == "10.3390/resources8020110", "Unexpected process DOI")
            audit.near("attributed_process_values", row["impact_per_100_crates_reconditioned"],
                       centre, f"{indicator}/{row['process']}", atol=0., rtol=3e-14)
            require(row["process"] not in process, "Duplicate source process")
            process[row["process"]] = ((centre - half) / unit, (centre + half) / unit)
        require(set(process) == {"return_transport", "electricity", "water_heating"},
                f"Unexpected processes: {indicator}")
        mat, bounds = np.array(matrix), np.array(rhs)
        stock = linprog([1., 1., 0., 0., 0.], A_ub=mat, b_ub=bounds,
                        bounds=[(None, None)] * 5, method="highs")
        remaining = linprog([0., 0., 1., 1., 1.], A_ub=mat, b_ub=bounds,
                            bounds=[(None, None)] * 5, method="highs")
        lp_solves += 2
        require(stock.success and remaining.success, f"Source minimization failed: {indicator}")
        upper_denominator = 1. + SU_PRECISION[indicator] / (2 * unit)
        minimum_remaining = float(remaining.fun) - process["return_transport"][1]
        audit.check("source_nonnegative_remaining", minimum_remaining > 0. and
                    process["return_transport"][0] > 0., f"Floor premise fails: {indicator}")
        support[indicator] = {
            "A": mat, "b": bounds, "single_use": unit, "process": process,
            "upper_denominator": upper_denominator,
            "lower_denominator": 1. - SU_PRECISION[indicator] / (2 * unit),
        }
        minima[indicator] = {
            "initial_terminal_minimum": float(stock.fun),
            "remaining_stage_minimum": minimum_remaining,
            "transport_minimum": process["return_transport"][0],
            "upper_single_use_denominator_ratio": upper_denominator,
        }
    return support, minima, lp_solves


def weight_vector(name):
    require(name in WEIGHT_NAMES, f"Unknown weight profile: {name}")
    numbers = [1.] * 14 if name == "uniform" else [float(j + 1) for j in range(14)]
    if name == "reverse_ramp":
        numbers.reverse()
    denominator = math.fsum(numbers)
    return [value / denominator for value in numbers]


def independent_profile(name, compression, origin, tau):
    """Explicit scalar reconstruction, including zero pre-horizon occupancy."""
    weights = weight_vector(name)
    phases = [origin + (1. - compression) * 2 * math.pi * j / 14 for j in range(14)]
    demand = [[35000. * weights[j] *
               (1. + .5 * math.sin(2 * math.pi * t / 12 + phases[j]))
               for j in range(14)] for t in range(12)]
    require(all(math.isfinite(value) and value >= 0. for wave in demand for value in wave),
            "Invalid independently calculated demand")
    total = math.fsum(value for wave in demand for value in wave)
    require(abs(total - 420000.) < 1e-7, "Independent demand does not close the full period")
    owner_maximum = [0.] * 14
    aggregate_maximum = 0.
    for t in range(12):
        owner_rolling = []
        for j in range(14):
            value = 0.
            for previous in range(max(0, t - tau), t + 1):
                value += demand[previous][j]
            owner_rolling.append(value)
            owner_maximum[j] = max(owner_maximum[j], value)
        aggregate_maximum = max(aggregate_maximum, math.fsum(owner_rolling))
    dedicated = math.fsum(owner_maximum)
    coherence = math.hypot(math.fsum(weights[j] * math.cos(phases[j]) for j in range(14)),
                           math.fsum(weights[j] * math.sin(phases[j]) for j in range(14)))
    return {"weights": weights, "phases": phases, "demand": demand,
            "phase_coherence": coherence, "delivered_services": total,
            "pooled_stock_lower_bound": aggregate_maximum,
            "dedicated_stock_lower_bound": dedicated,
            "shareable_reserve_lower_bound": dedicated - aggregate_maximum,
            "relative_stock_gap": (dedicated - aggregate_maximum) / dedicated}


def origin_index(value):
    value = float(value)
    matches = [i for i, origin in enumerate(ORIGINS) if abs(origin - value) < 1e-12]
    require(len(matches) == 1, f"Origin outside declared grid: {value}")
    return matches[0]


def compression_index(value):
    value = float(value)
    step = round(value * 20)
    require(0 <= step <= 20 and abs(value - step / 20) < 1e-12,
            f"Compression outside declared grid: {value}")
    return step


def screen_key(row):
    tau = float(row["return_delay_weeks"])
    require(tau in DELAYS, f"Delay outside declared grid: {tau}")
    return (row["weights"], origin_index(row["phase_origin_radians"]), int(tau),
            compression_index(row["compression"]))


def verify_screen(root, audit, minima):
    table = rows(root / "phase_stock_screen.csv")
    require(table and set(table[0]) == SCREEN_FIELDS, "Unexpected screen schema")
    expected_keys = {(weight, origin, tau, step) for weight in WEIGHT_NAMES
                     for origin in range(3) for tau in DELAYS for step in range(21)}
    indexed = {}
    for row in table:
        key = screen_key(row)
        require(key not in indexed, f"Duplicate screen grid key: {key}")
        indexed[key] = row
    require(set(indexed) == expected_keys,
            f"Incomplete screen grid: {len(table)} rows, expected {len(expected_keys)}; "
            f"missing={sorted(expected_keys - set(indexed))[:8]}")
    # The originally requested 1701 count contradicts the stated Cartesian grid.
    # Preserve the contradiction as a failed claim check until protocol is fixed.
    protocol = (root / "analysis_protocol.md").read_text(encoding="utf-8")
    claims = [int(value) for value in re.findall(r"\b(\d+)\s+rows\b", protocol)]
    erratum = protocol.partition("## Arithmetic-count erratum")[2]
    correction = bool(erratum and re.search(r"\b567\b", erratum) and
                      "no coordinate, policy, parameter or retention rule changed" in erratum)
    audit.check("protocol_cardinality_claims", bool(claims) or correction,
                "Protocol has no explicit screen row count or documented correction")
    for claimed in claims:
        audit.check("protocol_cardinality_claims", claimed == len(expected_keys) or correction,
                    f"Uncorrected protocol count {claimed} conflicts with the complete {len(expected_keys)}-key grid")
    climate = minima["Climate change"]
    cert_counts = Counter()
    independent_values = {}
    for key, row in indexed.items():
        weight, origin, tau, step = key
        rebuilt = independent_profile(weight, step / 20, ORIGINS[origin], tau)
        factor = climate["initial_terminal_minimum"] / (
            rebuilt["delivered_services"] * climate["upper_single_use_denominator_ratio"])
        expected = {field: rebuilt[field] for field in (
            "phase_coherence", "delivered_services", "pooled_stock_lower_bound",
            "dedicated_stock_lower_bound", "shareable_reserve_lower_bound", "relative_stock_gap")}
        expected.update(pooled_climate_floor=rebuilt["pooled_stock_lower_bound"] * factor,
                        dedicated_climate_floor=rebuilt["dedicated_stock_lower_bound"] * factor)
        independent_values[key] = expected
        for field, value in expected.items():
            audit.near("screen_numeric_values", row[field], value, f"{key}/{field}",
                       atol=1e-6 if "stock" in field or "reserve" in field else 2e-9)
        for mode in ("pooled", "dedicated"):
            certificate = expected[f"{mode}_climate_floor"] >= 1.
            audit.check("screen_certificate_values", parse_bool(row[f"{mode}_all_routing_failure"]) ==
                        certificate, f"{key}/{mode}: invalid all-routing certificate")
            cert_counts[mode] += int(certificate)
        if step == 20:
            audit.near("aligned_bound_controls", rebuilt["pooled_stock_lower_bound"],
                       rebuilt["dedicated_stock_lower_bound"], str(key), atol=1e-6)
        if weight == "uniform" and step == 0:
            audit.near("uniform_stagger_bound_controls", rebuilt["pooled_stock_lower_bound"],
                       35000. * (tau + 1), str(key), atol=1e-6)
        audit.counts["screen_demand_cells_reconstructed"] += 12 * 14
        audit.counts["screen_rows_checked"] += 1
    examples, differing_classifications, compared, same_coherence_different_stock = [], [], 0, 0
    for weight in WEIGHT_NAMES:
        for tau in DELAYS:
            for step in range(21):
                for left, right in ((0, 1), (0, 2), (1, 2)):
                    a = independent_values[(weight, left, tau, step)]
                    b = independent_values[(weight, right, tau, step)]
                    label = f"{weight}/tau{tau}/c{step / 20}/origins{left}:{right}"
                    audit.near("origin_invariant_coherence", a["phase_coherence"], b["phase_coherence"],
                               label, atol=2e-14)
                    compared += 1
                    differences = {mode: b[f"{mode}_stock_lower_bound"] - a[f"{mode}_stock_lower_bound"]
                                   for mode in ("pooled", "dedicated")}
                    if max(abs(value) for value in differences.values()) <= 1e-6:
                        continue
                    same_coherence_different_stock += 1
                    example = {"weights": weight, "return_delay_weeks": tau, "compression": step / 20,
                               "phase_origin_radians": [ORIGINS[left], ORIGINS[right]],
                               "independent_phase_coherence": [a["phase_coherence"], b["phase_coherence"]],
                               "pooled_stock_lower_bounds": [a["pooled_stock_lower_bound"], b["pooled_stock_lower_bound"]],
                               "dedicated_stock_lower_bounds": [a["dedicated_stock_lower_bound"], b["dedicated_stock_lower_bound"]],
                               "pooled_climate_floors": [a["pooled_climate_floor"], b["pooled_climate_floor"]],
                               "dedicated_climate_floors": [a["dedicated_climate_floor"], b["dedicated_climate_floor"]],
                               "pooled_all_routing_climate_failure": [a["pooled_climate_floor"] >= 1., b["pooled_climate_floor"] >= 1.],
                               "dedicated_all_routing_climate_failure": [a["dedicated_climate_floor"] >= 1., b["dedicated_climate_floor"] >= 1.]}
                    examples.append(example)
                    if (a["pooled_climate_floor"] >= 1.) != (b["pooled_climate_floor"] >= 1.) or \
                            (a["dedicated_climate_floor"] >= 1.) != (b["dedicated_climate_floor"] >= 1.):
                        differing_classifications.append(example)
    return indexed, {"observed_rows": len(table), "complete_declared_grid_rows": len(expected_keys),
                     "authoritative_row_count": len(expected_keys), "documented_count_erratum": correction,
                     "all_routing_climate_failure_rows": dict(cert_counts),
                     "same_coherence_calendar_checks": {
                         "origin_pairs_compared": compared,
                         "same_coherence_different_finite_stock_pairs": same_coherence_different_stock,
                         "different_all_routing_climate_classification_pairs": len(differing_classifications),
                         "classification_can_differ": bool(differing_classifications),
                         "different_stock_examples": examples[:3],
                         "different_classification_examples": differing_classifications[:3],
                         "interpretation": "Changing origin rotates the scalar phasor without changing its magnitude; "
                                           "the finite, zero-pre-horizon campaign samples different occupancy windows. "
                                           "These are stock-floor classifications, not constructive LP classifications."}}


def verify_geometry(root, package, audit):
    meta = read_json(root / "geometry_provenance.json")
    original = read_json(package / "inputs" / "geometry" / "geometry_provenance.json")
    audit.check("geometry_provenance", meta == original, "Phase geometry provenance differs from frozen input")
    require(len(meta["sorters"]) == 14 and len(set(meta["sorters"])) == 14 and
            len(meta["washers"]) == 35 and len(set(meta["washers"])) == 35,
            "Geometry requires 14 unique customers and 35 unique washers")
    frozen = rows(package / "inputs" / "geometry" / "geometry_numeric_facts.csv")
    current = rows(root / "geometry_numeric_facts.csv")

    def index(values):
        result = {}
        for row in values:
            key = (row["origin"], row["destination"])
            require(key not in result, f"Duplicate geometry edge: {key}")
            value = float(row["distance_km"])
            require(math.isfinite(value) and value >= 0., "Invalid distance")
            result[key] = value
        return result

    target, actual = index(frozen), index(current)
    edges = {(j, i) for j in meta["sorters"] for i in meta["washers"]}
    require(set(actual) == set(target) == edges and len(actual) == 490, "Incomplete geometry")
    for key in target:
        audit.near("geometry_distance_values", actual[key], target[key], str(key), atol=0., rtol=0.)
    return [[actual[(j, i)] for i in meta["washers"]] for j in meta["sorters"]]


def check_scenario(scenario, audit):
    fixed = {"weeks": 12, "mean_weekly_demand": 35000., "seasonal_amplitude": .5,
             "seasonal_period_weeks": 12, "collection_fraction": 1., "wash_capacity_ratio": 2.,
             "total_rejection_fraction": .0055, "transport_burden_factor": 1.,
             "electricity_burden_factor": 1., "heating_burden_factor": 1., "clean_leg_factor": 1.,
             "cost_per_truck_km": .5, "declared_empty_crates_per_truck": 1000.,
             "manufacturing_eur_per_crate": 4.55, "wash_eur_per_crate": .07,
             "holding_eur_per_crate_week": .03}
    for field, expected in fixed.items():
        audit.near("network_scenario_constants", scenario[field], expected,
                   f"{scenario['name']}/{field}", atol=0., rtol=0.)
    audit.check("network_scenario_constants", scenario["seasonal_phase"] == "explicit_phase_size",
                "Wrong demand scenario marker")


def verify_network(root, audit, screen):
    document = read_json(root / "calendar_network_results.json")
    outputs = document["outputs"]
    require(len(outputs) == 32, f"Expected all 32 policies, found {len(outputs)}")
    expected = {(weight, tau, step, mode) for weight in ("uniform", "ramp")
                for tau in (1, 3) for step in (0, 10, 18, 20) for mode in MODES}
    policies, profile_by_name = {}, {}
    for result in outputs:
        s, metadata, mode = result["scenario"], result["phase_metadata"], result["mode"]
        step = compression_index(metadata["compression"])
        require(origin_index(metadata["phase_origin_radians"]) == 0, "Unexpected network origin")
        key = (metadata["weights"], s["return_delay_weeks"], step, mode)
        require(key in expected and key not in policies, f"Unexpected/duplicate policy: {key}")
        require(result["success"] is True, f"Retained unsuccessful policy cannot be certified: {key}")
        check_scenario(s, audit)
        rebuilt = independent_profile(metadata["weights"], step / 20, 0., s["return_delay_weeks"])
        require(len(metadata["node_weights"]) == len(metadata["node_phases"]) == 14,
                f"Incorrect phase metadata dimensions: {key}")
        for j in range(14):
            audit.near("explicit_metadata_values", metadata["node_weights"][j], rebuilt["weights"][j],
                       f"{key}/weight{j}", atol=2e-14)
            audit.near("explicit_metadata_values", metadata["node_phases"][j], rebuilt["phases"][j],
                       f"{key}/phase{j}", atol=2e-13)
        audit.near("explicit_metadata_values", metadata["phase_coherence"], rebuilt["phase_coherence"],
                   str(key), atol=2e-14)
        dem = result["calendar_demands"]
        require(len(dem) == 12 and all(len(wave) == 14 for wave in dem), "Wrong demand matrix dimensions")
        for t in range(12):
            for j in range(14):
                value = 35000. * float(metadata["node_weights"][j]) * (
                    1. + .5 * math.sin(2 * math.pi * t / 12 + float(metadata["node_phases"][j])))
                audit.near("network_demand_cells", dem[t][j], value, f"{key}/t{t}/j{j}", atol=1e-8)
        bound = rebuilt["dedicated_stock_lower_bound" if mode == "dedicated_adaptive" else
                        "pooled_stock_lower_bound"]
        audit.near("network_stock_bounds", result["stock_lower_bound"], bound, str(key), atol=1e-6)
        audit.check("network_stock_bounds", result["counts"]["initial_assets"] >= bound - 1e-5,
                    f"Initial stock violates owner-specific bound: {key}")
        audit.near("network_demand_totals", result["counts"]["delivered_services"],
                   rebuilt["delivered_services"], str(key), atol=1e-6)
        require(set(result["impact_ratios"]) == set(INDICATORS), "Missing impact indicator, including ozone")
        uppers = {k: float(result["impact_ratios"][k]["upper"]) for k in INDICATORS}
        lowers = {k: float(result["impact_ratios"][k]["lower"]) for k in INDICATORS}
        for indicator in INDICATORS:
            audit.check("network_impact_intervals", math.isfinite(lowers[indicator]) and
                        math.isfinite(uppers[indicator]) and lowers[indicator] <= uppers[indicator] + 1e-10,
                        f"Invalid source interval: {key}/{indicator}")
        outcome = ("certified_source_advantage" if max(uppers.values()) < 1. else
                   "certified_source_failure" if max(lowers.values()) >= 1. else "unresolved")
        audit.check("selected_policy_outcomes", result["environmental_outcome"] == outcome,
                    f"Incorrect selected-policy classification: {key}")
        audit.near("selected_policy_outcomes", result["maximum_upper_ratio"], max(uppers.values()),
                   str(key), atol=1e-10)
        audit.check("selected_policy_outcomes", uppers[result["limiting_indicator"]] == max(uppers.values()),
                    f"Incorrect limiting indicator: {key}")
        if "environmental_constraints" in result:
            audit.check("nominal_policy_scope", result["environmental_constraints"] is False,
                        "Unexpected robust optimization in nominal policies")
        audit.check("nominal_policy_scope", result["cuts"] == result["cut_iterations"] == 0,
                    "Unexpected environmental constraint cuts")
        if s["name"] in profile_by_name:
            require(profile_by_name[s["name"]] == (metadata["weights"], s["return_delay_weeks"], step),
                    "Scenario name collides across phase profiles")
        profile_by_name[s["name"]] = (metadata["weights"], s["return_delay_weeks"], step)
        policies[key] = result
    require(set(policies) == expected and len(profile_by_name) == 16, "Incomplete policy matrix")
    table = rows(root / "phase_network_comparisons.csv")
    require(len(table) == 16 and set(table[0]) == PAIR_FIELDS, "Wrong paired-comparison shape/schema")
    seen, interpretation = set(), []
    for row in table:
        key = (row["weights"], int(row["return_delay_weeks"]), compression_index(row["compression"]))
        require(key not in seen, "Duplicate paired-comparison row")
        seen.add(key)
        d, p = policies[(*key, "dedicated_adaptive")], policies[(*key, "pooled")]
        audit.check("comparison_categorical_values", row["case"] == d["scenario"]["name"] ==
                    p["scenario"]["name"], f"Comparison scenario mismatch: {key}")
        for tag, result in (("dedicated", d), ("pooled", p)):
            audit.check("comparison_categorical_values", parse_bool(row[f"{tag}_success"]) is True,
                        f"Wrong comparison solve status: {key}/{tag}")
            audit.check("comparison_categorical_values", row[f"{tag}_selected_policy_outcome"] ==
                        result["environmental_outcome"], f"Wrong selected outcome: {key}/{tag}")
        screen_row = screen[(key[0], 0, key[1], key[2])]
        expected_values = {
            "phase_coherence": d["phase_metadata"]["phase_coherence"],
            "delivered_services": d["counts"]["delivered_services"],
            "dedicated_stock_lower_bound": d["stock_lower_bound"],
            "pooled_stock_lower_bound": p["stock_lower_bound"],
            "dedicated_initial_stock": d["counts"]["initial_assets"],
            "pooled_initial_stock": p["counts"]["initial_assets"],
            "reference_partial_cost_reduction_fraction": 1. - p["cost_eur"] / d["cost_eur"],
            "dedicated_maximum_upper_ratio": d["maximum_upper_ratio"],
            "pooled_maximum_upper_ratio": p["maximum_upper_ratio"],
        }
        for field, expected_value in expected_values.items():
            audit.near("comparison_numeric_values", row[field], expected_value, f"{key}/{field}", atol=1e-7)
        audit.check("comparison_cost_dominance", p["cost_eur"] <= d["cost_eur"] + 1e-5,
                    f"Pooled nominal cost exceeds feasible dedicated benchmark: {key}")
        interpretation.append({
            "case": row["case"],
            "dedicated_all_routing_climate_failure": parse_bool(screen_row["dedicated_all_routing_failure"]),
            "pooled_all_routing_climate_failure": parse_bool(screen_row["pooled_all_routing_failure"]),
            "dedicated_selected_policy_outcome": d["environmental_outcome"],
            "pooled_selected_policy_outcome": p["environmental_outcome"],
            "certificate_logic": "All-routing flags follow mandatory-stock climate floors only; "
                                 "selected-policy outcomes concern the reported schedules.",
        })
        if key[2] == 20:
            audit.near("aligned_network_stock_controls", d["stock_lower_bound"],
                       p["stock_lower_bound"], str(key), atol=1e-6)
    return policies, interpretation


def validate_flow_rows(root, policies, audit):
    flows = rows(root / "calendar_positive_flows.csv")
    lookup = {(r["scenario"]["name"], r["mode"]): r for r in policies.values()}
    grouped, seen = defaultdict(dict), set()
    for row in flows:
        identity = (row["scenario"], row["mode"])
        require(identity in lookup, f"Orphan flow identity: {identity}")
        result = lookup[identity]
        variable = row["variable"]
        coordinates = tuple(int(v) for v in row["coordinates"].split(":"))
        slots = len(result["asset_slots"])
        if variable == "initial":
            valid = len(coordinates) == 1 and 0 <= coordinates[0] < slots
        elif variable in ("ready", "new"):
            valid = len(coordinates) == 2 and 0 <= coordinates[0] < 12 and 0 <= coordinates[1] < slots
        elif variable == "dirty":
            valid = len(coordinates) == 2 and 0 <= coordinates[0] < 12 and 0 <= coordinates[1] < 14
        elif variable in ("dispatch", "wash"):
            valid = (len(coordinates) == 3 and 0 <= coordinates[0] < 12 and
                     0 <= coordinates[1] < 14 and 0 <= coordinates[2] < 35)
        else:
            valid = False
        require(valid, f"Invalid flow coordinates: {row}")
        key = (variable, *coordinates)
        require((*identity, *key) not in seen, f"Duplicate flow: {identity}/{key}")
        seen.add((*identity, *key))
        value = float(row["value"])
        require(math.isfinite(value) and value > 0., f"Nonpositive/nonfinite positive flow: {identity}/{key}")
        grouped[identity][key] = value
        audit.counts["positive_flow_rows_validated"] += 1
    require(set(grouped) == set(lookup), "Not all 32 policies have positive flow ledgers")
    for identity, result in lookup.items():
        slots = result["asset_slots"]
        if result["mode"] == "dedicated_adaptive":
            require(len(slots) == 490 and {tuple(slot) for slot in slots} ==
                    {(j, i) for j in range(14) for i in range(35)}, "Incomplete dedicated owner slots")
        else:
            require(slots == list(range(35)), "Unexpected pooled slots")
        audit.counts["flow_slot_schemas_checked"] += 1
    return grouped


def supplemental_flow_counts(policies, flows, distance, audit):
    """Check reported counters beyond those covered by the existing verifier."""
    for result in policies.values():
        identity = (result["scenario"]["name"], result["mode"])
        ledger = flows[identity]
        totals = {tag: math.fsum(value for key, value in ledger.items() if key[0] == tag)
                  for tag in ("initial", "new", "ready", "dirty", "wash", "dispatch")}
        dirty_km = math.fsum(value * distance[key[2]][key[3]] for key, value in ledger.items()
                             if key[0] == "wash")
        clean_km = math.fsum(value * distance[key[2]][key[3]] for key, value in ledger.items()
                             if key[0] == "dispatch")
        s = result["scenario"]
        virgin = totals["initial"] + totals["new"]
        rejects = .0055 * totals["wash"]
        expected = {
            "initial_assets": totals["initial"], "replacement_virgin": totals["new"],
            "virgin_assets": virgin, "terminal_recovered_assets": virgin,
            "closed_loop_rejections": rejects, "wash_inputs": totals["wash"],
            "ready_crate_weeks": totals["ready"], "dirty_crate_weeks": totals["dirty"],
            "delivered_services": totals["dispatch"], "missing_assets": 0.,
            "empty_dirty_crate_km": dirty_km, "empty_clean_crate_km": clean_km,
            "mean_dirty_distance_km": dirty_km / totals["wash"],
            "initial_fleet_service_utilization": 420000. / totals["initial"],
            "source_mass_recycled_kg": 1.49 * (virgin + rejects),
            "source_washing_water_m3": .00055 * totals["wash"],
            "source_washing_electricity_kwh": .0248 * totals["wash"],
            "source_washing_heat_mj": .322 * totals["wash"],
        }
        require(set(result["counts"]) == set(expected), f"Unverified counter schema: {identity}")
        for field, expected_value in expected.items():
            audit.near("flow_counter_recalculations", result["counts"][field], expected_value,
                       f"{identity}/{field}", atol=1e-5, rtol=2e-12)
        components = {
            "manufacturing": 4.55 * (virgin + rejects), "reconditioning": .07 * totals["wash"],
            "holding": .03 * (totals["ready"] + totals["dirty"]),
            "empty_transport": .5 / 1000. * (dirty_km + clean_km),
        }
        require(set(result["cost_components_eur"]) == set(components), "Wrong nominal cost component schema")
        for field, value in components.items():
            audit.near("nominal_cost_component_recalculations", result["cost_components_eur"][field],
                       value, f"{identity}/{field}", atol=1e-6)
        cost = math.fsum(components.values())
        audit.near("nominal_total_cost_recalculations", result["cost_eur"], cost, str(identity), atol=1e-6)
        audit.near("nominal_total_cost_recalculations", result["cost_eur_per_service"],
                   cost / 420000., str(identity), atol=1e-10)
        audit.near("no_midcampaign_virgin_supply", totals["new"], 0., str(identity), atol=1e-8, rtol=0.)
        audit.near("reported_dual_objectives", result["dual_objective_eur"], cost,
                   str(identity), atol=1e-4, rtol=1e-8)
        audit.check("nominal_policy_scope", result["collection_fate_scope"] ==
                    "all assets collected; source terminal treatment retained", str(identity))


def existing_event_audit(root, package, audit):
    """Execute the independent frozen ledger, capturing its only write in memory.

    No result payloads, demand, geometry, flows or source facts are substituted.
    The only read redirection is the frozen calibration when the phase producer
    does not emit that unchanged file. This avoids modifying a colleague's files.
    """
    require(sys.flags.optimize == 0, "Run without -O: existing independent ledger uses active assertions")
    script = package / "code" / "verify_calendar_pool_v23.py"
    require(script.is_file(), f"Missing existing independent verifier: {script}")
    ns = runpy.run_path(str(script), run_name="independent_event_ledger")
    original_open = Path.open
    calibration = root / "source_calendar_calibration.csv"
    fallback = package / "results" / "calendar_network_v23_adaptive" / calibration.name
    actual_calibration = calibration if calibration.is_file() else fallback
    require(actual_calibration.is_file(), "No frozen calibration for the existing ledger")
    captured = {}

    def safe_write(path, content, *args, **kwargs):
        require(path.resolve() == (root / "independent_calendar_verification.json").resolve(),
                f"Existing verifier attempted an unexpected write: {path}")
        captured.update(json.loads(content))
        return len(content)

    def safe_open(path, *args, **kwargs):
        if path.resolve() == calibration.resolve():
            path = actual_calibration
        return original_open(path, *args, **kwargs)

    arguments = [str(script), "--results-dir", str(root), "--source-grid",
                 str(package / "inputs" / "tua_reported_grid.csv"), "--stage-facts",
                 str(package / "inputs" / "replacement_stage_facts_v21.json"), "--process-facts",
                 str(package / "inputs" / "tua_process_contributions.csv")]
    with patch.object(sys, "argv", arguments), patch.object(Path, "write_text", safe_write), \
            patch.object(Path, "open", safe_open), redirect_stdout(io.StringIO()):
        ns["main"]()
    require(captured.get("all_passed") is True, "Existing event audit produced no pass report")
    require(captured["successful_network_schedules"] == len(captured["checks"]) == 32,
            "Existing event audit did not check all 32 policies")
    require(all(r["impact_endpoint_checks"] == 28 for r in captured["checks"]),
            "Existing event audit did not check all 14 indicators at both endpoints")
    endpoints = sum(r["impact_endpoint_checks"] for r in captured["checks"])
    require(endpoints == 896, "Incomplete event-based source endpoints")
    audit.counts["existing_event_schedules_checked"] = len(captured["checks"])
    audit.counts["existing_event_impact_endpoint_checks"] = endpoints
    audit.counts["existing_source_calibration_checks"] = captured["source_calendar_checks"]
    return {"script": str(script), "script_sha256": fingerprint(script),
            "report_captured_in_memory_no_external_write": True,
            "calibration_source": str(actual_calibration),
            "all_32_phase_ledgers_checked": True,
            "source_endpoint_checks": endpoints, "details": captured}


def verify_frozen_controls(package, policies, current_flows, audit):
    folder = package / "results" / "calendar_network_v23_adaptive"
    frozen = read_json(folder / "calendar_network_results.json")["outputs"]
    frozen_flows = rows(folder / "calendar_positive_flows.csv")
    controls = []
    for tau in (1, 3):
        for mode in MODES:
            result = policies[("uniform", tau, 0, mode)]
            historical = next(r for r in frozen if r["scenario"]["name"] == f"staggered_tau{tau}" and
                              r["mode"] == mode)
            require(historical["success"] is True, "Frozen nominal control was not successful")
            identity = (result["scenario"]["name"], mode)
            for t in range(12):
                for j in range(14):
                    audit.near("frozen_control_demand_cells", result["calendar_demands"][t][j],
                               historical["calendar_demands"][t][j], f"{identity}/{t}/{j}", atol=1e-8)
            for field in ("cost_eur", "cost_eur_per_service", "maximum_upper_ratio"):
                audit.near("frozen_control_summary_values", result[field], historical[field],
                           f"{identity}/{field}", atol=1e-5)
            for field in historical["counts"]:
                audit.near("frozen_control_summary_values", result["counts"][field], historical["counts"][field],
                           f"{identity}/{field}", atol=1e-5)
            for indicator in INDICATORS:
                for endpoint in ("lower", "upper"):
                    audit.near("frozen_control_impact_endpoints", result["impact_ratios"][indicator][endpoint],
                               historical["impact_ratios"][indicator][endpoint],
                               f"{identity}/{indicator}/{endpoint}", atol=1e-9)
            audit.check("frozen_control_outcomes", result["environmental_outcome"] ==
                        historical["environmental_outcome"], str(identity))
            scope = [row for row in frozen_flows if row["scenario"] == historical["scenario"]["name"] and
                     row["mode"] == mode]
            require(scope, "Frozen control has no flow schedule")
            totals = {tag: math.fsum(float(row["value"]) for row in scope if row["variable"] == tag)
                      for tag in ("initial", "new", "ready", "dirty", "wash", "dispatch")}
            ledger = current_flows[identity]
            for tag in totals:
                actual = math.fsum(value for key, value in ledger.items() if key[0] == tag)
                audit.near("frozen_control_flow_aggregates", actual, totals[tag],
                           f"{identity}/{tag}", atol=1e-5)
            exact_arc_comparisons = 0
            old = {(row["variable"], *map(int, row["coordinates"].split(":"))): float(row["value"])
                   for row in scope}
            # Arc equality is recorded, not required: LP optimal schedules may
            # be degenerate, while costs/counts/endpoints are the retained controls.
            for key in set(old) | set(ledger):
                if abs(old.get(key, 0.) - ledger.get(key, 0.)) <= 1e-5:
                    exact_arc_comparisons += 1
            controls.append({"phase_scenario": identity[0], "mode": mode,
                             "frozen_nominal_scenario": historical["scenario"]["name"],
                             "frozen_positive_flow_rows": len(scope),
                             "arc_values_equal_within_1e_5": exact_arc_comparisons,
                             "arc_values_compared_for_diagnostic": len(set(old) | set(ledger)),
                             "required_match": "demand, nominal costs, all counters, all source endpoints and flow totals"})
            audit.counts["frozen_nominal_control_schedules_checked"] += 1
    return controls


def verify_summary(root, audit, minima, policies, screen_count):
    summary = read_json(root / "phase_analysis_summary.json")
    for field, expected in (("stock_screen_rows", screen_count), ("new_nominal_lp_solves", 32),
                            ("feasible_policies", 32), ("paired_cases", 16)):
        audit.near("summary_values", summary[field], expected, field, atol=0., rtol=0.)
    audit.near("summary_values", summary["minimum_climate_initial_terminal_coefficient"],
               minima["Climate change"]["initial_terminal_minimum"], "climate coefficient", atol=1e-10)
    values = summary["source_nonnegative_checks"]
    require(len(values) == 14 and {r["indicator"] for r in values} == set(INDICATORS),
            "Summary excludes indicators from source nonnegative checks")
    for row in values:
        for field in ("remaining_stage_minimum", "transport_minimum"):
            audit.near("summary_source_minima", row[field], minima[row["indicator"]][field],
                       f"{row['indicator']}/{field}", atol=1e-10)
    for field, expected in (("bounds_only_screen_no_constructive_claim", True),
                            ("all_results_retained", True), ("empirical_demand_or_current_LCC", False)):
        audit.check("summary_scope", summary[field] is expected, f"Incorrect summary scope: {field}")


def stock_attainment_and_requested_cases(policies, minima, audit):
    gaps, requested = [], []
    for key, result in policies.items():
        weight, tau, step, mode = key
        profile = independent_profile(weight, step / 20, 0., tau)
        bound = profile["dedicated_stock_lower_bound" if mode == "dedicated_adaptive" else
                        "pooled_stock_lower_bound"]
        actual = float(result["counts"]["initial_assets"])
        gaps.append({"scenario": result["scenario"]["name"], "mode": mode,
                     "independent_stock_lower_bound": bound, "actual_initial_stock": actual,
                     "actual_minus_bound": actual - bound, "attained_within_1e_5": abs(actual - bound) <= 1e-5})
    climate = minima["Climate change"]
    for weight in ("uniform", "ramp"):
        profile = independent_profile(weight, .5, 0., 3)
        factor = climate["initial_terminal_minimum"] / (
            profile["delivered_services"] * climate["upper_single_use_denominator_ratio"])
        d, p = policies[(weight, 3, 10, "dedicated_adaptive")], policies[(weight, 3, 10, "pooled")]
        dedicated_floor = profile["dedicated_stock_lower_bound"] * factor
        pooled_floor = profile["pooled_stock_lower_bound"] * factor
        requested.append({
            "scenario": p["scenario"]["name"], "weights": weight, "return_delay_weeks": 3, "compression": .5,
            "phase_coherence": profile["phase_coherence"],
            "independent_pooled_stock_lower_bound": profile["pooled_stock_lower_bound"],
            "independent_dedicated_stock_lower_bound": profile["dedicated_stock_lower_bound"],
            "independent_pooled_climate_floor": pooled_floor,
            "independent_dedicated_climate_floor": dedicated_floor,
            "pooled_all_routing_climate_failure": pooled_floor >= 1.,
            "dedicated_all_routing_climate_failure": dedicated_floor >= 1.,
            "independent_climate_initial_terminal_minimum": climate["initial_terminal_minimum"],
            "independent_climate_remaining_stage_minimum": climate["remaining_stage_minimum"],
            "independent_climate_transport_minimum": climate["transport_minimum"],
            "pooled_selected_policy_outcome": p["environmental_outcome"],
            "pooled_maximum_upper_ratio": p["maximum_upper_ratio"],
            "pooled_limiting_indicator": p["limiting_indicator"],
            "pooled_climate_ratio": p["impact_ratios"]["Climate change"],
            "pooled_ozone_depletion_ratio": p["impact_ratios"]["Ozone depletion"],
            "dedicated_selected_policy_outcome": d["environmental_outcome"],
            "all_14_pooled_endpoint_intervals_verified": True,
            "certificate_interpretation": "Climate floor >= 1 plus nonnegative source remaining terms "
                                          "certifies failure for all routing policies in that ownership model. "
                                          "Floor < 1 carries no advantage guarantee; selected outcome is separate.",
        })
        audit.counts["requested_tau3_c05_cases_checked"] += 1
    return {"policy_count": len(gaps), "attained_count_within_1e_5": sum(r["attained_within_1e_5"] for r in gaps),
            "all_32_bounds_attained_within_1e_5": all(r["attained_within_1e_5"] for r in gaps),
            "maximum_actual_minus_bound_stock_gap": max(r["actual_minus_bound"] for r in gaps),
            "minimum_actual_minus_bound_stock_gap": min(r["actual_minus_bound"] for r in gaps),
            "maximum_absolute_actual_minus_bound_stock_gap": max(abs(r["actual_minus_bound"]) for r in gaps),
            "policies": gaps}, requested


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--package-root", type=Path, required=True)
    args = parser.parse_args()
    root, package = args.results_dir.resolve(), args.package_root.resolve()
    require(root.is_absolute(), "Choose an absolute package directory")
    require(root.is_dir(), f"Results directory does not exist: {root}")
    output = root / "independent_phase_verification.json"
    audit = Audit()
    report = {"all_passed": False, "execution_completed": False,
              "verified_at_utc": datetime.now(timezone.utc).isoformat(),
              "results_dir": str(root), "package_root": str(package),
              "scope": "Declared phase/size demand, owner stock bounds, positive-flow physical/cost/source "
                       "ledgers, frozen nominal controls and climate floor logic; not observed demand, "
                       "current operator LCC, representative EU performance or global LP optimality.",
              "independence": "No producer imports. Scalar trigonometry and explicit rolling loops; "
                              "source LPs reconstructed from attributed printed facts. Existing frozen "
                              "independent verifier executed with only its calibration path and report "
                              "write redirected, without changing producer artifacts."}
    phase_names = ("analysis_protocol.md", "phase_stock_screen.csv", "phase_network_comparisons.csv",
                   "calendar_network_results.json", "calendar_positive_flows.csv",
                   "geometry_numeric_facts.csv", "geometry_provenance.json", "phase_analysis_summary.json")
    paths = [root / name for name in phase_names]
    paths += [package / "inputs" / name for name in (
        "tua_reported_grid.csv", "replacement_stage_facts_v21.json", "tua_process_contributions.csv",
        "geometry/geometry_numeric_facts.csv", "geometry/geometry_provenance.json")]
    frozen = package / "results" / "calendar_network_v23_adaptive"
    paths += [frozen / name for name in (
        "calendar_network_results.json", "calendar_positive_flows.csv", "source_calendar_calibration.csv")]
    paths.append(package / "code" / "verify_calendar_pool_v23.py")
    if (root / "source_calendar_calibration.csv").is_file():
        paths.append(root / "source_calendar_calibration.csv")
    before = {}
    try:
        missing = [str(path) for path in paths if not path.is_file()]
        require(not missing, "Artifacts not ready: " + ", ".join(missing))
        before = {str(path): fingerprint(path) for path in paths}
        support, minima, lp_solves = independent_source(package / "inputs", audit)
        report["independent_source_lp_solves"] = lp_solves
        report["independent_source_minima"] = minima
        screen, report["screen"] = verify_screen(root, audit, minima)
        distance = verify_geometry(root, package, audit)
        policies, report["certificate_vs_selected_policy"] = verify_network(root, audit, screen)
        flows = validate_flow_rows(root, policies, audit)
        supplemental_flow_counts(policies, flows, distance, audit)
        report["existing_event_ledger"] = existing_event_audit(root, package, audit)
        report["frozen_nominal_controls"] = verify_frozen_controls(package, policies, flows, audit)
        verify_summary(root, audit, minima, policies, len(screen))
        report["stock_bound_attainment"], report["requested_tau3_c05_cases"] = \
            stock_attainment_and_requested_cases(policies, minima, audit)
        for path in paths:
            audit.check("artifact_stability", fingerprint(path) == before[str(path)],
                        f"Artifact changed during verification: {path}")
        report["execution_completed"] = True
        report["all_passed"] = not audit.issues
        report["numerical_and_ledger_checks_passed"] = not any(
            issue["group"] != "protocol_cardinality_claims" for issue in audit.issues)
        report["network_policies_checked"] = len(policies)
        report["paired_cases_checked"] = len(report["certificate_vs_selected_policy"])
        report["selected_policy_outcome_counts"] = dict(Counter(
            r["environmental_outcome"] for r in policies.values()))
        report["remaining_boundaries"] = [
            "Harmonic coherence alone does not determine finite-campaign dedicated stock or route/capacity feasibility.",
            "A selected nominal-cost policy failure does not establish failure of all feasible routing policies.",
            "All-routing climate failure uses the stock floor plus source nonnegative remaining-stage and transport premises.",
            "Aligned controls require equal stock lower bounds, not equal network costs.",
            "Reported solver dual diagnostics are checked; the full nominal network LP is not independently solved.",
            "The calibration fallback verifies frozen source calibration, not a new calibration generated in this phase run.",
        ]
    except Exception as exc:
        audit.issues.append({"severity": "HIGH", "group": "execution_or_structure",
                             "message": f"{type(exc).__name__}: {exc}"})
    report["actual_checkcounts"] = dict(sorted(audit.counts.items()))
    report["maximum_absolute_errors"] = dict(sorted(audit.max_absolute_error.items()))
    report["issues"] = audit.issues
    report["issue_count"] = len(audit.issues)
    report["input_fingerprints_before_verification"] = before
    report["verifier_sha256"] = fingerprint(Path(__file__))
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"all_passed": report["all_passed"], "execution_completed": report["execution_completed"],
                      "issue_count": len(audit.issues), "issues": audit.issues[:20],
                      "actual_checkcounts": report["actual_checkcounts"], "report": str(output)}, indent=2))
    return 0 if report["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
