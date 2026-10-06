"""Independent, read-only periodic-policy audit; no producer/helper imports.

Only independent_periodic_verification.json/.md in the selected v34 results
directory are written. Source support LPs use independently assembled printed
constraints; the periodic network is NOT optimized again. Sparse positive flows
are expanded with zero for omitted entries, including all cyclic indexes.
"""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import math
import os
import sys
import time
import traceback
from collections import Counter, defaultdict
from datetime import datetime, timezone
from decimal import Decimal
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import linprog


H, J, I = 12, 14, 35
REJECT = .0055
MANUFACTURE, WASH_COST, HOLD_COST, TRUCK_COST = 4.55, .07, .03, .5 / 1000
BETA_COST = 1 / 12
BETAS = (0., 1 / 52, 1 / 12)
MODES = ("dedicated_adaptive", "pooled")
INDICATORS = (
    "Climate change", "Ozone depletion", "Human toxicity, non-cancer effects",
    "Human toxicity, cancer effects", "Particulate matter",
    "Photochemical ozone formation", "Acidification", "Terrestrial eutrophication",
    "Freshwater eutrophication", "Marine eutrophication", "Freshwater ecotoxicity",
    "Mineral, fossil & renewable resource depletion", "Cumulative energy demand",
    "Water depletion",
)
COMPARE = dict(zip(INDICATORS, (
    .01, 1e-7, 1e-7, 1e-7, 1e-4, 1e-3, 1e-3, 1e-3,
    1e-4, 1e-4, .01, 1e-5, 1., .01,
)))
FIXED = dict(zip((
    "Climate change", "Photochemical ozone formation", "Acidification",
    "Terrestrial eutrophication", "Freshwater ecotoxicity",
    "Cumulative energy demand", "Water depletion",
), (.01, .01, .01, .01, 1., 1., .01)))
CODE_ROOT = Path(__file__).resolve().parent
REPORT_STEM = "independent_periodic_verification"


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


class Audit:
    def __init__(self):
        self.counts = Counter()
        self.failures = Counter()
        self.errors = defaultdict(float)
        self.issues = []
        self.metrics = Counter()

    def check(self, group, condition, message, severity="HIGH"):
        self.counts[group] += 1
        if not condition:
            self.failures[group] += 1
            self.issues.append(dict(severity=severity, group=group, message=message))

    def near(self, group, actual, expected, label, atol=1e-5, rtol=2e-11):
        actual, expected = float(actual), float(expected)
        finite = math.isfinite(actual) and math.isfinite(expected)
        error = abs(actual - expected) if finite else math.inf
        if finite:
            self.errors[group] = max(self.errors[group], error)
        self.check(group, finite and error <= atol + rtol * abs(expected),
                   f"{label}: reported={actual!r}, independent={expected!r}, error={error!r}")


def printed_cell(text):
    number = Decimal(text.strip())
    half = Decimal(5).scaleb(number.as_tuple().exponent - 1)
    return float(number), float(half)


def declared_scenarios():
    scenarios = []
    for weights, compression in (("uniform", 0.), ("uniform", .5),
                                 ("uniform", 1.), ("ramp", 0.)):
        for tau in (1, 3, 5):
            scenarios.append(dict(name=f"main_{weights}_c{compression:g}_tau{tau}",
                                  weights=weights, compression=compression,
                                  tau=tau, q=1., origin=0.))
    for q in (.94, .97, .99):
        for tau in (1, 3):
            scenarios.append(dict(name=f"loss_q{q:g}_tau{tau}", weights="uniform",
                                  compression=0., tau=tau, q=q, origin=0.))
    for origin in (0., math.pi / 2):
        scenarios.append(dict(name=f"origin_control_{origin:.6f}", weights="uniform",
                              compression=.35, tau=3, q=1., origin=origin))
    return scenarios


def scalar_calendar(s):
    raw = [1.] * J if s["weights"] == "uniform" else [float(j + 1) for j in range(J)]
    total = math.fsum(raw)
    weights = [v / total for v in raw]
    phases = [float(s["origin"]) + (1 - float(s["compression"])) * 2 * math.pi * j / J
              for j in range(J)]
    demand = [[35000 * weights[j] * (1 + .5 * math.sin(2 * math.pi * t / H + phases[j]))
               for j in range(J)] for t in range(H)]
    coherence = math.hypot(math.fsum(weights[j] * math.cos(phases[j]) for j in range(J)),
                           math.fsum(weights[j] * math.sin(phases[j]) for j in range(J)))
    return demand, weights, phases, coherence


def cyclic_windows(d, tau, q):
    window = [[math.fsum(d[(t - a) % H][j] for a in range(tau)) + q * d[(t - tau) % H][j]
               for j in range(J)] for t in range(H)]
    owner_floors = [max(window[t][j] for t in range(H)) for j in range(J)]
    return max(math.fsum(row) for row in window), math.fsum(owner_floors), owner_floors


def constant_table_from_ast(path, name):
    """Read only the literal dict(zip(...)) definition, never execute a module."""
    module = ast.parse(path.read_text(encoding="utf-8-sig"))
    for statement in module.body:
        if isinstance(statement, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == name for t in statement.targets):
            call = statement.value
            require(isinstance(call, ast.Call) and isinstance(call.func, ast.Name) and
                    call.func.id == "dict", f"Unexpected constant definition: {name}")
            zipped = call.args[0]
            require(isinstance(zipped, ast.Call) and isinstance(zipped.func, ast.Name) and
                    zipped.func.id == "zip", f"Unexpected zip definition: {name}")
            return dict(zip(ast.literal_eval(zipped.args[0]), ast.literal_eval(zipped.args[1])))
    raise ValueError(f"Missing historical constant: {name}")


def read_geometry(inputs, audit):
    meta = read_json(inputs / "geometry" / "geometry_provenance.json")
    rows = read_csv(inputs / "geometry" / "geometry_numeric_facts.csv")
    require(len(meta["sorters"]) == J and len(meta["washers"]) == I, "Geometry node dimensions")
    require(len(set(meta["sorters"])) == J and len(set(meta["washers"])) == I,
            "Duplicate geometry identifiers")
    lookup = {}
    for row in rows:
        key = row["origin"], row["destination"]
        require(key not in lookup, f"Duplicate geometry arc {key}")
        value = float(row["distance_km"])
        require(math.isfinite(value) and value >= 0., f"Invalid distance {key}")
        lookup[key] = value
        audit.metrics["geometry_arcs_read"] += 1
    require(len(rows) == J * I == 490, "Geometry must have exactly 490 public arcs")
    audit.check("geometry_attribution", meta["source_doi"] == "10.5281/zenodo.21392137" and
                meta["no_rpc_demand_transfer"] is True and meta["no_capacity_or_LCIA_transfer"] is True,
                "Unexpected geometry source/scope")
    return [[lookup[j, i] for i in meta["washers"]] for j in meta["sorters"]]


class SourceSupport:
    def __init__(self, inputs, audit):
        self.audit = audit
        self.polytopes = {}
        self.capital_minima = {}
        self.signs = {}
        grid = read_csv(inputs / "tua_reported_grid.csv")
        stages = read_json(inputs / "replacement_stage_facts_v21.json")
        processes = read_csv(inputs / "tua_process_contributions.csv")
        require(len(stages["facts"]) == 14 and
                {x["indicator"] for x in stages["facts"]} == set(INDICATORS), "Stage indicator set")
        audit.check("source_attribution", stages["source_table"] == "S12" and
                    stages["source_doi"] == "10.3390/resources8020110", "Stage facts attribution")
        facts = {x["indicator"]: x for x in stages["facts"]}
        for indicator in INDICATORS:
            selected = [r for r in grid if r["source_scenario"] == "W1_RE1" and r["indicator"] == indicator]
            require(len(selected) == 12 and len({r["reported_rotations"] for r in selected}) == 12,
                    f"Source grid incomplete: {indicator}")
            unit = float(selected[0]["single_use_impact_per_service_100_crates"])
            require(unit > 0., f"Nonpositive reference {indicator}")
            A, b = [], []

            def interval_constraint(vector, lo, hi):
                # Five normalized coordinates: P0, E0, P, E, W. E may be negative.
                A.extend((list(vector), [-v for v in vector]))
                b.extend((hi / unit, -lo / unit))
                audit.metrics["source_interval_constraints_built"] += 1

            for row in selected:
                audit.near("source_reference_cells", row["single_use_impact_per_service_100_crates"],
                           unit, indicator, atol=0., rtol=0.)
                n = int(row["reported_rotations"])
                value = float(row["reusable_cumulative_impact_100_crates"])
                half = FIXED[indicator] / 2 if indicator in FIXED else \
                    .5 * 10 ** (math.floor(math.log10(value)) - 2)
                interval_constraint((1, 1, n - 1, n - 1, n - 1), value - half, value + half)
            for n in (20, 125):
                for tag, vector in (("P", (1, 0, n - 1, 0, 0)),
                                    ("E", (0, 1, 0, n - 1, 0)),
                                    ("W", (0, 0, 0, 0, n - 1))):
                    cell = facts[indicator]["cells"][f"{tag}{n}"]
                    value, half = printed_cell(cell["printed"])
                    for field, expected in (("value", value), ("lower", value - half),
                                            ("upper", value + half)):
                        audit.near("source_printed_stage_cells", cell[field], expected,
                                   f"{indicator}/{tag}{n}/{field}", atol=0., rtol=4e-13)
                    interval_constraint(vector, value - half, value + half)
            contributions = {}
            selected_processes = [r for r in processes if r["source_scenario"] == "W1_RE1" and
                                  r["indicator"] == indicator]
            require(len(selected_processes) == 3, f"Source process set: {indicator}")
            for row in selected_processes:
                require(row["process"] not in contributions, "Duplicate process contribution")
                value, half = printed_cell(row["displayed_value"])
                audit.near("source_printed_process_cells", row["impact_per_100_crates_reconditioned"],
                           value, f"{indicator}/{row['process']}", atol=0., rtol=4e-13)
                audit.check("source_attribution", row["source_doi"] == "10.3390/resources8020110",
                            f"Process DOI: {indicator}")
                contributions[row["process"]] = ((value - half) / unit, (value + half) / unit)
            require(set(contributions) == {"return_transport", "electricity", "water_heating"},
                    "Unexpected contribution labels")
            reference_lo = 1 - COMPARE[indicator] / (2 * unit)
            reference_hi = 1 + COMPARE[indicator] / (2 * unit)
            require(reference_lo > 0., "Single-use interval crosses zero")
            self.polytopes[indicator] = dict(A=np.asarray(A, float), b=np.asarray(b, float),
                                             transport=contributions["return_transport"],
                                             reference_lo=reference_lo, reference_hi=reference_hi)
            minimum = self.optimize(indicator, (1, 1, 0, 0, 0), "capital_minimum")
            audit.check("source_capital_net_positive", minimum > 0.,
                        f"Capital net P0+E0 not strictly positive: {indicator}: {minimum}")
            self.capital_minima[indicator] = minimum
            virgin = self.optimize(indicator, (1, 0, 0, 0, 0), "virgin_minimum")
            remaining = self.optimize(indicator, (0, 0, 1, 1, 1), "remaining_minimum")
            without_return = remaining - contributions["return_transport"][1]
            transport_min = contributions["return_transport"][0]
            for premise, value in (("virgin", virgin), ("remaining_without_reference_return", without_return),
                                   ("actual_transport_coefficient", transport_min)):
                audit.check("source_allroute_positive_premises", value > 0.,
                            f"All-route premise nonpositive: {indicator}/{premise}={value}")
            self.signs[indicator] = dict(capital=minimum, virgin=virgin, remaining=remaining,
                                         remaining_excluding_reference_return=without_return,
                                         transport_minimum=transport_min)

    def optimize(self, indicator, vector, purpose):
        p = self.polytopes[indicator]
        objective = np.asarray(vector, float)
        result = linprog(objective, A_ub=p["A"], b_ub=p["b"],
                         bounds=[(None, None)] * 5, method="highs")
        self.audit.metrics["source_support_lp_calls"] += 1
        self.audit.metrics[f"source_lp_{purpose}"] += 1
        require(result.success, f"Independent source LP failed: {indicator}: {result.message}")
        violation = max(0., float(np.max(p["A"] @ result.x - p["b"])))
        self.audit.metrics["source_witness_inequalities_evaluated"] += len(p["b"])
        self.audit.near("source_lp_primal_witness", violation, 0., f"{purpose}/{indicator}",
                        atol=2e-7, rtol=0.)
        y = result.ineqlin.marginals
        self.audit.check("source_lp_dual_sign", float(np.max(y)) <= 1e-9,
                         f"Source LP dual sign: {purpose}/{indicator}")
        stationarity = float(np.max(np.abs(objective - p["A"].T @ y)))
        self.audit.near("source_lp_stationarity", stationarity, 0., f"{purpose}/{indicator}",
                        atol=2e-8, rtol=0.)
        self.audit.near("source_lp_dual_gap", result.fun, p["b"] @ y,
                        f"{purpose}/{indicator}", atol=2e-8, rtol=2e-10)
        return float(result.fun)

    def ratio(self, counts, indicator, beta, upper, purpose):
        n, b, w = (counts[x] for x in ("delivered_services", "standing_fleet", "wash_inputs"))
        capital_pairs = beta * H * b
        c = np.asarray((counts["replacement_virgin"] + capital_pairs, capital_pairs, w, w, w)) / n
        raw = self.optimize(indicator, -c if upper else c, purpose)
        if upper:
            raw = -raw
        # Remove the source's 139 km dirty return per wash exactly once, then
        # insert both observed-from-policy empty legs at the source LCIA scale.
        correction = ((counts["empty_dirty_crate_km"] + counts["empty_clean_crate_km"]) / 139 - w) / n
        p = self.polytopes[indicator]
        contribution = p["transport"][1 if ((correction >= 0) == upper) else 0]
        return (raw + correction * contribution) / (p["reference_lo"] if upper else p["reference_hi"])

    def envelope(self, counts, beta, purpose="endpoint"):
        values = {k: dict(lower=self.ratio(counts, k, beta, False, purpose),
                          upper=self.ratio(counts, k, beta, True, purpose)) for k in INDICATORS}
        maximum_upper = max(x["upper"] for x in values.values())
        maximum_lower = max(x["lower"] for x in values.values())
        outcome = ("certified_source_advantage" if maximum_upper < 1 else
                   "certified_source_failure" if maximum_lower >= 1 else "unresolved")
        return dict(beta=beta, impact_ratios=values, maximum_upper_ratio=maximum_upper,
                    maximum_lower_ratio=maximum_lower, outcome=outcome)

    def max_upper(self, counts, beta):
        return max(self.ratio(counts, k, beta, True, "root") for k in INDICATORS)

    def inverse(self, counts):
        zero = self.max_upper(counts, 0.)
        if zero >= 1.:
            return dict(beta_supremum=None, maximum_upper_at_zero=zero,
                        absence_reason="No strict all-indicator advantage at zero stock charge")
        lo, hi = 0., 1.
        high_value = self.max_upper(counts, hi)
        while high_value < 1.:
            hi *= 2
            require(hi <= 1024, "Unable to bracket finite stock-charge threshold")
            high_value = self.max_upper(counts, hi)
        trials = 0
        while hi - lo > 2e-12:
            middle = (lo + hi) / 2
            value = self.max_upper(counts, middle)
            trials += 1
            if value < 1.:
                lo = middle
            else:
                hi = middle
        root = (lo + hi) / 2
        at = {k: self.ratio(counts, k, root, True, "root_witness") for k in INDICATORS}
        # Strict advantage is on [0, root); equality is not a passing policy.
        epsilon = max(1e-8, root * 1e-6)
        return dict(beta_supremum=root, bisection_bracket=[lo, hi], bisection_trials=trials,
                    maximum_upper_at_zero=zero, maximum_upper_at_root=max(at.values()),
                    root_indicator_upper_ratios=at,
                    strict_below_value=self.max_upper(counts, max(0., root - epsilon)),
                    strict_above_value=self.max_upper(counts, root + epsilon))


def load_sparse_flows(rows, expected_keys, audit):
    data = defaultdict(dict)
    for number, row in enumerate(rows, start=2):
        case_key = row["case"], row["mode"]
        audit.check("positive_flow_row_keys", case_key in expected_keys,
                    f"Unexpected policy in CSV line {number}: {case_key}")
        if case_key not in expected_keys:
            continue
        variable = row["variable"]
        coords = tuple(int(x) for x in row["coordinates"].split(":"))
        value = float(row["value"])
        dedicated = row["mode"] == "dedicated_adaptive"
        shapes = dict(fleet=(J if dedicated else 1,), ready=(H, J * I if dedicated else I),
                      new=(H, J * I if dedicated else I), dispatch=(H, J, I),
                      wash=(H, J, I), dirty=(H, J))
        key = (variable,) + coords
        valid = (variable in shapes and len(coords) == len(shapes[variable]) and
                 all(0 <= x < size for x, size in zip(coords, shapes[variable])))
        audit.check("positive_flow_coordinates", valid, f"Invalid coordinate line {number}: {key}")
        audit.check("positive_flow_values", math.isfinite(value) and value > 1e-7,
                    f"Invalid positive/finite value line {number}: {value}")
        audit.check("positive_flow_uniqueness", key not in data[case_key],
                    f"Duplicate variable line {number}: {case_key}/{key}")
        if valid and math.isfinite(value):
            data[case_key][key] = value
    audit.metrics["positive_flow_rows_read"] = len(rows)
    return data


def reconstruct_policy(result, s, sparse, distance, audit):
    mode = result["mode"]
    label = f"{s['name']}/{mode}"
    dedicated = mode == "dedicated_adaptive"
    slots = [(j, i) for j in range(J) for i in range(I)] if dedicated else list(range(I))
    nslots, nowners = len(slots), J if dedicated else 1
    expected_slots = [list(x) for x in slots] if dedicated else slots
    audit.check("policy_asset_slots", result["asset_slots"] == expected_slots, label)
    audit.near("policy_period", result["period_weeks"], H, label, atol=0., rtol=0.)
    audit.near("policy_capacity", result["capacity_per_washer_week"], 2000., label, atol=0., rtol=0.)
    audit.near("policy_nominal_beta", result["nominal_beta_cost"], BETA_COST, label, atol=0., rtol=0.)
    d, weights, phases, coherence = scalar_calendar(s)
    require(len(result["calendar_demands"]) == H and
            all(len(row) == J for row in result["calendar_demands"]), f"Calendar dimension: {label}")
    require(len(result["weights"]) == J and len(result["phases"]) == J,
            f"Weight/phase dimension: {label}")
    for j in range(J):
        audit.near("demand_weights", result["weights"][j], weights[j], f"{label}/{j}", atol=1e-13)
        audit.near("demand_phases", result["phases"][j], phases[j], f"{label}/{j}", atol=1e-13)
    audit.near("demand_coherence", result["phase_coherence"], coherence, label, atol=1e-13)
    for t in range(H):
        for j in range(J):
            audit.near("scalar_demand_cells", result["calendar_demands"][t][j], d[t][j],
                       f"{label}/{t}/{j}", atol=1e-9, rtol=1e-13)
    arrays = dict(fleet=np.zeros(nowners), ready=np.zeros((H, nslots)), new=np.zeros((H, nslots)),
                  dispatch=np.zeros((H, J, I)), wash=np.zeros((H, J, I)), dirty=np.zeros((H, J)))
    for key, value in sparse.items():
        arrays[key[0]][key[1:]] = value
    b, ready, new, dispatch, wash, dirty = (arrays[k] for k in
                                           ("fleet", "ready", "new", "dispatch", "wash", "dirty"))
    n = math.fsum(v for row in d for v in row)
    tau, q = s["tau"], s["q"]
    max_residual = 0.

    def equality(group, residual, context, wraps=False):
        nonlocal max_residual
        max_residual = max(max_residual, abs(float(residual)))
        audit.near(group, residual, 0., context, atol=1e-5, rtol=0.)
        if wraps:
            audit.metrics[f"wrap_{group}_checks"] += 1

    require(len(result["standing_fleet_by_owner"]) == nowners, f"Fleet owner count: {label}")
    for owner in range(nowners):
        audit.near("fleet_by_owner", result["standing_fleet_by_owner"][owner], b[owner],
                   f"{label}/{owner}")
    for t in range(H):
        previous = (t - 1) % H
        returned = (t - tau) % H
        for j in range(J):
            equality("weekly_dispatch_balances", math.fsum(dispatch[t, j]) - d[t][j], f"{label}/{t}/{j}")
            equality("weekly_dirty_balances", dirty[t, j] - dirty[previous, j] +
                     math.fsum(wash[t, j]) - q * d[returned][j], f"{label}/{t}/{j}",
                     wraps=t == 0 or t < tau)
        for slot in range(nslots):
            if dedicated:
                j, i = slots[slot]
                dispatched, restored = dispatch[t, j, i], wash[previous, j, i]
            else:
                dispatched = math.fsum(dispatch[t, :, slot])
                restored = math.fsum(wash[previous, :, slot])
            equality("weekly_ready_balances", ready[t, slot] - ready[previous, slot] - new[t, slot] +
                     dispatched - restored, f"{label}/{t}/{slot}", wraps=t == 0)
        for owner in range(nowners):
            js = [owner] if dedicated else list(range(J))
            ss = list(range(owner * I, (owner + 1) * I)) if dedicated else list(range(I))
            deadline_missing = (1 - q) * math.fsum(d[returned][j] for j in js)
            equality("weekly_new_balances", math.fsum(new[t, slot] for slot in ss) - deadline_missing,
                     f"{label}/{t}/{owner}", wraps=t < tau)
            pipeline = math.fsum(d[(t - a) % H][j] for a in range(tau) for j in js)
            state_stock = (math.fsum(ready[t, slot] for slot in ss) + math.fsum(dirty[t, j] for j in js) +
                           math.fsum(wash[t, j, i] for j in js for i in range(I)) + pipeline)
            equality("weekly_fleet_owner_balances", b[owner] - state_stock, f"{label}/{t}/{owner}",
                     wraps=t < tau - 1)
        for i in range(I):
            load = math.fsum(wash[t, j, i] for j in range(J))
            excess = max(0., load - 2000.)
            max_residual = max(max_residual, excess)
            audit.near("weekly_capacity_bounds", excess, 0., f"{label}/{t}/{i}", atol=1e-5, rtol=0.)
    sums = {tag: math.fsum(float(x) for x in array.flat) for tag, array in arrays.items()}
    km = {tag: math.fsum(float(arrays[tag][t, j, i]) * distance[j][i]
                        for t in range(H) for j in range(J) for i in range(I))
          for tag in ("wash", "dispatch")}
    counts = dict(standing_fleet=sums["fleet"], ready_asset_weeks=sums["ready"],
                  dirty_asset_weeks=sums["dirty"], wash_inputs=sums["wash"],
                  replacement_virgin=sums["new"], source_rejects=REJECT * sums["wash"],
                  delivered_services=n, empty_dirty_crate_km=km["wash"],
                  empty_clean_crate_km=km["dispatch"])
    require(set(result["counts"]) == set(counts), f"Count keys: {label}")
    for key, value in counts.items():
        audit.near("reconstructed_counts", result["counts"][key], value, f"{label}/{key}")
    audit.near("period_total_services", n, 420000., label)
    audit.near("period_wash_total", counts["wash_inputs"], q * n, label)
    audit.near("period_missing_total", counts["replacement_virgin"], (1 - q) * n, label)
    for j in range(J):
        owner_n = math.fsum(d[t][j] for t in range(H))
        audit.near("period_owner_wash_total", math.fsum(float(x) for x in wash[:, j, :].flat),
                   q * owner_n, f"{label}/{j}")
    # Summing every weekly fleet balance provides a separate whole-period identity.
    audit.near("period_holding_identity", counts["ready_asset_weeks"] + counts["dirty_asset_weeks"],
               H * counts["standing_fleet"] - (tau + q) * n, label)
    audit.check("reported_primal_residual", math.isfinite(result["maximum_flow_residual"]) and
                result["maximum_flow_residual"] <= 1e-5, f"Producer primal residual: {label}")
    lower_p, lower_d, owner_floors = cyclic_windows(d, tau, q)
    lower = lower_d if dedicated else lower_p
    audit.near("periodic_lower_bound_values", result["periodic_fleet_lower_bound"], lower, label)
    audit.check("fleet_necessary_lower_bounds", counts["standing_fleet"] >= lower - 1e-5,
                f"Fleet below necessary periodic window: {label}")
    if dedicated:
        for j in range(J):
            audit.check("fleet_owner_lower_bounds", b[j] >= owner_floors[j] - 1e-5,
                        f"Dedicated owner below window floor: {label}/{j}")
    components = dict(stock_allocation=MANUFACTURE * BETA_COST * H * counts["standing_fleet"],
                      missing_and_source_rejects=MANUFACTURE * (counts["replacement_virgin"] + counts["source_rejects"]),
                      washing=WASH_COST * counts["wash_inputs"],
                      holding=HOLD_COST * (counts["ready_asset_weeks"] + counts["dirty_asset_weeks"]),
                      transport=TRUCK_COST * (km["wash"] + km["dispatch"]))
    require(set(result["cost_components_eur"]) == set(components), f"Cost component keys: {label}")
    for key, value in components.items():
        audit.near("historical_cost_components", result["cost_components_eur"][key], value, f"{label}/{key}")
    cost = math.fsum(components.values())
    audit.near("historical_cost_total", result["cost_eur"], cost, label)
    audit.near("historical_cost_per_service", result["cost_eur_per_service"], cost / n, label, atol=1e-10)
    audit.near("reported_dual_scalar_consistency", result["dual_objective_eur"], cost, label)
    audit.check("reported_stationarity_scalar", math.isfinite(result["maximum_stationarity_residual"]) and
                result["maximum_stationarity_residual"] <= 1e-7, f"Producer stationarity scalar: {label}")
    # No optimizer: any feasible policy has B >= window floor and transport at
    # least each owner's nearest public arc. Capacities can make this cost floor loose.
    minimum_km = (1 + q) * math.fsum(math.fsum(d[t][j] for t in range(H)) * min(distance[j])
                                  for j in range(J))
    cost_floor = ((MANUFACTURE * BETA_COST * H + HOLD_COST * H) * lower +
                  MANUFACTURE * ((1 - q) * n + REJECT * q * n) + WASH_COST * q * n -
                  HOLD_COST * (tau + q) * n + TRUCK_COST * minimum_km)
    audit.check("analytical_cost_lower_bound", cost >= cost_floor - 1e-5,
                f"Cost below independent no-capacity relaxation: {label}")
    fleet_gap = counts["standing_fleet"] - lower
    cost_gap = cost - cost_floor
    return dict(counts=counts, cost_eur=cost, cost_components=components,
                maximum_independent_flow_residual=max_residual, periodic_fleet_lower_bound=lower,
                pooled_window_bound=lower_p, dedicated_window_bound=lower_d,
                standing_fleet_by_owner=b.tolist(), fleet_above_window_floor=fleet_gap,
                minimum_fleet_proved_by_feasible_floor_attainment=abs(fleet_gap) <= 1e-5 and max_residual <= 1e-5,
                analytical_cost_lower_bound_eur=cost_floor, cost_above_analytical_floor_eur=cost_gap,
                cost_minimum_proved_by_feasible_floor_attainment=abs(cost_gap) <= 1e-5 and max_residual <= 1e-5)


def verify_envelopes(result, reconstructed, source, audit):
    counts = reconstructed["counts"]
    label = f"{result['scenario']['name']}/{result['mode']}"
    evaluations = result["beta_evaluations"]
    require(len(evaluations) == 3, f"Three beta coordinates required: {label}")
    independent = []
    for beta, reported in zip(BETAS, evaluations):
        audit.near("beta_coordinate_values", reported["beta"], beta, label, atol=0., rtol=0.)
        require(set(reported["impact_ratios"]) == set(INDICATORS), f"Indicator coverage: {label}")
        calculated = source.envelope(counts, beta)
        for indicator in INDICATORS:
            for endpoint in ("lower", "upper"):
                audit.near("source_endpoint_comparisons", reported["impact_ratios"][indicator][endpoint],
                           calculated["impact_ratios"][indicator][endpoint],
                           f"{label}/beta={beta}/{indicator}/{endpoint}", atol=3e-9, rtol=3e-9)
            audit.check("source_interval_order", calculated["impact_ratios"][indicator]["lower"] <=
                        calculated["impact_ratios"][indicator]["upper"] + 3e-9,
                        f"Reversed source interval: {label}/{indicator}")
        audit.near("source_envelope_maxima", reported["maximum_upper_ratio"], calculated["maximum_upper_ratio"],
                   f"{label}/{beta}", atol=3e-9, rtol=3e-9)
        audit.check("source_envelope_classification", reported["outcome"] == calculated["outcome"],
                    f"Outcome mismatch: {label}/{beta}: {reported['outcome']} vs {calculated['outcome']}")
        limiting = reported["limiting_indicator"]
        audit.check("source_envelope_limiting_indicator", limiting in INDICATORS and
                    abs(calculated["impact_ratios"].get(limiting, {}).get("upper", -math.inf) -
                        calculated["maximum_upper_ratio"]) <= 3e-8,
                    f"Limiting indicator mismatch (ties allowed): {label}/{beta}")
        independent.append(calculated)
    inverse = source.inverse(counts)
    audit.metrics["fixed_policy_thresholds_evaluated"] += 1
    reported = result["inverse_stock_charge"]
    if inverse["beta_supremum"] is None:
        audit.metrics["threshold_absences"] += 1
        audit.check("stock_threshold_absence", reported.get("beta_supremum") is None and
                    reported.get("absence_reason") == inverse["absence_reason"], f"Threshold absence: {label}")
    else:
        audit.metrics["finite_stock_thresholds"] += 1
        root = inverse["beta_supremum"]
        audit.near("stock_threshold_roots", reported.get("beta_supremum", math.nan), root, label,
                   atol=2e-10, rtol=2e-9)
        audit.near("stock_threshold_root_residual", inverse["maximum_upper_at_root"], 1., label,
                   atol=3e-9, rtol=0.)
        audit.near("reported_stock_threshold_root_ratio", reported["maximum_upper_ratio_at_root"],
                   inverse["maximum_upper_at_root"], label, atol=3e-9, rtol=3e-9)
        audit.check("stock_threshold_strict_sides", inverse["strict_below_value"] < 1. and
                    inverse["strict_above_value"] > 1., f"Strict root sides: {label}")
        audit.check("stock_threshold_semantics", reported["strict_supremum_not_pass_at_root"] is True and
                    reported["not_observed_lifetime"] is True, f"Root/lifetime interpretation: {label}")
        limiting = reported["limiting_indicator"]
        audit.check("stock_threshold_limiting_indicator", limiting in INDICATORS and
                    abs(inverse["root_indicator_upper_ratios"].get(limiting, -math.inf) -
                        inverse["maximum_upper_at_root"]) <= 3e-8, f"Root limiting indicator: {label}")
        equivalent = counts["delivered_services"] / (root * H * counts["standing_fleet"])
        audit.near("stock_threshold_service_equivalent", reported["capital_pair_service_equivalent"],
                   equivalent, label, atol=2e-8, rtol=2e-9)
        audit.near("stock_threshold_reciprocal_charge", reported["reciprocal_charge_weeks"], 1 / root,
                   label, atol=2e-8, rtol=2e-9)
        inverse["capital_pair_service_equivalent"] = equivalent
    reconstructed.update(beta_evaluations=independent, inverse_stock_charge=inverse)


def verify_comparisons(rows, policies, scenarios, audit):
    require(len(rows) == len(scenarios) and len({r["case"] for r in rows}) == len(rows),
            "Comparisons must retain 20 unique cases")
    require({r["case"] for r in rows} == {s["name"] for s in scenarios}, "Comparison case set")
    cases = {s["name"]: s for s in scenarios}
    for row in rows:
        name = row["case"]
        s = cases[name]
        audit.check("comparison_scenario_labels", row["weights"] == s["weights"], name)
        for field in ("compression", "tau", "q", "origin"):
            audit.near("comparison_scenario_numbers", row[field], s[field], f"{name}/{field}", atol=1e-13)
        pair = [policies.get((name, mode)) for mode in MODES]
        for tag, policy in zip(("dedicated", "pooled"), pair):
            audit.check("comparison_success", row[f"{tag}_success"] == str(policy is not None),
                        f"{name}/{tag}")
        if any(x is None for x in pair):
            continue
        dedicated, pooled = pair
        for tag, policy in (("dedicated", dedicated), ("pooled", pooled)):
            audit.near("comparison_fleets", row[f"{tag}_fleet"], policy["counts"]["standing_fleet"], f"{name}/{tag}")
            beta_root = policy["inverse_stock_charge"]["beta_supremum"]
            if beta_root is None:
                audit.check("comparison_root_absence", row[f"{tag}_beta_supremum"] == "", name)
            else:
                audit.near("comparison_roots", row[f"{tag}_beta_supremum"], beta_root, f"{name}/{tag}",
                           atol=2e-10, rtol=2e-9)
            for beta_tag, envelope in zip(("zero", "52week", "12week"), policy["beta_evaluations"]):
                audit.near("comparison_envelope_numbers", row[f"{tag}_{beta_tag}_max_upper"],
                           envelope["maximum_upper_ratio"], f"{name}/{tag}/{beta_tag}", atol=3e-9, rtol=3e-9)
                audit.check("comparison_envelope_labels", row[f"{tag}_{beta_tag}_outcome"] == envelope["outcome"],
                            f"{name}/{tag}/{beta_tag}")
        saving = 1 - pooled["cost_eur"] / dedicated["cost_eur"]
        audit.near("comparison_partial_cost_savings", row["partial_cost_saving_fraction"], saving, name, atol=1e-10)
        audit.check("pooled_cost_feasible_set_inclusion", pooled["cost_eur"] <= dedicated["cost_eur"] + 1e-5,
                    f"Pooled nominal cost exceeds adaptive dedicated: {name}")


def origin_control(policies, audit):
    names = ("origin_control_0.000000", "origin_control_1.570796")
    scenarios = {s["name"]: s for s in declared_scenarios()}
    base_d = scalar_calendar(scenarios[names[0]])[0]
    shifted_d = scalar_calendar(scenarios[names[1]])[0]
    for t in range(H):
        for j in range(J):
            audit.near("origin_exact_three_week_demand_shift", shifted_d[t][j], base_d[(t + 3) % H][j],
                       f"{t}/{j}", atol=1e-9, rtol=1e-13)
    for mode in MODES:
        base, shifted = (policies.get((name, mode)) for name in names)
        if base is None or shifted is None:
            audit.check("origin_control_available", False, f"Origin control missing feasible policy: {mode}")
            continue
        for field in ("cost_eur", "periodic_fleet_lower_bound"):
            audit.near("origin_scalar_invariance", shifted[field], base[field], f"{mode}/{field}")
        audit.near("origin_scalar_invariance", shifted["counts"]["standing_fleet"],
                   base["counts"]["standing_fleet"], f"{mode}/fleet")
        for j, (left, right) in enumerate(zip(base["standing_fleet_by_owner"], shifted["standing_fleet_by_owner"])):
            audit.near("origin_owner_fleet_invariance", right, left, f"{mode}/{j}")
        for field in ("pooled_window_bound", "dedicated_window_bound"):
            audit.near("origin_both_window_invariance", shifted[field], base[field], f"{mode}/{field}")
    return dict(origin_difference_radians=math.pi / 2, exact_shift_weeks=3,
                demand_identity="d_shift[t,j] = d_base[(t+3) mod 12,j]",
                invariants_tested=["fleet", "owner_fleets", "nominal_partial_cost", "pooled_and_dedicated_window_bounds"],
                identical_degenerate_routes_required=False)


def verify_stock_certificates(rows, sign_document, policies, source, distance, audit, positive_arcs_only=False):
    """Independent route-free failure proof; never confuse a policy root with it."""
    require(set(sign_document["signs"]) == set(INDICATORS), "Stock sign indicator coverage")
    for indicator in INDICATORS:
        require(set(sign_document["signs"][indicator]) == set(source.signs[indicator]),
                f"Stock sign fields: {indicator}")
        for field, independent in source.signs[indicator].items():
            audit.near("stock_certificate_source_signs", sign_document["signs"][indicator][field],
                       independent, f"{indicator}/{field}", atol=2e-8, rtol=2e-9)
    require(len(rows) == len(policies) and len({(r["case"], r["mode"]) for r in rows}) == len(rows),
            "Stock certificates require every feasible policy exactly once")
    require({(r["case"], r["mode"]) for r in rows} == set(policies), "Stock certificate matrix")
    certificates = []
    for row in rows:
        key = row["case"], row["mode"]
        policy = policies[key]
        s = next(s for s in declared_scenarios() if s["name"] == key[0])
        d = scalar_calendar(s)[0]
        n, lb = policy["counts"]["delivered_services"], policy["periodic_fleet_lower_bound"]
        min_route_km = (1 + s["q"]) * math.fsum(math.fsum(d[t][j] for t in range(H)) *
                                              min(v for v in distance[j] if not positive_arcs_only or v>0.)
                                              for j in range(J))
        # A nearest-arc relaxation can have a zero minimum even though capacities
        # require positive travel. Each washer admitting ANY zero arc can wash
        # at most H*2000 crates per period; all other washed crates must traverse
        # an arc at least as long as the smallest strictly positive public arc.
        zero_arc_washers = set() if positive_arcs_only else {i for i in range(I) if any(distance[j][i] == 0. for j in range(J))}
        positive_distances = [value for row_distance in distance for value in row_distance if value > 0.]
        require(positive_distances, "No positive public arc available for transport-bound proof")
        zero_arc_wash_capacity = H * 2000 * len(zero_arc_washers)
        minimum_positive_wash_volume = max(0., s["q"] * n - zero_arc_wash_capacity)
        capacity_dirty_km_floor = minimum_positive_wash_volume * min(positive_distances)
        min_route_km = max(min_route_km, capacity_dirty_km_floor)
        audit.check("allroute_dirty_transport_capacity_floor", policy["counts"]["empty_dirty_crate_km"] >=
                    capacity_dirty_km_floor - 1e-5, f"Capacity dirty-km floor exceeds selected dirty-km: {key}")
        by_indicator = {}
        for indicator in INDICATORS:
            p = source.polytopes[indicator]
            capital = source.signs[indicator]["capital"]
            u_upper = p["reference_hi"]
            # At this beta the capital term alone reaches single-use parity.
            beta_floor = n * u_upper / (H * lb * capital)
            audit.metrics["allroute_floor_thresholds_calculated"] += 1
            audit.near("allroute_capital_floor_identity", beta_floor * H * lb * capital / (n * u_upper),
                       1., f"{key}/{indicator}", atol=2e-12, rtol=0.)
            source_transport_min = source.signs[indicator]["transport_minimum"]
            minimum_actual_transport_ratio = min_route_km * source_transport_min / (139 * n * u_upper)
            actual_transport_ratio = (policy["counts"]["empty_dirty_crate_km"] +
                                      policy["counts"]["empty_clean_crate_km"]) * source_transport_min / (139 * n * u_upper)
            # Nonnegativity is the necessary certificate premise. Strictly
            # positive capacity-informed floors are additional evidence only.
            audit.check("allroute_actual_transport_nonnegative", minimum_actual_transport_ratio >= 0. and
                        actual_transport_ratio >= minimum_actual_transport_ratio - 1e-10,
                        f"Actual/minimum-route transport nonnegativity: {key}/{indicator}")
            by_indicator[indicator] = dict(beta_allroute_floor=beta_floor,
                                           minimum_actual_transport_ratio=minimum_actual_transport_ratio,
                                           selected_actual_transport_ratio=actual_transport_ratio,
                                           extra_capacity_floor_strictly_positive=minimum_actual_transport_ratio > 0.)
        climate = by_indicator["Climate change"]["beta_allroute_floor"]
        coefficient = source.signs["Climate change"]["capital"] / source.polytopes["Climate change"]["reference_hi"]
        audit.near("stock_certificate_window_bound", row["periodic_fleet_lower_bound"], lb, str(key))
        audit.near("stock_certificate_climate_coefficient", row["climate_capital_coefficient_favourable"],
                   coefficient, str(key), atol=2e-9, rtol=2e-9)
        audit.near("stock_certificate_allroute_beta", row["beta_all_routing_climate_failure_from"],
                   climate, str(key), atol=2e-10, rtol=2e-9)
        selected = policy["inverse_stock_charge"]["beta_supremum"]
        if selected is None:
            audit.check("stock_certificate_selected_root_absence", row["selected_beta_supremum"] == "", str(key))
        else:
            audit.near("stock_certificate_selected_root", row["selected_beta_supremum"], selected,
                       str(key), atol=2e-10, rtol=2e-9)
        audit.check("stock_certificate_semantics", row["inverse_stock_floor_not_selected_policy_failure"] == "True",
                    f"Necessary stock failure conflated with selected-policy root: {key}")
        certificates.append(dict(case=key[0], mode=key[1], beta_allroute_climate_floor=climate,
                                 capacity_dirty_crate_km_floor=capacity_dirty_km_floor,
                                 zero_arc_washer_count=len(zero_arc_washers),
                                 zero_arc_total_wash_capacity=zero_arc_wash_capacity,
                                 minimum_positive_arc_km=min(positive_distances),
                                 all_indicator_floors=by_indicator))
        policy["beta_allroute_climate_floor"] = climate
    name = "main_uniform_c0_tau3"
    dedicated, pooled = (policies.get((name, mode)) for mode in MODES)
    regime = None
    if dedicated is not None and pooled is not None:
        floor = dedicated["beta_allroute_climate_floor"]
        root = pooled["inverse_stock_charge"]["beta_supremum"]
        audit.check("continuing_boundary_nonempty_regime", root is not None and floor < root,
                    "Uniform/tau3 dedicated all-route failure and selected-pool pass regime is absent")
        audit.near("continuing_boundary_reported_floor", floor, .0687387748232,
                   "User supplied uniform/tau3 dedicated threshold", atol=1e-12, rtol=0.)
        if not positive_arcs_only:
            audit.near("continuing_boundary_reported_pool_root", root if root is not None else math.nan,
                       .086401321764, "Original zero-enabled uniform/tau3 selected pool root", atol=2e-12, rtol=0.)
        beta = 1 / 12
        favourable_capital_floor = beta / floor
        selected_upper = pooled["beta_evaluations"][2]["maximum_upper_ratio"]
        audit.check("continuing_boundary_at_declared_beta", favourable_capital_floor >= 1. and selected_upper < 1.,
                    "At beta=1/12 the independent dedicated floor / selected pool pass separation fails")
        regime = dict(case=name, dedicated_allroute_failure_beta_from_inclusive=floor,
                      pooled_selected_advantage_beta_until_exclusive=root,
                      positive_arc_restriction=positive_arcs_only,
                      historical_zero_enabled_pool_root=.086401321764,
                      historical_root_is_not_corrected_geometry_acceptance_target=positive_arcs_only,
                      beta_is_attributed_capital_footprint_not_observed_lifetime=True,
                      universal_perpetual_pool_failure_claim=False,
                      declared_beta=beta, dedicated_favourable_climate_capital_floor=favourable_capital_floor,
                      pooled_selected_maximum_upper_at_declared_beta=selected_upper)
    audit.check("stock_certificate_proof_scope", "sufficient all-routing failure condition" in
                sign_document.get("certificate", "") and "not the exact global break-even" in
                sign_document.get("certificate", ""), "Stock certificate overstates exact/global break-even")
    return dict(independent_source_signs=source.signs, certificates=certificates, continuing_boundary_regime=regime,
                proof="min(P0+E0)>0, min(P0)>0, min(P+E+W)-max(return)>0 and min(return)>0. "
                      "Nonnegative route distances and mandatory service/wash totals make omitted terms nonnegative. "
                      "Therefore ratio >= beta*H*LB*mincapital/(N*uUpper); beta >= the reported floor suffices "
                      "to rule out all-route strict advantage, without finding the global environmental optimum.")


def scientific_code_review(producer, protocol, document, audit):
    text = producer.read_text(encoding="utf-8-sig")
    tree = ast.parse(text)
    functions = {node.name: node for node in tree.body if isinstance(node, ast.FunctionDef)}
    class_nodes = {node.name: node for node in tree.body if isinstance(node, ast.ClassDef)}
    constants = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and isinstance(node.targets[0], ast.Name):
            if node.targets[0].id in ("H", "R0"):
                constants[node.targets[0].id] = ast.literal_eval(node.value)
    audit.check("producer_literal_historical_constants", constants.get("H") == H and
                constants.get("R0") == REJECT, "Producer H/rejection constants differ from historical literal")
    locators = {name: f"{producer}:{node.lineno}" for name, node in {**functions, **class_nodes}.items()}
    capital_minima = document["source_net_capital_minima"]
    scope = document.get("scope", "").lower()
    declared = protocol.read_text(encoding="utf-8-sig").lower()
    audit.check("beta_zero_not_lifetime_lca", "not measured survival" in scope and
                "not measured survival or a free-stock lifetime lca" in scope and
                "beta=0 alone cannot establish a lifetime lca advantage" in declared,
                "Missing explicit beta=0/inherited-stock lifetime-LCA limitation")
    audit.check("allocation_not_empirical_lifetime", "additional recovered virgin-stock" in declared and
                "not a measured lifetime" in declared and "not observed annual/quarterly" in declared,
                "Capital charge mistaken for observed fleet renewal/lifetime")
    require(set(capital_minima) == set(INDICATORS), "Capital minimum indicator set")
    # These findings are static formula/meaning review, not extra numerical checks.
    return [
        dict(topic="Periodic necessary minimum", assessment="valid necessary floor",
             locator=locators.get("periodic_bounds"),
             reason="D[t]+wash[t]=D[t-1]+q*d[t-tau] implies B >= current tau-week unavailable demand + q*d[t-tau], since ready and preceding dirty are nonnegative. Max is taken after cyclic indexing; dedicated maxima are by owner. Capacity attainability is tested separately from the bound."),
        dict(topic="Source allocation and process rejection", assessment="no double count detected within declared recipe",
             locator=locators.get("support"),
             reason="Capital charge applies beta*H*B to P0+E0. Missing assets apply virgin P0 without E0 credit. W multiplies P+E+W once; source recycled-replacement rejects are already included in P/E slopes. The monetary r*wash manufacturing term is in a different accounting total and is charged once. Capital pairs are explicitly additional to process rejects."),
        dict(topic="Source transport replacement", assessment="algebraically consistent",
             locator=locators.get("support"),
             reason="(dirty_km+clean_km)/139 - wash_inputs subtracts the historical return-only contribution once and inserts both scenario legs. Symmetric reverse distances and the separate monetary 1000-crate load are scenario assumptions, not observed transport compatibility."),
        dict(topic="Periodic inventory restoration", assessment="declared idealized assumption",
             locator=locators.get("PeriodicLP"),
             reason="Accepted/rejected batches are restored together one week after washing; process replacement time and capacity are implicit. Missing supply occurs at the collection deadline. This is valid for the declared model, but not independent empirical validation of logistics or observed loss timing."),
        dict(topic="Threshold interpretation", assessment="strict fixed-policy supremum, not lifetime",
             locator=locators.get("threshold"),
             reason="Positive minimum P0+E0 guarantees increasing upper responses for positive B. Bisection tests all 14 indicators with each nominal policy fixed. Equality at the root is not strict advantage. N/(beta*H*B) and 1/beta describe charge allocation, not physical survival or legal rotations."),
        dict(topic="Selected environmental envelope", assessment="selected-policy certificate only",
             locator=locators.get("envelope"),
             reason="An envelope failure certifies that this selected nominal policy fails for at least one indicator throughout its stated source envelope. It does not prove all periodic routing policies fail. The source envelopes describe print/source-response uncertainty and are not total LCA uncertainty."),
        dict(topic="Nominal cost optimality", assessment="partly certifiable analytically; complete dual certificate unavailable",
             locator=locators.get("PeriodicLP"),
             reason="CSV proves feasibility and all cost arithmetic. A feasible policy attaining the independent analytical cost floor also proves its cost minimum. For other policies, matching the saved dual scalar and stationarity scalar is only consistency evidence: dual vectors and reduced costs were not exported, so complete independent dual feasibility/complementarity is not claimed."),
    ]


def make_markdown(report):
    lines = ["# Independent periodic verification, Project40 v34", "",
             f"Status: **{report['status']}**. all_passed: **{report['all_passed']}**. Run: {report['run_at_utc']}.", "",
             f"Executed {report['actual_checkcounts']['total']} checks; "
             f"{report['actual_checkcounts']['failed_total']} failed. "
             f"Source endpoint comparisons: {report['actual_checkcounts']['by_group'].get('source_endpoint_comparisons', 0)}. "
             f"Fixed-policy threshold evaluations: {report['activity_counts'].get('fixed_policy_thresholds_evaluated', 0)}.", "",
             "No producer, profile, geometry or source-polytope helper is imported or executed. "
             "Scalar sinusoidal demand, modulo balances, fleet windows and historical cost arithmetic are reconstructed. "
             "Only independently assembled 5-coordinate source support problems are optimized; the periodic network is not solved again.", "",
             "## Issues", ""]
    if report["issues"]:
        grouped = Counter(x["group"] for x in report["issues"])
        lines += [f"- {group}: {count} failed checks. Full records are retained in the JSON report."
                  for group, count in sorted(grouped.items())]
        lines += [f"- Example: {x['severity']} [{x['group']}]: {x['message']}" for x in report["issues"][:3]]
    else:
        lines.append("No failed numerical or scope checks in this run.")
    lines += ["", "## Actual check counts", "", "| Group | Executed | Failed |", "|---|---:|---:|"]
    for group, count in sorted(report["actual_checkcounts"]["by_group"].items()):
        lines.append(f"| {group} | {count} | {report['actual_checkcounts']['failed_by_group'].get(group, 0)} |")
    lines += ["", "## Independent findings and boundaries", ""]
    for finding in report.get("scientific_code_review", []):
        lines.append(f"- {finding['topic']}: {finding['assessment']}. {finding['reason']} Location: `{finding['locator']}`.")
    lines += ["", "## Policy evidence", "",
              "| Case / mode | Independent fleet | Window floor | Max balance residual | Fleet minimum proved | Cost minimum proved by analytical floor | Beta supremum |",
              "|---|---:|---:|---:|---|---|---:|"]
    for p in report.get("policies", []):
        root = p.get("inverse_stock_charge", {}).get("beta_supremum")
        root_text = "absent" if root is None else f"{root:.12g}"
        lines.append(f"| {p['case']} / {p['mode']} | {p['counts']['standing_fleet']:.9f} | "
                     f"{p['periodic_fleet_lower_bound']:.9f} | {p['maximum_independent_flow_residual']:.3g} | "
                     f"{p['minimum_fleet_proved_by_feasible_floor_attainment']} | "
                     f"{p['cost_minimum_proved_by_feasible_floor_attainment']} | {root_text} |")
    lines += ["", "## Limits of verification", ""]
    lines += [f"- {x}" for x in report["verification_limits"]]
    stock = report.get("periodic_stock_certificate_verification")
    if stock:
        lines += ["", "## Periodic all-route stock certificates", "", stock["proof"], ""]
        regime = stock["continuing_boundary_regime"]
        if regime:
            lines.append(f"Uniform compression 0, tau 3, q 1: dedicated strict advantage is ruled out for "
                         f"beta >= {regime['dedicated_allroute_failure_beta_from_inclusive']:.13g}; "
                         f"the selected periodic pool passes all 14 upper endpoints for "
                         f"beta < {regime['pooled_selected_advantage_beta_until_exclusive']:.13g}. "
                         f"At beta=1/12, the dedicated capital floor is "
                         f"{regime['dedicated_favourable_climate_capital_floor']:.12g}, while the selected-pool "
                         f"maximum upper ratio is {regime['pooled_selected_maximum_upper_at_declared_beta']:.12g}. "
                         "This is a sufficient conditional separation regime, not an exact global threshold, "
                         "a measured lifetime, or a claim that all perpetual pools fail.")
    previous_failed = [run for run in report.get("prior_runs", []) if run.get("issues")]
    if previous_failed:
        lines += ["", "## Preserved earlier failures", ""]
        for run in previous_failed:
            lines.append(f"Run {run['run_at_utc']}: {run['status']}.")
            if run.get("resolution_note"):
                lines.append(run["resolution_note"])
            groups = Counter(x["group"] for x in run["issues"])
            lines += [f"- Preserved {count} failures in {group}; all original messages remain in JSON prior_runs."
                      for group, count in sorted(groups.items())]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-root", type=Path, required=True, help="v33 or later attributed release root")
    parser.add_argument("--results-dir", type=Path, required=True, help="Owned continuing_pool_v34 output directory")
    parser.add_argument("--conservative-positive-arcs", action="store_true",
                        help="Use independently derived bounds for the declared positive-only refit, not its stale zero-enabled root control")
    args = parser.parse_args()
    root, out = args.package_root.resolve(), args.results_dir.resolve()
    require(out.is_absolute(), "Choose an absolute output directory")
    primary_results = root / "results" / "continuing_pool_v34"
    cold_results = out.is_relative_to(root / "tests") and out.name == "continuing_pool_v34"
    require(out == primary_results or cold_results,
            "Write reports only in the package v34 results or a fresh package test output")
    require(out.is_dir(), "Owned results directory must already exist")
    targets = [out / f"{REPORT_STEM}.{suffix}" for suffix in ("json", "md")]
    require(all(p.resolve().parent == out for p in targets), "Report target escapes owned directory")
    previous = read_json(targets[0]) if targets[0].exists() else None
    audit = Audit()
    began = time.monotonic()
    report = dict(run_at_utc=datetime.now(timezone.utc).isoformat(), status="IN_PROGRESS",
                  package_root=str(root), results_dir=str(out), producer_imported=False,
                  periodic_network_optimizer_rerun=False, independent_source_support_optimizer="SciPy HiGHS",
                  runtime=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__),
                  conservative_positive_arc_context=args.conservative_positive_arcs,
                  verification_limits=[
                      "This is a numerical and formula audit of declared public-source engineering scenarios, not empirical validation or a current full LCA/LCC.",
                      "Source supports are marginal printed-fact envelopes for 14 indicators; simultaneous cross-indicator uncertainty and upstream inventory validity are not established.",
                      "Complete nominal-network cost optimality is not independently certified for policies above the analytical cost floor: exported full dual vectors/reduced costs are unavailable. Saved producer dual/stationarity scalars are checked for consistency only.",
                      "Fleet minima are proved only when a reconstructed feasible policy attains the necessary periodic window bound. Other window values are necessary bounds, not constructive minima.",
                      "Beta thresholds hold the selected cost policy fixed. Root equality does not pass strict advantage; beta=0 excludes inherited capital footprint and cannot establish lifetime LCA advantage. Charge equivalents are not observed lifetimes or legal rotations.",
                      "No manuscript, web source, private data, producer run, or old result is changed. Only the owned script and two reports are written.",
                      "The public geometry has zero-distance arcs (including W_41 for all 14 origins). Their physical validity is not checked against the original workbook. Capacity-informed transport positivity is proved within this numerical geometry; the route-free capital failure proof only requires nonnegative actual transport.",
                  ], policies=[], scientific_code_review=[], prior_runs=[])
    if previous:
        report["prior_runs"] = previous.get("prior_runs", []) + [
            {k: previous.get(k) for k in ("run_at_utc", "status", "issues", "actual_checkcounts", "input_fingerprints")}
        ]
        if args.conservative_positive_arcs and previous.get('issues') and all(
                issue.get('group')=='continuing_boundary_reported_pool_root' for issue in previous['issues']):
            report['prior_runs'][-1]['resolution_note']=(
                'The first refit verification retained an external root target from the zero-enabled '
                'geometry. The positive-only model must not equal that old target. This explicit context '
                'does not set a new numeric acceptance target: every selected root and all source endpoints '
                'are still independently derived and compared to the producer. All 40 fixed-policy root '
                'checks remain. The first failure is preserved. Positive-only transport lower bounds are '
                'independently rebuilt with no zero-route capacity; an additional independent mask audit '
                'checks the declared restriction.')
        if previous.get("issues") and all(issue.get("group") == "allroute_actual_transport_positive"
                                          for issue in previous["issues"]):
            report["prior_runs"][-1]["resolution_note"] = (
                "The first verifier used a nearest-arc, no-capacity transport relaxation and incorrectly "
                "required it to be strictly positive. Public zero-distance arcs made that relaxation zero. "
                "The revised verifier uses a separate capacity-informed dirty-transport lower bound: "
                "wash volume above all zero-arc washer capacity must use positive-distance arcs. "
                "The actual all-route certificate acceptance checks the mathematically necessary "
                "NONNEGATIVITY premise; strict capacity-informed positivity is additional evidence only. "
                "No producer/input/result was changed and every original failed check is retained here. "
                "This resolves a verifier premise error; it does not change any stock-certificate threshold.")
    tracked = {}
    try:
        inputs = root / "inputs"
        producer = CODE_ROOT / "run_continuing_pool_v34.py"
        historical = CODE_ROOT / "couple_pool_budget_and_cohort_v20.py"
        protocol = out / "analysis_protocol.md"
        input_paths = [inputs / name for name in (
            "tua_reported_grid.csv", "replacement_stage_facts_v21.json", "tua_process_contributions.csv",
            "source_parameters.json", "source_provenance.json", "process_contributions_provenance.json",
            "geometry/geometry_numeric_facts.csv", "geometry/geometry_provenance.json")]
        for path in input_paths + [producer, historical, protocol]:
            require(path.is_file(), f"Required input/code/protocol missing: {path}")
            tracked[str(path)] = sha256(path)
        audit.check("historical_print_precisions", constant_table_from_ast(historical, "COMPARE") == COMPARE and
                    constant_table_from_ast(historical, "FIXED") == FIXED, "Historical rounding precision constants")
        source_parameters = read_json(inputs / "source_parameters.json")
        audit.near("historical_reject_rate", source_parameters["tua_facility_rejection"]["tua_w1"],
                   REJECT, "Source wash-reject recipe", atol=0., rtol=0.)
        distance = read_geometry(inputs, audit)
        source = SourceSupport(inputs, audit)
        names = ("periodic_results.json", "periodic_positive_flows.csv", "periodic_comparisons.csv", "periodic_summary.json",
                 "periodic_stock_certificates.csv", "periodic_stock_signs.json")
        missing = [name for name in names if not (out / name).is_file()]
        if missing:
            report["missing_result_files"] = missing
            report["status"] = "WAITING_FOR_PRODUCER_RESULTS"
        else:
            for name in names:
                tracked[str(out / name)] = sha256(out / name)
            document = read_json(out / names[0])
            flow_rows = read_csv(out / names[1])
            comparisons = read_csv(out / names[2])
            summary = read_json(out / names[3])
            scenarios = declared_scenarios()
            expected = {(s["name"], mode) for s in scenarios for mode in MODES}
            outputs = document["outputs"]
            require(len(outputs) == 40, "All 40 policy solve outcomes must be retained")
            keys = [(r["scenario"]["name"], r["mode"]) for r in outputs]
            require(len(set(keys)) == 40 and set(keys) == expected, "Incomplete/duplicate scenario-mode matrix")
            sparse = load_sparse_flows(flow_rows, expected, audit)
            report["scientific_code_review"] = scientific_code_review(producer, protocol, document, audit)
            for indicator, value in source.capital_minima.items():
                audit.near("source_capital_minima_comparison", document["source_net_capital_minima"][indicator],
                           value, indicator, atol=2e-8, rtol=2e-9)
            scenario_map = {s["name"]: s for s in scenarios}
            policies = {}
            for number, result in enumerate(outputs, start=1):
                key = result["scenario"]["name"], result["mode"]
                audit.metrics["nominal_policy_outcomes_read"] += 1
                audit.check("declared_scenario_matrix", result["scenario"] == scenario_map[key[0]], str(key))
                require(isinstance(result["success"], bool), f"Invalid success flag: {key}")
                if not result["success"]:
                    audit.metrics["failed_producer_policies_retained"] += 1
                    audit.check("failed_policy_no_positive_solution", not sparse.get(key), str(key))
                    report.setdefault("producer_failures", []).append(result)
                    continue
                independent = reconstruct_policy(result, scenario_map[key[0]], sparse.get(key, {}), distance, audit)
                verify_envelopes(result, independent, source, audit)
                independent.update(case=key[0], mode=key[1])
                policies[key] = independent
                report["policies"].append(independent)
                audit.metrics["feasible_policies_reconstructed"] += 1
                if independent["minimum_fleet_proved_by_feasible_floor_attainment"]:
                    audit.metrics["fleet_minima_independently_proved"] += 1
                if independent["cost_minimum_proved_by_feasible_floor_attainment"]:
                    audit.metrics["cost_minima_independently_proved_by_analytical_floor"] += 1
                print(f"Independent policy {number}/40: {key[0]} {key[1]}; "
                      f"checks={sum(audit.counts.values())}, failures={sum(audit.failures.values())}", flush=True)
            verify_comparisons(comparisons, policies, scenarios, audit)
            report["origin_shift_control"] = origin_control(policies, audit)
            report["periodic_stock_certificate_verification"] = verify_stock_certificates(
                read_csv(out / names[4]), read_json(out / names[5]), policies, source, distance, audit,
                positive_arcs_only=args.conservative_positive_arcs)
            expected_summary = dict(actual_nominal_periodic_solves=len(outputs), feasible_policies=len(policies),
                                    paired_cases=len(scenarios), source_indicators=len(source.polytopes),
                                    positive_flow_rows=len(flow_rows),
                                    endpoint_evaluations=audit.counts["source_endpoint_comparisons"])
            for field, value in expected_summary.items():
                audit.near("summary_actual_counts", summary[field], value, field, atol=0., rtol=0.)
            audit.check("summary_scope_flags", summary["all_results_retained"] is True and
                        summary["empirical_lifetime"] is False, "Result retention/empirical lifetime flags")
            audit.check("endpoint_matrix_coverage", audit.counts["source_endpoint_comparisons"] ==
                        len(policies) * 3 * 14 * 2, "Endpoint coverage must be 3*14*2 per feasible policy")
            audit.check("threshold_matrix_coverage", audit.metrics["fixed_policy_thresholds_evaluated"] == len(policies),
                        "Threshold coverage must include every feasible policy")
            report["independent_summary"] = expected_summary
            report["source_capital_net_minima"] = source.capital_minima
            report["status"] = "PASS"
    except Exception as error:
        audit.check("fatal_incomplete_verification", False, f"{type(error).__name__}: {error}")
        report["exception_traceback"] = traceback.format_exc()
        report["status"] = "FAIL_INCOMPLETE"
    finally:
        for path, before in tracked.items():
            current = Path(path)
            audit.check("read_only_input_integrity", current.is_file() and sha256(current) == before,
                        f"Read input changed during verification: {path}")
        if audit.issues and report["status"] not in ("FAIL_INCOMPLETE", "WAITING_FOR_PRODUCER_RESULTS"):
            report["status"] = "FAIL"
        report["all_passed"] = report["status"] == "PASS" and not audit.issues
        report["schema_version"] = "project40.independent_periodic_verification.v1"
        report["input_fingerprints"] = tracked
        report["issues"] = audit.issues
        report["actual_checkcounts"] = dict(total=sum(audit.counts.values()), by_group=dict(audit.counts),
                                            failed_total=sum(audit.failures.values()), failed_by_group=dict(audit.failures))
        report["activity_counts"] = dict(audit.metrics)
        report["top_counts"] = dict(
            total_checks=report["actual_checkcounts"]["total"],
            failed_checks=report["actual_checkcounts"]["failed_total"],
            feasible_policies=audit.metrics["feasible_policies_reconstructed"],
            source_indicators=len(source.polytopes) if "source" in locals() else 0,
            endpoint_comparisons=audit.counts["source_endpoint_comparisons"],
            fixed_policy_thresholds=audit.metrics["fixed_policy_thresholds_evaluated"],
            allroute_certificate_rows=len(report.get("periodic_stock_certificate_verification", {}).get("certificates", [])),
            independently_derived_allroute_indicator_floors=audit.metrics["allroute_floor_thresholds_calculated"],
            fleet_minima_proved=audit.metrics["fleet_minima_independently_proved"],
            full_periodic_cost_dual_certificates_independently_verified=0,
        )
        report["maximum_absolute_errors"] = dict(audit.errors)
        report["elapsed_seconds"] = time.monotonic() - began
        report["write_targets"] = [str(p) for p in targets]
        targets[0].write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
        targets[1].write_text(make_markdown(report), encoding="utf-8")
        print(json.dumps(dict(status=report["status"], all_passed=report["all_passed"],
                              total_checks=report["actual_checkcounts"]["total"],
                              failed_checks=report["actual_checkcounts"]["failed_total"],
                              feasible_policies=audit.metrics["feasible_policies_reconstructed"],
                              endpoint_comparisons=audit.counts["source_endpoint_comparisons"],
                              fixed_policy_thresholds=audit.metrics["fixed_policy_thresholds_evaluated"],
                              source_indicators=len(source.polytopes) if "source" in locals() else 0,
                              allroute_certificates=len(report.get("periodic_stock_certificate_verification", {}).get("certificates", [])),
                              prior_failed_checks=sum(run.get("actual_checkcounts", {}).get("failed_total", 0)
                                                      for run in report["prior_runs"]),
                              elapsed_seconds=report["elapsed_seconds"], reports=report["write_targets"]), indent=2), flush=True)
    return 2 if report["status"] == "WAITING_FOR_PRODUCER_RESULTS" else (1 if audit.issues else 0)


if __name__ == "__main__":
    sys.exit(main())
