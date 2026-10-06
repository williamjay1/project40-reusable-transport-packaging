"""Matched joint energy-dependence and support guarantees for frozen policies.

Uses independently verified v31 affine coefficients, not new empirical inputs.
Rank stresses are NOT extremal copula bounds. Full support is a strict event
guarantee on the declared energy rectangle, not an empirical probability.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from scipy.optimize import brentq
from integrate_joint_prior_v30 import integrate
from run_joint_travel_v31 import CONTEXT_FIELDS

TARGET = .95


def rows(path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def normalize(planes):
    ans = []
    for a, b, c in planes:
        scale = max(abs(a), abs(b), abs(c))
        ans.append((a / scale, b / scale, c / scale) if scale else (0., 0., 0.))
    return ans


def triangular_cdf(x):
    x = min(1., max(0., x))
    return 2 * x - x * x


def positive_rank(planes):
    lower, upper = 0., 1.
    for a, b, c in normalize(planes):
        d = a + b
        if d > 0:
            upper = min(upper, c / d)
        elif d < 0:
            lower = max(lower, c / d)
        elif c <= 0:
            return 0.
    if upper <= lower:
        return 0.
    return max(0., triangular_cdf(upper) - triangular_cdf(lower))


def negative_rank(planes):
    """Uniform-rank interval lengths, allowing disconnected event sets."""
    planes = normalize(planes)
    cuts = {0., 1.}
    for a, b, c in planes:
        if a == b == 0:
            if c <= 0:
                return 0.
            continue

        def margin(z):
            return c - a * (1 - math.sqrt(max(0., 1 - z))) - b * (1 - math.sqrt(max(0., z)))

        sections = [0., 1.]
        if a * b > 0:
            critical = b * b / (a * a + b * b)
            if 0 < critical < 1:
                sections.insert(1, critical)
                cuts.add(critical)
        for left, right in zip(sections, sections[1:]):
            fl, fr = margin(left), margin(right)
            if fl == 0:
                cuts.add(left)
            if fr == 0:
                cuts.add(right)
            if fl * fr < 0:
                root = brentq(margin, left, right, xtol=5e-16, rtol=4e-15)
                cuts.add(root)
    cuts = sorted(cuts)
    probability = 0.
    for left, right in zip(cuts, cuts[1:]):
        midpoint = (left + right) / 2
        u = 1 - math.sqrt(max(0., 1 - midpoint))
        h = 1 - math.sqrt(max(0., midpoint))
        if all(a * u + b * h < c for a, b, c in planes):
            probability += right - left
    return min(1., max(0., probability))


def coverage_root(probability):
    if probability(0.) < TARGET:
        return None
    low, high = 0., 1.
    for _ in range(40):
        if probability(high) < TARGET:
            break
        low, high = high, 2 * high
    else:
        raise ValueError("No finite bracket: do not assert infinite headroom")
    for _ in range(70):
        middle = (low + high) / 2
        if probability(middle) >= TARGET:
            low = middle
        else:
            high = middle
    return (low + high) / 2


def support_margin(planes):
    return min(c - a * u - b * h for a, b, c in planes
               for u, h in ((0., 0.), (0., 1.), (1., 0.), (1., 1.)))


def support_root(plane_function):
    # Strict pass at every support point; a failed origin is not clipped to 0.
    if support_margin(plane_function(0.)) <= 0:
        return None
    low, high = 0., 1.
    for _ in range(40):
        if support_margin(plane_function(high)) <= 0:
            break
        low, high = high, 2 * high
    else:
        raise ValueError("No finite full-support bracket")
    for _ in range(70):
        middle = (low + high) / 2
        if support_margin(plane_function(middle)) > 0:
            low = middle
        else:
            high = middle
    return (low + high) / 2


def controls():
    checks = []
    for scale in (1., 1e-20):
        examples = [([], 1., 1.), ([(0., 0., 0.)], 0., 0.),
                    ([(0., 0., 1.)], 1., 1.), ([(1., 0., .5)], .75, .75),
                    ([(1., 1., 1.)], .75, 1.), ([(1., -1., 0.)], 0., .5),
                    ([(-1., -1., -.8)], .36, 1 - math.sqrt(.8064)),
                    ([(1., 0., .5), (-1., 0., -.5)], 0., 0.)]
        for ps, positive, negative in examples:
            ps = [(scale * a, scale * b, scale * c) for a, b, c in ps]
            pp, pn = positive_rank(ps), negative_rank(ps)
            assert abs(pp - positive) < 1e-12, (ps, pp, positive)
            assert abs(pn - negative) < 1e-12, (ps, pn, negative)
            checks.append(dict(scale=scale, positive_probability=pp, negative_probability=pn))
    for method in (positive_rank, negative_rank):
        root = coverage_root(lambda g: method([(1., 0., 1 - g)]))
        assert abs(root - math.sqrt(.05)) < 1e-13
        assert coverage_root(lambda g: method([(0., 0., -g)])) is None
        assert abs(coverage_root(lambda g: method([(0., 0., 1 - g)])) - 1.) < 1e-14
    assert abs(support_root(lambda g: [(1., 0., 2 - g)]) - 1.) < 1e-14
    assert support_root(lambda g: [(1., 0., 1 - g)]) is None
    return checks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--travel-dir", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    if (out / "joint_dependence_thresholds.csv").exists():
        raise ValueError("Fresh output required; frozen results are not overwritten")
    inputs = [args.travel_dir / n for n in ("joint_travel_thresholds.csv", "continuous_environment_coefficients.csv",
                                           "independent_joint_travel_verification.json")]
    protected = {str(p): sha(p) for p in inputs}
    assert json.loads(inputs[2].read_text())["all_passed"] is True
    old = rows(inputs[0])
    coefficients = {}
    for r in rows(inputs[1]):
        key = r["case"], float(r["missing_fate_virgin_equivalent_rho"])
        coefficients.setdefault(key, []).append(r)
    assert len(old) == 432 and len(coefficients) == 36
    analytical = controls()
    output, baseline_differences = [], []
    probes = (0., .25, 1., 3., 10., 30.)
    for r in old:
        selected = coefficients[r["case"], float(r["missing_fate_virgin_equivalent_rho"])]
        assert len(selected) == 14
        kappa = float(r["extra_net_cost_eur_per_discounted_service"])
        cash0 = float(r["cost_saving_at_gamma_zero_eur_per_discounted_service"]) - kappa
        slope = float(r["differential_transport_cost_slope_eur_per_discounted_service"])
        e, h = float(r["cost_energy_slope_e"]), float(r["cost_energy_slope_h"])
        assert slope <= 0

        def planes(gamma):
            ps = []
            for c in selected:
                t = (gamma * float(c["total_empty_crate_km"]) / 139 - float(c["wash_inputs"])) / 100
                low, high = float(c["transport_low"]), float(c["transport_high"])
                assert 0 < low <= high
                den = float(c["denominator"])
                bound = 1 - (float(c["source_support_upper"]) + max(t * low, t * high)) / den
                ps.append((float(c["electricity_slope"]), float(c["heat_slope"]), bound))
            return [*ps, (e, h, cash0 + gamma * slope)]

        result = {k: r[k] for k in [*CONTEXT_FIELDS, "extra_net_cost_eur_per_discounted_service"]}
        intensity = float(r["empty_crate_km_per_delivered_service"])
        result["empty_crate_km_per_delivered_service"] = intensity
        independent_root = None if not r["joint_gamma95"] else float(r["joint_gamma95"])
        result["independent_joint_gamma95"] = independent_root
        result["independent_joint_status"] = r["joint_status"]
        for gamma, field in ((1., "joint_coverage_at_gamma_one"), (3., "joint_coverage_at_gamma_three")):
            diff = abs(integrate(planes(gamma)) - float(r[field]))
            assert diff < 1e-12
            baseline_differences.append(diff)
        for name, method in (("positive_rank", positive_rank), ("negative_rank", negative_rank)):
            coverage = lambda gamma: method(planes(gamma))
            values = [coverage(gamma) for gamma in probes]
            assert all(a >= b - 1e-12 for a, b in zip(values, values[1:]))
            root = coverage_root(coverage)
            status = "absent_at_zero" if root is None else ("above_reference" if root > 1 else "below_reference")
            result[name + "_joint_gamma95"] = root
            result[name + "_status"] = status
            result[name + "_intensity_ceiling"] = None if root is None else root * intensity
            for gamma, word in ((0., "zero"), (1., "one"), (3., "three")):
                result[name + "_coverage_at_gamma_" + word] = coverage(gamma)
            result[name + "_coverage_at_supremum"] = None if root is None else coverage(root)
        robust = support_root(planes)
        result["full_support_joint_gamma_supremum"] = robust
        result["full_support_status"] = "absent_at_zero" if robust is None else ("above_reference" if robust > 1 else "below_reference")
        result["full_support_intensity_ceiling"] = None if robust is None else robust * intensity
        result["full_support_minimum_margin_at_gamma_zero"] = support_margin(planes(0.))
        result["full_support_minimum_margin_at_gamma_one"] = support_margin(planes(1.))
        if robust is not None:
            for root in (independent_root, result["positive_rank_joint_gamma95"], result["negative_rank_joint_gamma95"]):
                assert root is not None and robust <= root + 1e-9
            delta = max(1e-7, robust * 1e-6)
            assert support_margin(planes(max(0., robust - delta))) > 0
            assert support_margin(planes(robust + delta)) < 0
        output.append(result)
    with (out / "joint_dependence_thresholds.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)
    assert all(sha(Path(p)) == h for p, h in protected.items())
    report = dict(rows=len(output), analytical_controls=analytical, input_sha256=protected,
                  protected_inputs_unchanged=len(protected), maximum_independent_baseline_difference=max(baseline_differences),
                  classifications={name: {s: sum(r[name + "_status"] == s for r in output)
                                          for s in sorted({r[name + "_status"] for r in output})}
                                   for name in ("positive_rank", "negative_rank", "full_support")},
                  new_optimization_solves=0, independent_verification_pending=True,
                  scope="Declared triangular-marginal rank stresses and whole-rectangle guarantee, fixed policies; not empirical probability or extremal copula identification")
    (out / "joint_dependence_summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
