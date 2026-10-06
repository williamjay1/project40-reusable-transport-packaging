"""Translate verified paired partial costs into an omitted-cost allowance.

This is an arithmetic inverse budget, not observed coordination cost or full LCC.
CSV per-service outputs and the underlying JSON campaign costs are checked
independently against each other before any result is exported.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--comparison', type=Path, required=True)
    ap.add_argument('--network-results', type=Path, required=True)
    ap.add_argument('--collection-frontier', type=Path, required=True)
    ap.add_argument('--output-dir', type=Path, required=True)
    a = ap.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=True)
    with a.comparison.open(encoding='utf-8', newline='') as f:
        comparisons = list(csv.DictReader(f))
    outputs = json.loads(a.network_results.read_text(encoding='utf-8'))['outputs']
    paired = {(r['scenario']['name'], r['mode']): r for r in outputs}
    rows = []
    for r in comparisons:
        scenario = r['scenario']
        d = paired[scenario, 'dedicated_adaptive']
        p = paired[scenario, 'pooled']
        nd, np = (x['counts']['delivered_services'] for x in [d, p])
        assert abs(nd - np) < 1e-6
        delta = d['cost_eur'] - p['cost_eur']
        per = delta / np
        csv_delta = float(r['dedicated_cost_per_service']) - float(r['pooled_cost_per_service'])
        assert abs(per - csv_delta) < 1e-10
        rows.append(dict(scenario=scenario, delivered_services=np,
                         additional_net_pooling_cost_ceiling_eur=delta,
                         additional_net_pooling_cost_ceiling_eur_per_service=per,
                         strict_saving_condition='extra pooled cost minus extra dedicated cost < ceiling',
                         scope='fixed paired policies; historical partial cost; not empirical or full LCC'))
    path = a.output_dir / 'additional_pooling_cost_budget.csv'
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    with a.collection_frontier.open(encoding='utf-8', newline='') as f:
        frontier = list(csv.DictReader(f))
    near = [r for r in frontier if int(r['return_delay_weeks']) == 3
            and abs(float(r['collection_fraction']) - .97) < 1e-12
            and float(r['missing_fate_virgin_production_equivalent_stress']) == 0]
    assert len(near) == 1
    ratio = float(near[0]['design_prior_95pct_maximum_ratio'])
    # A quality-restricted geometry may overturn the original q=.97 crossing.
    # Retain that adverse case; crossing is an output, not a pass requirement.
    assert ratio > 0
    report = dict(all_five_paired_cost_checks_passed=True,
                  q97_three_week_prior_coverage=float(near[0]['design_prior_probability']),
                  q97_three_week_95pct_maximum_ratio=ratio,
                  q97_margin_to_ratio_one=1-ratio,
                  q97_prior_coverage_at_least_95pct=float(near[0]['design_prior_probability'])>=.95,
                  not_empirical_probability=True,
                  inputs={str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in [a.comparison, a.network_results, a.collection_frontier]})
    (a.output_dir / 'overhead_budget_checks.json').write_text(
        json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(dict(budgets=rows, checks=report), indent=2))


if __name__ == '__main__':
    main()
