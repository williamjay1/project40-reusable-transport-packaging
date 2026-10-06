"""Necessary periodic-stock certificate, separate from selected-policy roots."""
import argparse
import json
from pathlib import Path
from scipy.optimize import linprog
from calendar_pool_network_v23 import source_polytopes, dump_csv


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--package-root',type=Path,required=True);ap.add_argument('--results-dir',type=Path,required=True)
    a=ap.parse_args();inputs=a.package_root/'inputs';out=a.results_dir
    assert not (out/'periodic_stock_certificates.csv').exists()
    source=source_polytopes(inputs,inputs/'replacement_stage_facts_v21.json')
    signs={}
    for k,p in source.items():
        vals={}
        for tag,vec in [('capital',[1,1,0,0,0]),('virgin',[1,0,0,0,0]),('remaining',[0,0,1,1,1])]:
            sol=linprog(vec,A_ub=p['A'],b_ub=p['b'],bounds=[(None,None)]*5,method='highs')
            assert sol.success;vals[tag]=float(sol.fun)
        vals['remaining_excluding_reference_return']=vals['remaining']-p['process']['return_transport'][1]
        vals['transport_minimum']=p['process']['return_transport'][0]
        assert vals['capital']>0 and vals['virgin']>0 and vals['remaining_excluding_reference_return']>0 and vals['transport_minimum']>0
        signs[k]=vals
    climate=source['Climate change'];coef=signs['Climate change']['capital']/climate['u_upper_ratio']
    policies=json.loads((out/'periodic_results.json').read_text())['outputs'];rows=[]
    for r in policies:
        if not r['success']: continue
        N=r['counts']['delivered_services'];LB=r['periodic_fleet_lower_bound'];H=r['period_weeks']
        beta=N/(H*LB*coef)
        rows.append(dict(case=r['scenario']['name'],mode=r['mode'],periodic_fleet_lower_bound=LB,
                         climate_capital_coefficient_favourable=coef,
                         beta_all_routing_climate_failure_from=beta,
                         selected_beta_supremum=r['inverse_stock_charge']['beta_supremum'],
                         inverse_stock_floor_not_selected_policy_failure=True))
    dump_csv(out/'periodic_stock_certificates.csv',rows)
    (out/'periodic_stock_signs.json').write_text(json.dumps(dict(signs=signs,
        certificate='For beta >= N/(H*LB*kappa_min) the favourable climate capital floor is >=1. Independently positive remaining, missing-virgin and actual transport terms cannot rescue any route. This is a sufficient all-routing failure condition, not the exact global break-even.'),indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(periodic_stock_certificates=len(rows),indicators_signs_checked=len(signs)),indent=2))


if __name__=='__main__':main()
