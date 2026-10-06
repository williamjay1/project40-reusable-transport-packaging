"""Recompute joint budgets from saved reoptimized schedules, not a new LP run."""
import argparse
import json
from pathlib import Path
import numpy as np
from calendar_pool_network_v23 import source_polytopes,dump_csv
from calendar_financial_completion_v29 import upper_prior,joint_budget

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input-dir',type=Path,required=True);ap.add_argument('--stage-facts',type=Path,required=True)
    ap.add_argument('--results-dir',type=Path,required=True);a=ap.parse_args();root=a.results_dir
    meta=json.loads((root/'calendar_network_results.json').read_text());results=meta['outputs'];source=source_polytopes(a.input_dir,a.stage_facts)
    for r in results:
        if 'dual_objective_eur' in r:
            s=r['scenario'];terminal=(1+r['annual_discount_rate'])**(-(s['weeks']+s['return_delay_weeks']+1)/52)
            r['present_cost_dual_objective_eur']=r.pop('dual_objective_eur')+r['terminal_net_eur_per_asset']*terminal*r['counts']['missing_assets']
    prices=json.loads((root/'financial_input_audit.json').read_text())['ranges']['nrg_pc_205'];rng=np.random.default_rng(400023);le=rng.triangular(1,1,2,50000);lh=rng.triangular(1,1,2,50000);rows=[]
    for case in sorted({r['financial_case_id'] for r in results}):
        D=next(r for r in results if r['financial_case_id']==case and r['mode']=='dedicated_adaptive');P=next(r for r in results if r['financial_case_id']==case and r['mode']=='pooled')
        gap=D['included_present_cost_eur']-P['included_present_cost_eur'];transport_gap=D['present_components_eur']['empty_transport']-P['present_components_eur']['empty_transport']
        contexts=[(gamma,rho) for gamma in [1.,3.] for rho in ([0.,1.] if P['counts']['missing_assets'] else [0.])]
        for gamma,rho in contexts:
            env=upper_prior(P,source,le,lh,gamma,rho); mask=env<1
            for pe in [prices['lower'],prices['upper']]:
                for ph in [0.,.10,.20]:
                    unit=(le-1)*.0248*pe+(lh-1)*(.322/3.6)*ph
                    savings=gap+(gamma-1)*transport_gap-unit*(P['discounted_wash_equivalents']-D['discounted_wash_equivalents'])
                    budget,index,coverage=joint_budget(env,savings)
                    pool_cost=P['included_present_cost_eur']+(gamma-1)*P['present_components_eur']['empty_transport']+unit*P['discounted_wash_equivalents']
                    SU=float(np.sort(pool_cost[mask])[47499])/P['discounted_service_equivalents'] if mask.sum()>=47500 else None
                    rows.append(dict(case=case,parent_scenario=P['parent_scenario'],annual_discount_rate=P['annual_discount_rate'],terminal_net_eur_per_asset=P['terminal_net_eur_per_asset'],travel_factor=gamma,missing_fate_virgin_equivalent_rho=rho,electricity_price_context_eur_per_kwh=pe,declared_effective_heat_price_eur_per_kwh=ph,
                        environmental_prior_coverage=float(mask.mean()),joint_env_and_pool_cost_coverage_at_zero_extra_cost=coverage,
                        joint95_net_pooling_present_cost_supremum_eur=budget,joint95_budget_eur_per_discounted_service=None if budget is None else budget/P['discounted_service_equivalents'],budget_strictly_below_supremum=True,
                        joint95_required_singleuse_cost_eur_per_discounted_service=SU,SU_missing_cost_rule='add K_P/discounted_services; no observed single-use tariff claimed',draws=50000,seed=400023,cost_prior_link='same energy-input draws as pooled all-14-indicator environmental coverage; not product of marginal probabilities'))
    dump_csv(root/'joint_environment_financial_budgets.csv',rows)
    meta['monetary_indifference']='Absolute campaign savings below 1e-6 EUR classified as zero uniformly; no strict financial advantage at numerical equality.'
    meta['joint_budget_regeneration']='Same saved reoptimized schedules; joint evaluations include correlated missing-fate stress rho=0/1 when losses occur and uniformly suppress roundoff sign in exact flat control. Physical LP outputs not changed.'
    meta['joint_budget_rows']=len(rows)
    (root/'calendar_network_results.json').write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(recomputed_joint_rows=len(rows),reoptimized_again=False,roundoff_control='1e-6 EUR campaign neutrality'),indent=2))

if __name__=='__main__': main()
