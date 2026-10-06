"""Public-input collection/fate frontier and cost-structure stress scenarios.

Missing-fate allowance is an inverse budget, NOT a disposal inventory. A rho=1
stress adds one source virgin-production-equivalent burden per lost asset. No
observed loss probability or treatment pathway is inferred from these scenarios.
"""
import argparse
import json
from dataclasses import asdict
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from calendar_pool_network_v23 import Scenario, CalendarLP, cached_geometry, source_polytopes, dump_csv, calibrate


def evaluate_fate(result, source, le, lh, rho):
    c=result['counts']; s=result['scenario']; N=c['delivered_services']; W=c['wash_inputs']/100; L=c['missing_assets']
    D=(c['empty_dirty_crate_km']+s['clean_leg_factor']*c['empty_clean_crate_km'])/(139*100)
    vector=np.array([(c['virgin_assets']+rho*L)/100,c['terminal_recovered_assets']/100,
                     c['closed_loop_rejections']/.55,c['closed_loop_rejections']/.55,W])
    ratios=[]
    for k,p in source.items():
        lp=linprog(-vector,A_ub=p['A'],b_ub=p['b'],bounds=[(None,None)]*5,method='highs')
        assert lp.success
        base=vector@lp.x
        tt=p['process']['return_transport']; ee=p['process']['electricity']; hh=p['process']['water_heating']
        direct=base+max((D-W)*tt[0],(D-W)*tt[1])+np.maximum((le-1)*W*ee[0],(le-1)*W*ee[1])+np.maximum((lh-1)*W*hh[0],(lh-1)*W*hh[1])
        ratios.append(direct/(N/100*p['u_lower_ratio']))
    maximum=np.max(np.asarray(ratios),axis=0)
    return float(np.mean(maximum<1)),float(np.quantile(maximum,.95)),float(np.max(maximum))


def main():
    ap=argparse.ArgumentParser()
    for name in ['input-dir','stage-facts','geometry-cache','output-dir']:
        ap.add_argument('--'+name,type=Path,required=True)
    a=ap.parse_args(); out=a.output_dir; out.mkdir(parents=True,exist_ok=True)
    _,_,dist=cached_geometry(a.geometry_cache,out); source=source_polytopes(a.input_dir,a.stage_facts)
    dump_csv(out/'source_calendar_calibration.csv',calibrate(source,a.input_dir))
    rng=np.random.default_rng(400025); le=rng.triangular(1,1,2,50000); lh=rng.triangular(1,1,2,50000)
    outputs=[]; flows=[]; frontier=[]; allowances=[]
    for tau in [1,3]:
        for q in np.arange(70,101)/100:
            s=Scenario(f'staggered_tau{tau}_q{q:.2f}',return_delay_weeks=tau,collection_fraction=float(q),seasonal_amplitude=.5,seasonal_phase='staggered')
            model=CalendarLP(s,dist,'pooled'); r,x=model.solve(source)
            r.update(scenario=asdict(s),mode='pooled',environmental_constraints=False)
            outputs.append(r)
            if not r['success']: continue
            flows+=model.schedule(x,s.name)
            for rho in [0.,1.]:
                prob,tail,mx=evaluate_fate(r,source,le,lh,rho)
                _,rect,_=evaluate_fate(r,source,np.array([2.]),np.array([2.]),rho)
                frontier.append(dict(scenario=s.name,return_delay_weeks=tau,collection_fraction=q,
                                     missing_fate_virgin_production_equivalent_stress=rho,design_prior_probability=prob,
                                     design_prior_95pct_maximum_ratio=tail,full_energy_rectangle_maximum_ratio=rect,
                                     reference_cost_per_service=r['cost_eur_per_service'],source_scope='conditional engineering screen, not fitted EU collection probability'))
            L=r['counts']['missing_assets']
            if L>0:
                N=r['counts']['delivered_services']
                for k,p in source.items():
                    margin=(N/100*p['u_lower_ratio'])*(1-r['impact_ratios'][k]['upper'])*p['u']
                    allowances.append(dict(scenario=s.name,indicator=k,collection_fraction=q,missing_assets=L,
                                           maximum_extra_missing_fate_burden_per_lost_asset=margin/L,
                                           source_units='same indicator units as Tua supplement; negative means no nonnegative allowance',
                                           calculated_not_observed=True))
            print(s.name,r['maximum_upper_ratio'],flush=True)
    # Predeclared adverse financial/capacity structures; never pick only winners.
    stresses=[('transport_dominated',dict(manufacturing_eur_per_crate=2.275,cost_per_truck_km=2.,holding_eur_per_crate_week=0.)),
              ('asset_dominated',dict(manufacturing_eur_per_crate=9.10,cost_per_truck_km=.125,holding_eur_per_crate_week=.09)),
              ('capacity_tight',dict(wash_capacity_ratio=1.)),
              ('aligned_season',dict(seasonal_phase='aligned'))]
    comparisons=[]
    for label,changes in stresses:
        base=dict(seasonal_amplitude=.5,seasonal_phase='staggered'); base.update(changes)
        s=Scenario(label,**base); pair={}
        for mode in ['dedicated_adaptive','pooled']:
            model=CalendarLP(s,dist,mode); r,x=model.solve(source)
            r.update(scenario=asdict(s),mode=mode,environmental_constraints=False); outputs.append(r)
            if x is not None: flows+=model.schedule(x,s.name)
            pair[mode]=r
        d,p=pair['dedicated_adaptive'],pair['pooled']
        comparisons.append(dict(scenario=label,dedicated_success=d['success'],pooled_success=p['success'],
                                cost_reduction_fraction=1-p['cost_eur']/d['cost_eur'] if d['success'] and p['success'] else None,
                                stock_reduction_fraction=1-p['counts']['initial_assets']/d['counts']['initial_assets'] if d['success'] and p['success'] else None,
                                dedicated_maximum_upper_ratio=d.get('maximum_upper_ratio'),pooled_maximum_upper_ratio=p.get('maximum_upper_ratio')))
        print(label,comparisons[-1],flush=True)
    dump_csv(out/'calendar_positive_flows.csv',flows)
    dump_csv(out/'collection_energy_fate_frontier.csv',frontier)
    dump_csv(out/'missing_fate_allowances.csv',allowances)
    dump_csv(out/'financial_capacity_stress_comparisons.csv',comparisons)
    thresholds=[]
    for tau in [1,3]:
        for rho in [0.,1.]:
            subset=[r for r in frontier if r['return_delay_weeks']==tau and r['missing_fate_virgin_production_equivalent_stress']==rho]
            feasible=[r for r in subset if r['design_prior_probability']>=.95]
            probabilities=[r['design_prior_probability'] for r in subset]
            monotone=bool(np.all(np.diff(probabilities)>=-1e-12))
            q=min(r['collection_fraction'] for r in feasible) if feasible else None
            thresholds.append(dict(return_delay_weeks=tau,missing_fate_stress=rho,first_tested_collection_fraction_at_95pct_design_coverage=q,
                                   lower_bracket=q-.01 if q is not None and q>.70 else None,probability_monotone_on_tested_grid=monotone,
                                   exact_continuous_minimum_not_claimed=True))
    (out/'collection_frontier_summary.json').write_text(json.dumps(dict(seed=400025,draws=50000,prior='independent triangular(1,1,2) electricity and heat input multipliers',
            missing_fate_definition='rho=1 is one source virgin-production-equivalent burden per lost crate, not an actual waste inventory; rho=0 retains no missing recycling credit but no additional unknown fate charge',
            thresholds=thresholds,cost_stresses=comparisons),indent=2)+'\n',encoding='utf-8')
    (out/'calendar_network_results.json').write_text(json.dumps(dict(scope='conditional collection/fate and predeclared financial-capacity stress scenarios; no empirical EU rate',outputs=outputs),indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(schedules=len(outputs),thresholds=thresholds),indent=2),flush=True)


if __name__=='__main__': main()
