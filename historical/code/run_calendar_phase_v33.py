"""Partial synchronisation/unequal-size sensitivity. All demands are scenarios.

Uses optional explicit calendar demand; historical model/defaults are unchanged.
No private data, current EU operating probability or new background LCA.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import os
import time
from dataclasses import asdict
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from calendar_pool_network_v23 import CalendarLP, Scenario, cached_geometry, source_polytopes, dump_csv


def profile(weight_name, compression, origin, nodes=14, weeks=12):
    weights=np.ones(nodes) if weight_name=='uniform' else np.arange(1,nodes+1,dtype=float)
    if weight_name=='reverse_ramp':
        weights=weights[::-1]
    weights=weights/weights.sum()
    phases=origin+(1-compression)*2*np.pi*np.arange(nodes)/nodes
    demands=35000*weights[None,:]*(1+.5*np.sin(2*np.pi*np.arange(weeks)[:,None]/12+phases[None,:]))
    coherence=float(abs(np.sum(weights*np.exp(1j*phases))))
    return demands,weights,phases,coherence


def stock_bounds(demands,tau):
    rolling=np.asarray([demands[max(0,t-tau):t+1].sum(axis=0) for t in range(len(demands))])
    return float(np.max(np.sum(rolling,axis=1))),float(np.sum(np.max(rolling,axis=0)))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--package-root',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    a=ap.parse_args(); root=a.package_root.resolve(); out=a.output_dir.resolve()
    if False:
        raise ValueError('Work outputs must stay on D drive')
    out.mkdir(parents=True,exist_ok=True)
    assert (out/'analysis_protocol.md').is_file(),'Copy the declared protocol before execution'
    assert not (out/'phase_stock_screen.csv').exists(),'Fresh outputs required'
    began=time.monotonic()
    inputs=root/'inputs'
    source=source_polytopes(inputs,inputs/'replacement_stage_facts_v21.json')
    _,_,distance=cached_geometry(inputs/'geometry',out)
    climate=source['Climate change']
    minimum=linprog([1,1,0,0,0],A_ub=climate['A'],b_ub=climate['b'],bounds=[(None,None)]*5,method='highs')
    assert minimum.success
    nonnegative=[]
    for indicator,p in source.items():
        rest=linprog([0,0,1,1,1],A_ub=p['A'],b_ub=p['b'],bounds=[(None,None)]*5,method='highs')
        assert rest.success
        residual=float(rest.fun-p['process']['return_transport'][1])
        transport_min=p['process']['return_transport'][0]
        assert residual>0 and transport_min>0
        nonnegative.append(dict(indicator=indicator,remaining_stage_minimum=residual,transport_minimum=transport_min))
    rows=[]
    for weight in ['uniform','ramp','reverse_ramp']:
        for origin in [0.,math.pi/2,3*math.pi/2]:
            for tau in [1,3,5]:
                for step in range(21):
                    compression=step/20
                    d,w,p,R=profile(weight,compression,origin)
                    pooled,dedicated=stock_bounds(d,tau)
                    factor=float(minimum.fun/climate['u_upper_ratio']/d.sum())
                    rows.append(dict(weights=weight,phase_origin_radians=origin,compression=compression,
                                     return_delay_weeks=tau,phase_coherence=R,delivered_services=float(d.sum()),
                                     pooled_stock_lower_bound=pooled,dedicated_stock_lower_bound=dedicated,
                                     shareable_reserve_lower_bound=dedicated-pooled,
                                     relative_stock_gap=1-pooled/dedicated,
                                     pooled_climate_floor=pooled*factor,dedicated_climate_floor=dedicated*factor,
                                     pooled_all_routing_failure=pooled*factor>=1,
                                     dedicated_all_routing_failure=dedicated*factor>=1))
    dump_csv(out/'phase_stock_screen.csv',rows)
    outputs=[]; flows=[]; comparisons=[]
    for weight in ['uniform','ramp']:
        for tau in [1,3]:
            for compression in [0.,.5,.9,1.]:
                name=f'phase_{weight}_tau{tau}_c{compression:g}'
                d,w,p,R=profile(weight,compression,0.)
                bP,bD=stock_bounds(d,tau)
                pair={}
                for mode in ['dedicated_adaptive','pooled']:
                    s=Scenario(name,return_delay_weeks=tau,seasonal_amplitude=.5,seasonal_phase='explicit_phase_size')
                    lp=CalendarLP(s,distance,mode,demands=d)
                    result,x=lp.solve(source,robust=False)
                    result.update(scenario=asdict(s),mode=mode,
                                  phase_metadata=dict(weights=weight,compression=compression,phase_origin_radians=0.,
                                      node_weights=w.tolist(),node_phases=p.tolist(),phase_coherence=R),
                                  stock_lower_bound=bD if mode=='dedicated_adaptive' else bP)
                    outputs.append(result)
                    if result['success']:
                        assert result['counts']['initial_assets']>=result['stock_lower_bound']-1e-5
                        flows.extend(lp.schedule(x,name))
                    pair[mode]=result
                    print(f"{name} {mode} success={result['success']}",flush=True)
                D,P=pair['dedicated_adaptive'],pair['pooled']
                comparisons.append(dict(case=name,weights=weight,compression=compression,return_delay_weeks=tau,
                    phase_coherence=R,delivered_services=float(d.sum()),
                    dedicated_success=D['success'],pooled_success=P['success'],
                    dedicated_stock_lower_bound=bD,pooled_stock_lower_bound=bP,
                    dedicated_initial_stock=D['counts']['initial_assets'] if D['success'] else '',
                    pooled_initial_stock=P['counts']['initial_assets'] if P['success'] else '',
                    reference_partial_cost_reduction_fraction=1-P['cost_eur']/D['cost_eur'] if D['success'] and P['success'] else '',
                    dedicated_maximum_upper_ratio=D.get('maximum_upper_ratio',''),
                    pooled_maximum_upper_ratio=P.get('maximum_upper_ratio',''),
                    dedicated_selected_policy_outcome=D.get('environmental_outcome','not_solved'),
                    pooled_selected_policy_outcome=P.get('environmental_outcome','not_solved')))
    dump_csv(out/'calendar_positive_flows.csv',flows)
    dump_csv(out/'phase_network_comparisons.csv',comparisons)
    (out/'calendar_network_results.json').write_text(json.dumps(dict(outputs=outputs,
        scope='Declared phase and size sensitivity, all policies retained; source accounting, not measured operator histories'),indent=2)+'\n',encoding='utf-8')
    summary=dict(stock_screen_rows=len(rows),new_nominal_lp_solves=len(outputs),feasible_policies=sum(r['success'] for r in outputs),
                 paired_cases=len(comparisons),source_nonnegative_checks=nonnegative,
                 minimum_climate_initial_terminal_coefficient=float(minimum.fun),
                 bounds_only_screen_no_constructive_claim=True,all_results_retained=True,
                 empirical_demand_or_current_LCC=False,elapsed_seconds=time.monotonic()-began,
                 independent_verification_pending=True)
    (out/'phase_analysis_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
