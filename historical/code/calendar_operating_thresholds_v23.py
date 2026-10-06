"""Fixed-policy inverse travel budgets and declared-design-prior thresholds.

No distribution below is fitted to European observations. The reported 95%
threshold is coverage under the explicitly declared engineering prior, not an
empirical EU confidence limit or a claim about regulation. LCIA printed rounding
uses correlated set bounds, never Monte Carlo pseudo-observations.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from calendar_pool_network_v23 import source_polytopes, dump_csv

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--input-dir',type=Path,required=True)
    ap.add_argument('--stage-facts',type=Path,required=True)
    ap.add_argument('--network-results',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    args=ap.parse_args(); out=args.output_dir; out.mkdir(parents=True,exist_ok=True)
    source=source_polytopes(args.input_dir,args.stage_facts)
    results=json.loads(args.network_results.read_text(encoding='utf-8'))['outputs']
    seed=400023; draws=50000; rng=np.random.default_rng(seed)
    # Engineering-prior declaration committed in code before execution. Common
    # random numbers compare all schedules; no evidence-selected ranges.
    # Conservative design amendment: do not admit unverified energy cuts below
    # the source cleaning bundle. These distributions cover added demands only.
    le=rng.triangular(1.,1.,2.,draws); lh=rng.triangular(1.,1.,2.,draws)
    threshold_rows=[]; proof_rows=[]; certificates=[]
    for k,p in source.items():
        virgin_terminal=linprog([1,1,0,0,0],A_ub=p['A'],b_ub=p['b'],bounds=[(None,None)]*5,method='highs')
        rest=linprog([0,0,1,1,1],A_ub=p['A'],b_ub=p['b'],bounds=[(None,None)]*5,method='highs')
        assert virgin_terminal.success and rest.success
        residual=float(rest.fun-p['process']['return_transport'][1])
        assert residual>0,k
        proof_rows.append(dict(indicator=k,initial_net_ratio_lower=float(virgin_terminal.fun),
                               remaining_net_reconditioning_ratio_lower=residual))
    for result in results:
        if not result['success'] or result['mode']=='dedicated_nearest': continue
        c=result['counts']; s=result['scenario']; N=c['delivered_services']
        W=c['wash_inputs']/100
        D=(c['empty_dirty_crate_km']+s['clean_leg_factor']*c['empty_clean_crate_km'])/(139*100)
        vec=np.array([c['virgin_assets']/100,c['terminal_recovered_assets']/100,
                      c['closed_loop_rejections']/.55,c['closed_loop_rejections']/.55,W])
        per_indicator=[]; nominal=[]; positively_paired=[]; negatively_paired=[]; interval_minima=[]
        for k,p in source.items():
            sol=linprog(-vec,A_ub=p['A'],b_ub=p['b'],bounds=[(None,None)]*5,method='highs')
            assert sol.success
            base=float(vec@sol.x); target=N/100*p['u_lower_ratio']
            tl,th=p['process']['return_transport']; el,eh=p['process']['electricity']; hl,hh=p['process']['water_heating']
            energy=np.maximum((le-1)*W*el,(le-1)*W*eh)+np.maximum((lh-1)*W*hl,(lh-1)*W*hh)
            # max t*(gamma D-W) is piecewise affine; use the matching endpoint.
            high=(target-base-energy+th*W)/(th*D)
            low=(target-base-energy+tl*W)/(tl*D)
            root=np.where(high>=W/D,high,low)
            def inverse(ee,hh_):
                z=np.maximum((ee-1)*W*el,(ee-1)*W*eh)+np.maximum((hh_-1)*W*hl,(hh_-1)*W*hh)
                hi=(target-base-z+th*W)/(th*D); lo=(target-base-z+tl*W)/(tl*D)
                return np.where(hi>=W/D,hi,lo)
            positively_paired.append(inverse(np.sort(le),np.sort(lh)))
            negatively_paired.append(inverse(np.sort(le),np.sort(lh)[::-1]))
            interval_minima.append(float(inverse(2.,2.)))
            nh=(target-base+th*W)/(th*D); nl=(target-base+tl*W)/(tl*D)
            nominal.append((k,float(nh if nh>=W/D else nl)))
            per_indicator.append(root)
        roots=np.min(np.asarray(per_indicator),axis=0)
        positive_roots=np.min(np.asarray(positively_paired),axis=0)
        negative_roots=np.min(np.asarray(negatively_paired),axis=0)
        limiting,nominalroot=min(nominal,key=lambda v:v[1])
        threshold_rows.append(dict(scenario=s['name'],mode=result['mode'],
                                   fixed_policy_travel_multiplier_nominal=nominalroot,
                                   nominal_limiting_indicator=limiting,
                                   travel_multiplier_design_prior_95_coverage=float(np.quantile(roots,.05)),
                                   travel_multiplier_design_prior_median=float(np.median(roots)),
                                   probability_at_reported_distances=float(np.mean(roots>1)),
                                   probability_at_double_distances=float(np.mean(roots>2)),
                                   positive_rank_dependence_prior_95_coverage=float(np.quantile(positive_roots,.05)),
                                   negative_rank_dependence_prior_95_coverage=float(np.quantile(negative_roots,.05)),
                                   robust_travel_ceiling_for_full_energy_interval=min(interval_minima),
                                   quantile_convergence_5000=float(np.quantile(roots[:5000],.05)),
                                   quantile_convergence_10000=float(np.quantile(roots[:10000],.05)),
                                   negative_threshold_prior_fraction=float(np.mean(roots<=0)),
                                   draws=draws,declared_prior_not_empirical=True,
                                   reference_cost_eur_per_service=result['cost_eur_per_service'],
                                   environmental_scope='conditional lost-fate screen' if c['missing_assets'] else 'source-accounting terminal liquidation'))
        if c['missing_assets']==0:
            dem=np.array(result['calendar_demands']); ell=s['return_delay_weeks']+1
            rolling=np.array([dem[max(0,t-ell+1):t+1].sum(axis=0) for t in range(s['weeks'])])
            owner_dedicated=result['mode'].startswith('dedicated')
            opening_bound=float(rolling.max(axis=0).sum() if owner_dedicated else rolling.sum(axis=1).max())
            assert c['initial_assets']>=opening_bound-1e-5
            lower={r['indicator']:opening_bound/100*r['initial_net_ratio_lower']/(N/100*source[r['indicator']]['u_upper_ratio']) for r in proof_rows}
            certificates.append(dict(scenario=s['name'],mode=result['mode'],occupancy_stock_lower_bound=opening_bound,
                                     bound_type='sum of owner-specific rolling maxima' if owner_dedicated else 'maximum of aggregate rolling occupancy',
                                     actual_initial_stock=c['initial_assets'],minimum_opening_stock_impact_ratios=lower,
                                     source_conditional_global_failure=any(v>=1 for v in lower.values()),
                                     reason='source net manufacturing+terminal charge at mandatory occupied stock, even before nonnegative reconditioning and transport'))
    comparisons=[]
    names=sorted({r['scenario']['name'] for r in results if r['success']})
    for name in names:
        fixed=next(r for r in results if r['success'] and r['scenario']['name']==name and r['mode']==('dedicated_adaptive' if any(v['mode']=='dedicated_adaptive' for v in results) else 'dedicated_static'))
        pooled=next(r for r in results if r['success'] and r['scenario']['name']==name and r['mode']=='pooled')
        comparisons.append(dict(scenario=name,reference_cost_reduction_fraction=1-pooled['cost_eur']/fixed['cost_eur'],
                                initial_stock_reduction_fraction=1-pooled['counts']['initial_assets']/fixed['counts']['initial_assets'],
                                dedicated_initial_assets=fixed['counts']['initial_assets'],pooled_initial_assets=pooled['counts']['initial_assets'],
                                dedicated_cost_per_service=fixed['cost_eur_per_service'],pooled_cost_per_service=pooled['cost_eur_per_service'],
                                dedicated_maximum_upper_impact_ratio=fixed['maximum_upper_ratio'],pooled_maximum_upper_impact_ratio=pooled['maximum_upper_ratio'],
                                benchmark_mode=fixed['mode']))
    dump_csv(out/'fixed_policy_probabilistic_travel_thresholds.csv',threshold_rows)
    dump_csv(out/'strong_baseline_comparisons.csv',comparisons)
    dump_csv(out/'source_nonnegative_remaining_stage_checks.csv',proof_rows)
    (out/'occupancy_failure_certificates.json').write_text(json.dumps(certificates,indent=2)+'\n',encoding='utf-8')
    meta=dict(prior='Independent triangular(1,1,2) electricity/heat INPUT-QUANTITY multipliers on source processes; source energy is a hard lower bound, with no unverified hygiene-related reductions. These are declared design weights, not estimated EU distributions.',
              correlation_limit='Changing their dependence can alter probabilistic thresholds; robust interval and nominal results do not require this prior.',
              seed=seed,draws=draws,fixed_policy_no_reoptimization=True,
              financial_interpretation='Reference partial lifecycle operating/procurement cost, equivalently required single-use cost benchmark; no observed single-use tariff or complete LCC asserted.',
              travel_interpretation='Both specified empty route lengths multiplied uniformly under unchanged source vehicle/loading process. Not a maximum geographic network diameter.',
              threshold_rows=threshold_rows,comparisons=comparisons)
    (out/'operating_threshold_summary.json').write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(meta,indent=2))

if __name__=='__main__': main()
