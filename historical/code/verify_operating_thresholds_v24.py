"""Independent numerical inversion audit; never imports a model generator.

Rebuilds source rounding constraints from numeric facts, evaluates the direct
event ledger, and uses bisection rather than the producer's algebraic roots.
Declared priors are not fitted European probabilities.
"""
import argparse
import csv
import json
from fractions import Fraction
from pathlib import Path
import numpy as np
from scipy.optimize import linprog


def read_csv(path):
    with path.open(newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))


def reconstruct(grid, facts, processes):
    fixed = dict(zip(['Climate change', 'Photochemical ozone formation', 'Acidification',
                     'Terrestrial eutrophication', 'Freshwater ecotoxicity',
                     'Cumulative energy demand', 'Water depletion'], [.01,.01,.01,.01,1.,1.,.01]))
    order = ['Climate change','Ozone depletion','Human toxicity, non-cancer effects',
             'Human toxicity, cancer effects','Particulate matter','Photochemical ozone formation',
             'Acidification','Terrestrial eutrophication','Freshwater eutrophication',
             'Marine eutrophication','Freshwater ecotoxicity',
             'Mineral, fossil & renewable resource depletion','Cumulative energy demand','Water depletion']
    precision = dict(zip(order,[.01,1e-7,1e-7,1e-7,1e-4,1e-3,1e-3,1e-3,1e-4,1e-4,.01,1e-5,1.,.01]))
    polys={}
    for f in facts:
        k=f['indicator']; gs=[g for g in grid if g['indicator']==k and g['source_scenario']=='W1_RE1']
        u=float(gs[0]['single_use_impact_per_service_100_crates']); A=[]; b=[]
        def add(v,y,h):
            A.extend([v,[-x for x in v]]); b.extend([(y+h)/u,(-y+h)/u])
        for g in gs:
            n=int(g['reported_rotations']); y=float(g['reusable_cumulative_impact_100_crates'])
            h=fixed[k]/2 if k in fixed else 10**(np.floor(np.log10(y))-2)/2
            add([1,1,n-1,n-1,n-1],y,h)
        for n in [20,125]:
            for tag,v in [('P',[1,0,n-1,0,0]),('W',[0,0,0,0,n-1]),('E',[0,1,0,n-1,0])]:
                cell=f['cells'][f'{tag}{n}']; add(v,cell['value'],(cell['upper']-cell['lower'])/2)
        proc={}
        for p in [p for p in processes if p['indicator']==k]:
            m,_,ex=p['displayed_value'].upper().partition('E'); e=int(ex) if ex else 0
            h=Fraction(1,2)*Fraction(10)**(e-len(m.partition('.')[2])); y=Fraction(m)*Fraction(10)**e
            proc[p['process']]=np.array([float(y-h)/u,float(y+h)/u])
        polys[k]=(np.asarray(A),np.asarray(b),1-precision[k]/(2*u),proc)
    return polys


def main():
    ap=argparse.ArgumentParser()
    for name in ['source-grid','stage-facts','process-facts','network-results','thresholds','output']:
        ap.add_argument('--'+name,type=Path,required=True)
    a=ap.parse_args()
    polys=reconstruct(read_csv(a.source_grid),json.loads(a.stage_facts.read_text())['facts'],read_csv(a.process_facts))
    results=json.loads(a.network_results.read_text())['outputs']; table=read_csv(a.thresholds)
    rng=np.random.default_rng(400023); le=rng.triangular(1,1,2,50000); lh=rng.triangular(1,1,2,50000)
    checks=[]
    for row in table:
        r=next(r for r in results if r['success'] and r['mode']==row['mode'] and r['scenario']['name']==row['scenario'])
        c=r['counts']; s=r['scenario']; N=c['delivered_services']; W=c['wash_inputs']/100
        D=(c['empty_dirty_crate_km']+s['clean_leg_factor']*c['empty_clean_crate_km'])/(139*100)
        v=np.array([c['virgin_assets']/100,c['terminal_recovered_assets']/100,
                    c['closed_loop_rejections']/.55,c['closed_loop_rejections']/.55,W])
        support=[]
        for k,(A,b,u,proc) in polys.items():
            lp=linprog(-v,A_ub=A,b_ub=b,bounds=[(None,None)]*5,method='highs')
            assert lp.success
            assert np.all(proc['return_transport']>0)
            support.append((float(v@lp.x),N/100*u,proc))
        def residual(gamma,ee,hh):
            vals=[]
            for base,target,p in support:
                tr=np.maximum((gamma*D-W)*p['return_transport'][0],(gamma*D-W)*p['return_transport'][1])
                e=np.maximum((ee-1)*W*p['electricity'][0],(ee-1)*W*p['electricity'][1])
                h=np.maximum((hh-1)*W*p['water_heating'][0],(hh-1)*W*p['water_heating'][1])
                vals.append(base+tr+e+h-target)
            return np.max(np.asarray(vals),axis=0)
        def invert(ee,hh):
            shape=np.broadcast_arrays(np.asarray(ee),np.asarray(hh))[0].shape
            lo=np.full(shape,-10000.); hi=np.full(shape,10000.)
            assert np.all(residual(lo,ee,hh)<0) and np.all(residual(hi,ee,hh)>0)
            for _ in range(48):
                mid=(lo+hi)/2; feasible=residual(mid,ee,hh)<0
                lo=np.where(feasible,mid,lo); hi=np.where(feasible,hi,mid)
            return (lo+hi)/2
        samples=invert(le,lh); pos=invert(np.sort(le),np.sort(lh)); neg=invert(np.sort(le),np.sort(lh)[::-1])
        values={
            'fixed_policy_travel_multiplier_nominal':float(invert(1.,1.)),
            'travel_multiplier_design_prior_95_coverage':float(np.quantile(samples,.05)),
            'travel_multiplier_design_prior_median':float(np.median(samples)),
            'positive_rank_dependence_prior_95_coverage':float(np.quantile(pos,.05)),
            'negative_rank_dependence_prior_95_coverage':float(np.quantile(neg,.05)),
            'robust_travel_ceiling_for_full_energy_interval':float(invert(2.,2.)),
            'probability_at_reported_distances':float(np.mean(residual(1.,le,lh)<0)),
            'probability_at_double_distances':float(np.mean(residual(2.,le,lh)<0)),
            'negative_threshold_prior_fraction':float(np.mean(samples<=0)),
            'quantile_convergence_5000':float(np.quantile(samples[:5000],.05)),
            'quantile_convergence_10000':float(np.quantile(samples[:10000],.05))}
        errors={k:abs(val-float(row[k])) for k,val in values.items()}
        assert max(errors.values())<1e-7,(row['scenario'],row['mode'],errors)
        # A positive reported quantile must indeed be a physical coverage budget.
        gamma=values['travel_multiplier_design_prior_95_coverage']
        coverage=float(np.mean(residual(gamma,le,lh)<0)) if gamma>=0 else None
        if coverage is not None: assert abs(coverage-.95)<.0001
        checks.append({'scenario':row['scenario'],'mode':row['mode'],'checked_statistics':len(values),
                       'maximum_absolute_difference':max(errors.values()),'empirical_design_sample_coverage':coverage})
    occupancy_path=a.thresholds.parent/'occupancy_failure_certificates.json'
    certificates=json.loads(occupancy_path.read_text())
    for cert in certificates:
        r=next(r for r in results if r['success'] and r['mode']==cert['mode'] and r['scenario']['name']==cert['scenario'])
        d=np.asarray(r['calendar_demands']); tau=r['scenario']['return_delay_weeks']
        rolling=np.array([d[max(0,t-tau):t+1].sum(axis=0) for t in range(len(d))])
        bound=float(rolling.max(axis=0).sum() if r['mode'].startswith('dedicated') else rolling.sum(axis=1).max())
        assert abs(bound-cert['occupancy_stock_lower_bound'])<1e-6
        floors=[]
        for k,(A,b,u,p) in polys.items():
            init=linprog([1,1,0,0,0],A_ub=A,b_ub=b,bounds=[(None,None)]*5,method='highs')
            remain=linprog([0,0,1,1,1],A_ub=A,b_ub=b,bounds=[(None,None)]*5,method='highs')
            assert init.success and remain.success and remain.fun-max(p['return_transport'])>0
            floor=bound*init.fun/(d.sum()*(2-u)); floors.append(floor)
            assert abs(floor-cert['minimum_opening_stock_impact_ratios'][k])<1e-8
        assert bool(max(floors)>=1)==cert['source_conditional_global_failure']
    report={'all_passed':True,'independent_method':'raw-fact source polytope plus direct-ledger 48-step bisection; no generator imports',
            'draws_per_schedule':50000,'schedules':len(checks),'checks':checks,'owner_aware_occupancy_certificates_checked':len(certificates),
            'interpretation':'conditional coverage under declared engineering priors, not empirical EU probabilities'}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__': main()
