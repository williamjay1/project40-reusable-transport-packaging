"""Independent event/cost ledger verifier. Does not import the generator.

Checks positive schedules, their physical balances, source interval support and
terminal closure. LP optimality is not inferred from this verifier; the reported
HiGHS solution and stronger benchmark comparison are separate evidence.
"""
import argparse
import csv
import json
from fractions import Fraction
from pathlib import Path
import numpy as np
from scipy.optimize import linprog

def rows(path):
    with path.open(newline='',encoding='utf-8-sig') as f: return list(csv.DictReader(f))

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--results-dir',type=Path,required=True)
    ap.add_argument('--source-grid',type=Path,required=True)
    ap.add_argument('--stage-facts',type=Path,required=True)
    ap.add_argument('--process-facts',type=Path,required=True)
    a=ap.parse_args(); root=a.results_dir
    results=json.loads((root/'calendar_network_results.json').read_text())
    meta=json.loads((root/'geometry_provenance.json').read_text())
    geo=rows(root/'geometry_numeric_facts.csv'); J=len(meta['sorters']); I=len(meta['washers'])
    di=np.array([[float(next(r['distance_km'] for r in geo if r['origin']==j and r['destination']==i))
                  for i in meta['washers']] for j in meta['sorters']])
    fs=rows(root/'calendar_positive_flows.csv'); allchecks=[]
    grid=rows(a.source_grid); stages=json.loads(a.stage_facts.read_text())['facts']; pf=rows(a.process_facts)
    fixed=dict(zip(['Climate change','Photochemical ozone formation','Acidification','Terrestrial eutrophication','Freshwater ecotoxicity','Cumulative energy demand','Water depletion'],[.01,.01,.01,.01,1.,1.,.01]))
    comparative=dict(zip(['Climate change','Ozone depletion','Human toxicity, non-cancer effects','Human toxicity, cancer effects','Particulate matter','Photochemical ozone formation','Acidification','Terrestrial eutrophication','Freshwater eutrophication','Marine eutrophication','Freshwater ecotoxicity','Mineral, fossil & renewable resource depletion','Cumulative energy demand','Water depletion'],[.01,1e-7,1e-7,1e-7,1e-4,1e-3,1e-3,1e-3,1e-4,1e-4,.01,1e-5,1.,.01]))
    polytopes={}
    for sr in stages:
        k=sr['indicator']; gs=[g for g in grid if g['indicator']==k and g['source_scenario']=='W1_RE1']
        u=float(gs[0]['single_use_impact_per_service_100_crates']); aa=[]; bb=[]
        def constraint(v,y,h):
            aa.extend([v,[-x for x in v]]); bb.extend([(y+h)/u,(-y+h)/u])
        for g in gs:
            n=int(g['reported_rotations']); y=float(g['reusable_cumulative_impact_100_crates'])
            h=fixed[k]/2 if k in fixed else 10**(np.floor(np.log10(y))-2)/2
            constraint([1,1,n-1,n-1,n-1],y,h)
        for n in [20,125]:
            for tag,v in [('P',[1,0,n-1,0,0]),('W',[0,0,0,0,n-1]),('E',[0,1,0,n-1,0])]:
                cell=sr['cells'][f'{tag}{n}']; constraint(v,cell['value'],(cell['upper']-cell['lower'])/2)
        proc={}
        for row in [p for p in pf if p['indicator']==k]:
            mant,_,ex=row['displayed_value'].upper().partition('E'); e=int(ex) if ex else 0
            half=Fraction(1,2)*Fraction(10)**(e-len(mant.partition('.')[2]))
            val=Fraction(mant)*Fraction(10)**e
            proc[row['process']]=(float(val-half)/u,float(val+half)/u)
        polytopes[k]=(np.array(aa),np.array(bb),u,proc)
    for result in results['outputs']:
        if not result['success']: continue
        s=result['scenario']; name=s['name']; mode=result['mode']; H=s['weeks']
        scope=[r for r in fs if r['scenario']==name and r['mode']==mode]
        f={ (r['variable'],)+tuple(int(v) for v in r['coordinates'].split(':')):float(r['value']) for r in scope}
        assert len(scope)==len(f)
        def get(tag,*coord): return f.get((tag,)+coord,0.)
        stocklen=len(result['asset_slots']); initial=sum(get('initial',slot) for slot in range(stocklen))
        d=np.array(result['calendar_demands']); tau=s['return_delay_weeks']; q=s['collection_fraction']
        maxres=0.; totalwash=0.; DD=0.; DS=0.; QW=0.; DW=0.; new=0.
        def check(value):
            nonlocal maxres
            maxres=max(maxres,abs(value))
        for t in range(H):
            for j in range(J):
                check(sum(get('dispatch',t,j,i) for i in range(I))-d[t,j])
                check(get('dirty',t,j)-(get('dirty',t-1,j) if t else 0.)+
                      sum(get('wash',t,j,i) for i in range(I))-(q*d[t-tau,j] if t>=tau else 0.))
                DW+=get('dirty',t,j)
                for i in range(I):
                    v=get('wash',t,j,i); totalwash+=v; DD+=v*di[j,i]
                    DS+=get('dispatch',t,j,i)*di[j,i]
            for slot in range(stocklen):
                rhs=get('ready',t-1,slot) if t else get('initial',slot)
                if mode=='dedicated_nearest':
                    outbound=sum(get('dispatch',t,slot,i) for i in range(I))
                    inbound=sum(get('wash',t-1,slot,i) for i in range(I)) if t else 0.
                elif mode in ['dedicated_static','dedicated_adaptive']:
                    j,i=result['asset_slots'][slot]
                    outbound=get('dispatch',t,j,i); inbound=get('wash',t-1,j,i) if t else 0.
                else:
                    outbound=sum(get('dispatch',t,j,slot) for j in range(J))
                    inbound=sum(get('wash',t-1,j,slot) for j in range(J)) if t else 0.
                check(get('ready',t,slot)-rhs-get('new',t,slot)-inbound+outbound)
                QW+=get('ready',t,slot); new+=get('new',t,slot)
            for i in range(I):
                used=sum(get('wash',t,j,i) for j in range(J))
                assert used<=s['mean_weekly_demand']*s['wash_capacity_ratio']/I+1e-6
            cumulative_new=sum(get('new',tt,slot) for tt in range(t+1) for slot in range(stocklen))
            cumulative_loss=(1-q)*d[:max(t-tau+1,0)].sum()
            check(cumulative_new-cumulative_loss)
            inuse=d[max(t-tau+1,0):t+1].sum() if tau else 0.
            active=sum(get('ready',t,slot) for slot in range(stocklen))+sum(get('dirty',t,j) for j in range(J))+sum(get('wash',t,j,i) for j in range(J) for i in range(I))+inuse
            check(active-(initial+cumulative_new-cumulative_loss))
        assert maxres<1e-5,(name,mode,maxres)
        L=(1-q)*d.sum(); V=initial+new; T=V-L; R=s['total_rejection_fraction']*totalwash
        check(V+R-(T+R+L))
        c=result['counts']
        for tag,v in [('initial_assets',initial),('replacement_virgin',new),('terminal_recovered_assets',T),
                      ('wash_inputs',totalwash),('virgin_assets',V),('closed_loop_rejections',R),
                      ('empty_dirty_crate_km',DD),('empty_clean_crate_km',DS),('ready_crate_weeks',QW),('dirty_crate_weeks',DW)]:
            assert abs(c[tag]-v)<1e-5*max(1.,abs(v)),(name,mode,tag)
        cost=s['manufacturing_eur_per_crate']*(V+R)+s['wash_eur_per_crate']*totalwash+s['holding_eur_per_crate_week']*(QW+DW)+s['cost_per_truck_km']/s['declared_empty_crates_per_truck']*(DD+s['clean_leg_factor']*DS)
        assert abs(cost-result['cost_eur'])<1e-5*max(1.,cost)
        # Minimal opening stock lower bound: unreturned first tau+1 demand waves.
        startup=d[:min(tau+1,H)].sum()
        assert initial+sum(get('new',tt,slot) for tt in range(min(tau+1,H)) for slot in range(stocklen))>=startup-1e-5
        if q==1:
            rolling=np.array([d[max(0,t-tau):t+1].sum(axis=0) for t in range(H)])
            lower=float(rolling.max(axis=0).sum() if mode.startswith('dedicated') else rolling.sum(axis=1).max())
            assert initial>=lower-1e-5,(name,mode,initial,lower)
            if mode in ['dedicated_static','dedicated_adaptive']:
                for j in range(J):
                    slots=[slot for slot,(jj,i) in enumerate(result['asset_slots']) if jj==j]
                    assert sum(get('initial',slot) for slot in slots)>=rolling[:,j].max()-1e-5
                    for t in range(H):
                        # Owner stock includes ready, dirty, washing-in-process
                        # and not-yet-returned service, never another customer's.
                        owned=sum(get('ready',t,slot) for slot in slots)+get('dirty',t,j)+sum(get('wash',t,j,i) for i in range(I))+d[max(0,t-tau+1):t+1,j].sum()
                        assert abs(owned-sum(get('initial',slot) for slot in slots))<1e-5
        # Source-recipe feed feasibility; not an empirical full polymer MFA.
        secondary_available=R*1.49*.93
        secondary_needed=R*.91/.994
        assert secondary_available>=secondary_needed
        # Direct counter-based impact calculation, rather than the generator's
        # schedule-coefficient dot product. Source rounding remains correlated.
        v=np.array([V/100,T/100,R/.55,R/.55,totalwash/100])
        g=np.array([(s['transport_burden_factor']*(DD+s['clean_leg_factor']*DS)/139-totalwash)/100,
                    (s['electricity_burden_factor']-1)*totalwash/100,
                    (s['heating_burden_factor']-1)*totalwash/100])
        ratios=[]
        for k,(mat,rhs,u,proc) in polytopes.items():
            for upper in [False,True]:
                opt=linprog(-v if upper else v,A_ub=mat,b_ub=rhs,bounds=[(None,None)]*5,method='highs')
                assert opt.success
                value=float(v@opt.x)
                for tag,weight in zip(['return_transport','electricity','water_heating'],g):
                    options=[weight*z for z in proc[tag]]
                    value+=max(options) if upper else min(options)
                den=d.sum()/100*(1+(-1 if upper else 1)*comparative[k]/(2*u))
                ratio=value/den
                reported=result['impact_ratios'][k]['upper' if upper else 'lower']
                assert abs(ratio-reported)<1e-7*max(1,abs(ratio)),(name,mode,k,ratio,reported)
                ratios.append(ratio)
        assert result['relative_primal_dual_gap']<1e-8 and result['maximum_dual_stationarity_residual']<1e-7
        allchecks.append(dict(scenario=name,mode=mode,flow_and_asset_balance=True,
                              cost_recalculation=True,terminal_closure=True,
                              source_recipe_secondary_available=True,impact_envelope_recalculation=True,
                              impact_endpoint_checks=len(ratios),maximum_flow_residual=maxres))
    calibration=rows(root/'source_calendar_calibration.csv')
    for r in calibration:
        source=next(g for g in grid if g['source_scenario']=='W1_RE1' and g['indicator']==r['indicator'] and int(g['reported_rotations'])==int(r['services']))
        assert float(r['reported'])==float(source['reusable_cumulative_impact_100_crates'])
        assert abs(float(r['reconstructed'])-float(r['reported']))<=float(r['printed_half_width'])+1e-7*max(abs(float(r['reported'])),1e-20)
    out=dict(scope='independent source-calendar and event/cost/asset verification; not empirical validation or complete LCA/LCC',
             source_calendar_checks=len(calibration),successful_network_schedules=len(allchecks),checks=allchecks,
             all_passed=True,remaining='physical/source model completion and empirically representative operation are not established by arithmetic checks. HiGHS dual diagnostics are checked from reported output, not independently solved LP optimality.')
    (root/'independent_calendar_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__': main()
