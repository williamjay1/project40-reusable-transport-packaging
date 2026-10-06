"""Independent discounted-cash/event and same-draw joint-bound verifier.

No import of either financial generator or CalendarLP. Source polytope is
independently assembled; physical/source checks are also run by the separate
verify_calendar_pool_v23.py on every newly generated schedule.
"""
import argparse
import csv
import json
import math
from collections import defaultdict
from fractions import Fraction
from pathlib import Path
import numpy as np
from scipy.optimize import linprog

def rows(p):
    with p.open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))

def near(a,b):
    assert abs(a-b)<1e-8*max(1,abs(a),abs(b)),(a,b)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--results-dir',type=Path,required=True)
    ap.add_argument('--input-dir',type=Path,required=True);ap.add_argument('--stage-facts',type=Path,required=True)
    ap.add_argument('--electricity-raw',type=Path,required=True);a=ap.parse_args(); root=a.results_dir
    results=json.loads((root/'calendar_network_results.json').read_text())['outputs']
    flows=rows(root/'calendar_positive_flows.csv'); ledger=rows(root/'discounted_cashflow_ledger.csv')
    gm=json.loads((root/'geometry_provenance.json').read_text());georows=rows(root/'geometry_numeric_facts.csv')
    dist={(j,i):float(next(r['distance_km'] for r in georows if r['origin']==origin and r['destination']==dest))
        for j,origin in enumerate(gm['sorters']) for i,dest in enumerate(gm['washers'])}
    scalar_checks=0; computed={};cash_checks=0;max_cash_error=0
    for result in results:
        name=result['scenario']['name']; mode=result['mode']; s=result['scenario'];rate=result['annual_discount_rate'];v=result['terminal_net_eur_per_asset']
        factors={t:math.exp(-math.log1p(rate)*t/52) for t in range(s['weeks']+s['return_delay_weeks']+2)}
        expected=defaultdict(float); washes=0.;initial=0.;new=0.;Wpv=0.
        subset=[r for r in flows if r['scenario']==name and r['mode']==mode]
        assert len(subset)==len({(r['variable'],r['coordinates']) for r in subset})
        for r in subset:
            tag=r['variable']; c=tuple(map(int,r['coordinates'].split(':')));x=float(r['value']);assert x>=0
            if tag=='initial':expected['virgin_procurement',0]+=x*s['manufacturing_eur_per_crate'];initial+=x;continue
            t=c[0];week=t+1
            if tag=='new':expected['virgin_procurement',week]+=x*s['manufacturing_eur_per_crate'];new+=x
            elif tag in ['ready','dirty']:expected['holding',week]+=x*s['holding_eur_per_crate_week']
            elif tag=='wash':
                expected['washing',week]+=x*s['wash_eur_per_crate']
                expected['reject_recipe_procurement',week+1]+=x*s['total_rejection_fraction']*s['manufacturing_eur_per_crate']
                expected['empty_transport',week]+=x*dist[c[1],c[2]]*s['cost_per_truck_km']/s['declared_empty_crates_per_truck'];washes+=x;Wpv+=x*factors[week]
            elif tag=='dispatch':expected['empty_transport',week]+=x*dist[c[1],c[2]]*s['clean_leg_factor']*s['cost_per_truck_km']/s['declared_empty_crates_per_truck']
            else:raise AssertionError(tag)
        N=float(np.array(result['calendar_demands']).sum());T=initial+new-(1-s['collection_fraction'])*N;end=s['weeks']+s['return_delay_weeks']+1
        expected['terminal_net_value',end]=-v*T
        actual=defaultdict(float); actualPV=defaultdict(float)
        for r in ledger:
            if r['scenario']!=name or r['mode']!=mode:continue
            key=(r['component'],int(r['week']));actual[key]+=float(r['cost_eur']);actualPV[key]+=float(r['present_value_eur'])
        assert set(actual)==set(expected),(name,mode,set(actual)^set(expected))
        for key,cost in expected.items():
            near(cost,actual[key]);near(cost*factors[key[1]],actualPV[key]);cash_checks+=2
            max_cash_error=max(max_cash_error,abs(cost*factors[key[1]]-actualPV[key]))
        comps={k:sum(c*factors[t] for (kk,t),c in expected.items() if kk==k) for k in result['present_components_eur']}
        npv=sum(comps.values()); Npv=sum(sum(d)*factors[t+1] for t,d in enumerate(result['calendar_demands']))
        near(npv,result['included_present_cost_eur']);near(npv,result['PV_objective_eur']);near(Wpv,result['discounted_wash_equivalents']);near(Npv,result['discounted_service_equivalents']);scalar_checks+=4
        if 'present_cost_dual_objective_eur' in result:near(npv,result['present_cost_dual_objective_eur']);scalar_checks+=1
        for k,value in comps.items():near(value,result['present_components_eur'][k]);scalar_checks+=1
        if rate==0 and v==0:near(npv,result['cost_eur']);scalar_checks+=1
        computed[name,mode]=dict(NPV=npv,Wpv=Wpv,Npv=Npv,transport=comps['empty_transport'],terminal=T*factors[end],result=result)
    pairchecks=0
    for r in rows(root/'financial_completion_comparisons.csv'):
        D=computed[r['case'],'dedicated_adaptive'];P=computed[r['case'],'pooled'];gap=D['NPV']-P['NPV']
        near(gap,float(r['additional_net_pooling_present_cost_allowance_eur']));near(gap/P['Npv'],float(r['allowance_eur_per_discounted_service']));pairchecks+=2
        if r['fixed_policy_terminal_value_zero_gap_eur_per_asset']:
            critical=float(r['terminal_net_eur_per_asset'])+gap/(D['terminal']-P['terminal']);near(critical,float(r['fixed_policy_terminal_value_zero_gap_eur_per_asset']));pairchecks+=1
    # Price provenance/cutoff and public tax/band selection, not market fitting.
    raw=json.loads(a.electricity_raw.read_text());assert raw['updated'][:10]<='2026-08-29'
    dims=raw['id']; codes=[list(raw['dimension'][d]['category']['index']) for d in dims];prices=[]
    for index,value in raw['value'].items():
        if value is None:continue
        coords=np.unravel_index(int(index),raw['size']);cell={d:codes[k][coords[k]] for k,d in enumerate(dims)}
        if cell['tax']=='X_VAT' and cell['currency']=='EUR' and cell['unit']=='KWH' and cell['geo']=='EU27_2020' and not cell['nrg_cons'].startswith('TOT') and cell['time'] in ['2024-S1','2024-S2','2025-S1','2025-S2']:prices.append(float(value))
    assert len(prices)==28
    audit=json.loads((root/'financial_input_audit.json').read_text());near(min(prices),audit['ranges']['nrg_pc_205']['lower']);near(max(prices),audit['ranges']['nrg_pc_205']['upper']);assert 'nrg_pc_203' not in audit['ranges']
    # Assemble source constraints independently, retaining signed endpoints.
    grid=rows(a.input_dir/'tua_reported_grid.csv');stages=json.loads(a.stage_facts.read_text())['facts'];process=rows(a.input_dir/'tua_process_contributions.csv')
    fixed=dict(zip(['Climate change','Photochemical ozone formation','Acidification','Terrestrial eutrophication','Freshwater ecotoxicity','Cumulative energy demand','Water depletion'],[.01,.01,.01,.01,1,1,.01]))
    comparative=dict(zip(['Climate change','Ozone depletion','Human toxicity, non-cancer effects','Human toxicity, cancer effects','Particulate matter','Photochemical ozone formation','Acidification','Terrestrial eutrophication','Freshwater eutrophication','Marine eutrophication','Freshwater ecotoxicity','Mineral, fossil & renewable resource depletion','Cumulative energy demand','Water depletion'],[.01,1e-7,1e-7,1e-7,1e-4,1e-3,1e-3,1e-3,1e-4,1e-4,.01,1e-5,1.,.01]))
    polys={}
    for stage in stages:
        k=stage['indicator'];gs=[g for g in grid if g['indicator']==k and g['source_scenario']=='W1_RE1'];u=float(gs[0]['single_use_impact_per_service_100_crates']);A=[];b=[]
        def constrain(vec,value,half):
            A.extend([vec,[-z for z in vec]]);b.extend([(value+half)/u,(-value+half)/u])
        for g in gs:
            n=int(g['reported_rotations']);y=float(g['reusable_cumulative_impact_100_crates']);half=fixed.get(k,10**(math.floor(math.log10(y))-2))/2
            constrain([1,1,n-1,n-1,n-1],y,half)
        for n in [20,125]:
            for tag,vec in [('P',[1,0,n-1,0,0]),('E',[0,1,0,n-1,0]),('W',[0,0,0,0,n-1])]:
                c=stage['cells'][f'{tag}{n}'];constrain(vec,c['value'],(c['upper']-c['lower'])/2)
        proc={}
        for r in process:
            if r['indicator']!=k:continue
            m,_,e=r['displayed_value'].upper().partition('E'); exponent=int(e) if e else 0
            h=Fraction(1,2)*Fraction(10)**(exponent-len(m.partition('.')[2]));z=Fraction(m)*Fraction(10)**exponent
            proc[r['process']]=(float(z-h)/u,float(z+h)/u)
        polys[k]=(np.array(A),np.array(b),1-comparative[k]/(2*u),proc)
    rng=np.random.default_rng(400023);E=rng.triangular(1,1,2,50000);H=rng.triangular(1,1,2,50000); envcache={};jointchecks=0;strictchecks=0
    for r in rows(root/'joint_environment_financial_budgets.csv'):
        name=r['case'];gamma=float(r['travel_factor']);rho=float(r.get('missing_fate_virgin_equivalent_rho',0));D=computed[name,'dedicated_adaptive'];P=computed[name,'pooled'];key=(name,gamma,rho)
        if key not in envcache:
            c=P['result']['counts'];s=P['result']['scenario'];N=c['delivered_services'];W=c['wash_inputs']/100;weighted=np.array([(c['virgin_assets']+rho*c['missing_assets'])/100,c['terminal_recovered_assets']/100,c['closed_loop_rejections']/.55,c['closed_loop_rejections']/.55,W]); travel=(gamma*(c['empty_dirty_crate_km']+s['clean_leg_factor']*c['empty_clean_crate_km'])/139-c['wash_inputs'])/100
            impacts=[]
            for k,(A,b,lowU,proc) in polys.items():
                opt=linprog(-weighted,A_ub=A,b_ub=b,bounds=[(None,None)]*5,method='highs');assert opt.success
                transport=max(travel*z for z in proc['return_transport']);electric=np.maximum((E-1)*W*proc['electricity'][0],(E-1)*W*proc['electricity'][1]);heat=np.maximum((H-1)*W*proc['water_heating'][0],(H-1)*W*proc['water_heating'][1]);impacts.append((weighted@opt.x+transport+electric+heat)/(N/100*lowU))
            envcache[key]=np.max(impacts,axis=0)<1
        mask=envcache[key]; pe=float(r['electricity_price_context_eur_per_kwh']);ph=float(r['declared_effective_heat_price_eur_per_kwh']);assert pe in [min(prices),max(prices)]
        extra=.0248*pe*(E-1)+(.322/3.6)*ph*(H-1)
        savings=D['NPV']-P['NPV']+(gamma-1)*(D['transport']-P['transport'])-extra*(P['Wpv']-D['Wpv'])
        savings=np.where(np.abs(savings)<1e-6,0.,savings);required=47500
        near(mask.mean(),float(r['environmental_prior_coverage']));near(np.mean(mask&(savings>0)),float(r['joint_env_and_pool_cost_coverage_at_zero_extra_cost']));jointchecks+=2
        if mask.sum()<required:
            assert r['joint95_net_pooling_present_cost_supremum_eur']=='' and r['joint95_required_singleuse_cost_eur_per_discounted_service']=='';jointchecks+=2;continue
        # Descending order is independent of the producer's ascending index.
        budget=float(np.sort(savings[mask])[::-1][required-1]);near(budget,float(r['joint95_net_pooling_present_cost_supremum_eur']));near(budget/P['Npv'],float(r['joint95_budget_eur_per_discounted_service']));jointchecks+=2
        epsilon=max(1e-8,abs(budget)*1e-9)
        assert np.mean(mask&(savings>budget-epsilon))>=.95
        assert np.mean(mask&(savings>budget+epsilon))<=.95;strictchecks+=2
        pool_cost=P['NPV']+(gamma-1)*P['transport']+extra*P['Wpv']; SU=float(np.sort(pool_cost[mask])[required-1])/P['Npv'];near(SU,float(r['joint95_required_singleuse_cost_eur_per_discounted_service']));jointchecks+=1
    report=dict(all_passed=True,financial_schedules=len(results),cash_component_week_checks=cash_checks,financial_scalar_checks=scalar_checks,
        paired_budget_checks=pairchecks,industrial_price_cells_checked=len(prices),joint_scalar_or_missing_checks=jointchecks,strict_boundary_checks=strictchecks,
        independent_environmental_sample_vectors=len(envcache),maximum_absolute_cash_reaggregation_error_eur=max_cash_error,
        scope='Independent cashflow timing, PV/terminal accounting, source-price cutoff/tax selection and shared-draw order statistics. Physical/source endpoints checked separately. Not observed operator LCC or empirical confidence.',
        missing='Unknown K, terminal prices, effective heat prices and discount rates remain explicit coordinates; no acceptance claim')
    (root/'independent_financial_verification.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
