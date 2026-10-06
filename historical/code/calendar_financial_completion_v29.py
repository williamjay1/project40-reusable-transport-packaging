"""Discounted, terminal-aware cost completion and same-draw joint budgets.

The objective does NOT fill unknown facility/labour/coordination/exit costs with
zero facts. K collects their present values; results are inverse allowances for
K. No observed single-use price, terminal value or discount rate is asserted.
Only extra energy above the bundled reference wash is priced at current public
electricity contexts / declared effective-heat coordinates.
"""
import argparse
import csv
import json
import math
from dataclasses import asdict, replace
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from calendar_pool_network_v23 import CalendarLP, Scenario, cached_geometry, source_polytopes, dump_csv, calibrate

def df(rate,week): return (1+rate)**(-week/52)

def discounted_objective(model,rate,v):
    """Operating invoices at week end; reject replacement paid one week later."""
    s=model.s; c=np.zeros(model.N); terminal=df(rate,s.weeks+s.return_delay_weeks+1)
    for key,ix in model.index.items():
        tag=key[0]
        if tag=='initial': c[ix]=s.manufacturing_eur_per_crate-v*terminal
        else:
            t=key[1]; discount=df(rate,t+1)
            if tag=='new': c[ix]=s.manufacturing_eur_per_crate*discount-v*terminal
            elif tag in ['ready','dirty']: c[ix]=s.holding_eur_per_crate_week*discount
            elif tag=='wash': c[ix]=(s.wash_eur_per_crate+s.cost_per_truck_km/s.declared_empty_crates_per_truck*model.distance[key[2],key[3]])*discount+s.total_rejection_fraction*s.manufacturing_eur_per_crate*df(rate,t+2)
            elif tag=='dispatch': c[ix]=s.cost_per_truck_km/s.declared_empty_crates_per_truck*s.clean_leg_factor*model.distance[key[2],key[3]]*discount
    return c,v*terminal*model.total_missing

def cash_ledger(model,x,rate,v):
    s=model.s; rows=[]
    for key,ix in model.index.items():
        val=float(x[ix])
        if val<=1e-7: continue
        tag=key[0]
        if tag=='initial': costs=[('virgin_procurement',0,val*s.manufacturing_eur_per_crate)]
        else:
            t=key[1]; week=t+1
            if tag=='new': costs=[('virgin_procurement',week,val*s.manufacturing_eur_per_crate)]
            elif tag in ['ready','dirty']: costs=[('holding',week,val*s.holding_eur_per_crate_week)]
            elif tag=='wash': costs=[('washing',week,val*s.wash_eur_per_crate),
                ('reject_recipe_procurement',week+1,val*s.total_rejection_fraction*s.manufacturing_eur_per_crate),
                ('empty_transport',week,val*s.cost_per_truck_km/s.declared_empty_crates_per_truck*model.distance[key[2],key[3]])]
            else: costs=[('empty_transport',week,val*s.cost_per_truck_km/s.declared_empty_crates_per_truck*s.clean_leg_factor*model.distance[key[2],key[3]])]
        for component,week,cost in costs: rows.append(dict(component=component,week=week,cost_eur=cost,present_value_eur=cost*df(rate,week)))
    counts=model.counts(x); terminal_week=s.weeks+s.return_delay_weeks+1
    rows.append(dict(component='terminal_net_value',week=terminal_week,cost_eur=-v*counts['terminal_recovered_assets'],present_value_eur=-v*counts['terminal_recovered_assets']*df(rate,terminal_week)))
    components={k:sum(r['present_value_eur'] for r in rows if r['component']==k) for k in sorted({r['component'] for r in rows})}
    Wpv=sum(x[ix]*df(rate,key[1]+1) for key,ix in model.index.items() if key[0]=='wash')
    Npv=sum(model.dem[t].sum()*df(rate,t+1) for t in range(s.weeks))
    return rows,components,float(Wpv),float(Npv)

def upper_prior(result,source,le,lh,gamma,rho=0.):
    counts=result['counts']; s=result['scenario']; N=counts['delivered_services']; W=counts['wash_inputs']/100
    travel=(counts['empty_dirty_crate_km']+s['clean_leg_factor']*counts['empty_clean_crate_km'])/(100*139)
    vec=np.array([(counts['virgin_assets']+rho*counts['missing_assets'])/100,counts['terminal_recovered_assets']/100,counts['closed_loop_rejections']/.55,counts['closed_loop_rejections']/.55,W])
    ratios=[]
    for k,p in source.items():
        opt=linprog(-vec,A_ub=p['A'],b_ub=p['b'],bounds=[(None,None)]*5,method='highs'); assert opt.success
        base=float(vec@opt.x); g=gamma*travel-W
        trans=max(g*p['process']['return_transport'][0],g*p['process']['return_transport'][1])
        energy=(le-1)*W*p['process']['electricity'][1]+(lh-1)*W*p['process']['water_heating'][1]
        ratios.append((base+trans+energy)/(N/100*p['u_lower_ratio']))
    return np.max(np.asarray(ratios),axis=0)

def joint_budget(env,saving,prob=.95):
    """Supremum K for strictly beneficial cost AND all-indicator environment.

If m environmentally passing draws exist, at least ceil(prob*n) must also
have saving>K. K lies strictly below their kth-largest saving, not below an
independently multiplied quantile. No budget exists if m is too small.
"""
    # Exact flat controls can differ by ~1e-12 EUR after reaggregation. Never
    # promote solver roundoff to a 100% probability of financial advantage.
    saving=np.where(np.abs(saving)<1e-6,0.,saving)
    mask=env<1; required=math.ceil(prob*len(env)); vals=np.sort(saving[mask]); m=len(vals)
    if m<required: return None,None,float(np.mean(mask & (saving>0)))
    budget=float(vals[m-required])
    return budget,m-required,float(np.mean(mask & (saving>0)))

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--input-dir',type=Path,required=True)
    ap.add_argument('--stage-facts',type=Path,required=True); ap.add_argument('--geometry-cache',type=Path,required=True)
    ap.add_argument('--price-dir',type=Path,required=True); ap.add_argument('--output-dir',type=Path,required=True)
    ap.add_argument('--quick',action='store_true');ap.add_argument('--append-loss-tests',action='store_true')
    a=ap.parse_args(); out=a.output_dir;out.mkdir(parents=True,exist_ok=True)
    _,_,distance=cached_geometry(a.geometry_cache,out); source=source_polytopes(a.input_dir,a.stage_facts)
    prices=json.loads((a.price_dir/'financial_input_audit.json').read_text())['ranges']['nrg_pc_205']
    (out/'financial_input_audit.json').write_bytes((a.price_dir/'financial_input_audit.json').read_bytes())
    scenarios=[Scenario('flat_tau1'),Scenario('staggered_tau1',seasonal_amplitude=.5,seasonal_phase='staggered'),Scenario('staggered_tau3',seasonal_amplitude=.5,seasonal_phase='staggered',return_delay_weeks=3),Scenario('staggered_tau1_H24',weeks=24,seasonal_amplitude=.5,seasonal_phase='staggered'),Scenario('staggered_tau1_H48',weeks=48,seasonal_amplitude=.5,seasonal_phase='staggered')]
    losses=[Scenario('staggered_tau1_q94',seasonal_amplitude=.5,seasonal_phase='staggered',collection_fraction=.94),Scenario('staggered_tau3_q99',seasonal_amplitude=.5,seasonal_phase='staggered',return_delay_weeks=3,collection_fraction=.99)]
    scenarios=losses if a.append_loss_tests else scenarios+losses
    if a.quick: scenarios=scenarios[:1]
    outputs=[]; flows=[]; cash=[]; summaries=[]; joint=[]
    if a.append_loss_tests:
        outputs=json.loads((out/'calendar_network_results.json').read_text())['outputs']
        def read_csv(path):
            with path.open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f))
        flows=read_csv(out/'calendar_positive_flows.csv');cash=read_csv(out/'discounted_cashflow_ledger.csv')
        summaries=read_csv(out/'financial_completion_comparisons.csv')
        joint=read_csv(out/'joint_environment_financial_budgets.csv')
        assert not any(r['parent_scenario'] in [s.name for s in losses] for r in outputs),'Do not duplicate appended schedules'
    rng=np.random.default_rng(400023); le=rng.triangular(1,1,2,50000); lh=rng.triangular(1,1,2,50000)
    for original in scenarios:
        for rate in [0.,.08]:
            for v in [0.,4.55]:
                case=f'{original.name}_r{round(rate*100)}_v{v:g}'
                pair={}
                for mode in ['dedicated_adaptive','pooled']:
                    s=replace(original,name=case); model=CalendarLP(s,distance,mode)
                    objective,constant=discounted_objective(model,rate,v); model.cost=objective
                    result,x=model.solve(source); assert result['success'],(case,mode,result)
                    ledger,components,Wpv,Npv=cash_ledger(model,x,rate,v)
                    included=sum(components.values()); assert abs(included-result['cost_eur']-constant)<1e-5
                    reference=(s.manufacturing_eur_per_crate*(result['counts']['virgin_assets']+result['counts']['closed_loop_rejections'])+s.wash_eur_per_crate*result['counts']['wash_inputs']+s.holding_eur_per_crate_week*(result['counts']['ready_crate_weeks']+result['counts']['dirty_crate_weeks'])+s.cost_per_truck_km/s.declared_empty_crates_per_truck*(result['counts']['empty_dirty_crate_km']+s.clean_leg_factor*result['counts']['empty_clean_crate_km']))
                    # Keep independent physical/source verifier's reference-cost field.
                    # PV objective and PV optimality diagnostics are named separately.
                    present_dual=result.pop('dual_objective_eur')+constant
                    result.update(scenario=asdict(s),mode=mode,cost_eur=reference,cost_eur_per_service=reference/model.total_services,
                        financial_case_id=case,parent_scenario=original.name,annual_discount_rate=rate,terminal_net_eur_per_asset=v,
                        included_present_cost_eur=included,PV_objective_eur=float(objective@x+constant),
                        present_cost_dual_objective_eur=present_dual,
                        present_components_eur=components,discounted_wash_equivalents=Wpv,discounted_service_equivalents=Npv,
                        objective='discounted included cashflows less common-date terminal net value; unknown K kept external',
                        environmental_constraints=False)
                    outputs.append(result); flows+=model.schedule(x,case)
                    for r in ledger: r.update(scenario=case,mode=mode)
                    cash+=ledger; pair[mode]=result
                    print(case,mode,round(included,6),round(result['maximum_upper_ratio'],6),flush=True)
                D=pair['dedicated_adaptive']; P=pair['pooled']; gap=D['included_present_cost_eur']-P['included_present_cost_eur']
                terminal_gap=(D['counts']['terminal_recovered_assets']-P['counts']['terminal_recovered_assets'])*df(rate,original.weeks+original.return_delay_weeks+1)
                # Fixed paired policies: do not extrapolate their critical v to
                # the optimal-policy frontier across uncomputed terminal values.
                critical=v+gap/terminal_gap if terminal_gap>1e-8 else None
                summaries.append(dict(case=case,parent_scenario=original.name,annual_discount_rate=rate,terminal_net_eur_per_asset=v,
                    dedicated_included_present_cost_eur=D['included_present_cost_eur'],pooled_included_present_cost_eur=P['included_present_cost_eur'],
                    additional_net_pooling_present_cost_allowance_eur=gap,allowance_eur_per_discounted_service=gap/P['discounted_service_equivalents'],
                    fixed_policy_terminal_value_zero_gap_eur_per_asset=critical,additional_discounted_washes=P['discounted_wash_equivalents']-D['discounted_wash_equivalents'],
                    dedicated_initial_assets=D['counts']['initial_assets'],pooled_initial_assets=P['counts']['initial_assets'],
                    dedicated_upper_impact=D['maximum_upper_ratio'],pooled_upper_impact=P['maximum_upper_ratio']))
                for gamma in [1.,3.]:
                    env=upper_prior(P,source,le,lh,gamma)
                    envmask=env<1; required=math.ceil(.95*len(env)); envvals_count=int(envmask.sum())
                    transport_gap=D['present_components_eur']['empty_transport']-P['present_components_eur']['empty_transport']
                    for pe in [prices['lower'],prices['upper']]:
                        for ph in [0.,.10,.20]:
                            unit_extra=(le-1)*.0248*pe+(lh-1)*(.322/3.6)*ph
                            saving=gap+(gamma-1)*transport_gap-unit_extra*(P['discounted_wash_equivalents']-D['discounted_wash_equivalents'])
                            budget,index,coverage=joint_budget(env,saving)
                            pool_cost=P['included_present_cost_eur']+(gamma-1)*P['present_components_eur']['empty_transport']+unit_extra*P['discounted_wash_equivalents']
                            SU=(float(np.sort(pool_cost[envmask])[required-1])/P['discounted_service_equivalents']) if envvals_count>=required else None
                            joint.append(dict(case=case,parent_scenario=original.name,annual_discount_rate=rate,terminal_net_eur_per_asset=v,
                                travel_factor=gamma,electricity_price_context_eur_per_kwh=pe,declared_effective_heat_price_eur_per_kwh=ph,
                                environmental_prior_coverage=float(envmask.mean()),joint_env_and_pool_cost_coverage_at_zero_extra_cost=coverage,
                                joint95_net_pooling_present_cost_supremum_eur=budget,joint95_budget_eur_per_discounted_service=None if budget is None else budget/P['discounted_service_equivalents'],
                                budget_strictly_below_supremum=True,joint95_required_singleuse_cost_eur_per_discounted_service=SU,
                                SU_missing_cost_rule='add K_P/discounted_services; no observed single-use tariff claimed',draws=50000,seed=400023,
                                cost_prior_link='same energy-input draws as pooled all-14-indicator environmental coverage; not product of marginal probabilities'))
    dump_csv(out/'calendar_positive_flows.csv',flows); dump_csv(out/'discounted_cashflow_ledger.csv',cash)
    dump_csv(out/'financial_completion_comparisons.csv',summaries);dump_csv(out/'joint_environment_financial_budgets.csv',joint)
    dump_csv(out/'source_calendar_calibration.csv',calibrate(source,a.input_dir))
    meta=dict(scope='Parametrically complete differential lifecycle-cost decision, not fully observed current operator LCC',
        timing='Initial purchase at t=0; operating invoices at week end t+1; reject-source-recipe purchase at t+2; common terminal net value at H+tau+1; 52 weeks/year',
        unknown_costs='All omitted capital/site activation, extra labour, repair, coordination, irrecoverable taxes, end treatment fees/revenues and distribution differentials enter present-value K. Count only costs not already inside reference bundled unit charges; do not double count depreciation plus purchase or interest plus discount.',
        terminal='v is a net monetary recovery value per actually recovered asset, not a LCIA credit or measured scrap quote. v=4.55 is an intentionally extreme full-new-procurement-value stress; no continuing-stock environmental allocation is inferred.',
        rates='0 and 8% annual are declared sensitivity coordinates, not market estimates. All prices/reference coefficients are dimensional cost contexts, not a claim of a consistent 2026 price base.',
        energy_cost='Charge only incremental electricity/heat input above historical bundled baseline; heat price is a declared effective delivered-heat coordinate (0,0.10,0.20), not the excluded late gas API snapshot.',
        optimization='Each of 28 full-matrix financial cases reoptimizes BOTH ownership modes with identical coefficients/capacity/arcs. Joint gamma/energy screens keep these reference-energy policies fixed; no stochastic-optimality claim. q=.94/tau1 and q=.99/tau3 loss tests were added after the existing collection/fate grid established its first passing brackets; not represented as preregistered.',
        uncertainty='50k common draws, seed400023, independent triangular(1,1,2); source print polytopes remain set bounds. Financial and environment events share draws.',
        monetary_indifference='Absolute campaign savings below 1e-6 EUR are classified as zero, uniformly for all cases/draws; zero is not a strict advantage.',
        outputs=outputs,comparisons=summaries,joint_budget_rows=len(joint),environmental_endpoint_verification_pending=True)
    (out/'calendar_network_results.json').write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(successful_reoptimized_schedules=len(outputs),financial_pairs=len(summaries),joint_budget_rows=len(joint),positive_flows=len(flows)),indent=2))

if __name__=='__main__': main()
