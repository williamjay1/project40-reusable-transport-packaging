"""Deterministically integrate the declared independent triangular energy prior.

This keeps every v29 policy/input fixed. Exact polynomial polygon moments are
numerical integration under a design prior, never empirical EU probability.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from calendar_pool_network_v23 import source_polytopes,dump_csv

def clip_polygon(poly,a,b,c):
    """Strict half-plane differs only on measure-zero edges, not constant ties."""
    if not poly:return []
    if abs(a)+abs(b)<1e-15:return poly if c>0 else []
    result=[]
    for prev,cur in zip(poly[-1:]+poly[:-1],poly):
        f0=a*prev[0]+b*prev[1]-c;f1=a*cur[0]+b*cur[1]-c
        in0=f0<=0;in1=f1<=0
        if in0!=in1:
            t=f0/(f0-f1);result.append((prev[0]+t*(cur[0]-prev[0]),prev[1]+t*(cur[1]-prev[1])))
        if in1:result.append(cur)
    return result

def integrate(planes):
    poly=[(0.,0.),(1.,0.),(1.,1.),(0.,1.)]
    for a,b,c in planes:
        poly=clip_polygon(poly,a,b,c)
        if len(poly)<3:return 0.
    A=Mx=My=Mxy=0.
    for (x,y),(X,Y) in zip(poly,poly[1:]+poly[:1]):
        cross=x*Y-X*y;A+=cross/2;Mx+=(x+X)*cross/6;My+=(y+Y)*cross/6
        Mxy+=(2*x*y+x*Y+X*y+2*X*Y)*cross/24
    value=4*(A-Mx-My+Mxy)
    assert -1e-10<=value<=1+1e-10,value
    return min(1.,max(0.,float(value)))

def inverse_budget(planes,g,e,h,target=.95):
    probability=integrate(planes)
    if probability<target-1e-12:return None
    if abs(e)+abs(h)<1e-15:return g
    corners=[g-e*x-h*y for x in [0.,1.] for y in [0.,1.]]
    lo=min(corners)-1e-9;hi=max(corners)+1e-9
    for _ in range(64):
        mid=(lo+hi)/2;p=integrate([*planes,(e,h,g-mid)])
        if p>=target:lo=mid
        else:hi=mid
    return (lo+hi)/2

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--results-dir',type=Path,required=True);ap.add_argument('--input-dir',type=Path,required=True)
    ap.add_argument('--stage-facts',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True);a=ap.parse_args();out=a.output_dir;out.mkdir(parents=True,exist_ok=True)
    controls=dict(square=integrate([]),empty=integrate([(0.,0.,-1.)]),constant_strict_tie=integrate([(0.,0.,0.)]),half_axis=integrate([(1.,0.,.5)]),diagonal=integrate([(1.,1.,1.)]))
    for name,expected in [('square',1.),('empty',0.),('constant_strict_tie',0.),('half_axis',.75),('diagonal',5/6)]:assert abs(controls[name]-expected)<1e-12,(name,controls[name],expected)
    source=source_polytopes(a.input_dir,a.stage_facts);meta=json.loads((a.results_dir/'calendar_network_results.json').read_text());res=meta['outputs'];idx={(r['financial_case_id'],r['mode']):r for r in res}
    with (a.results_dir/'joint_environment_financial_budgets.csv').open(newline='',encoding='utf-8') as f:old=list(csv.DictReader(f))
    coefficients={};rows=[];budgetdiff=[];probdiff=[];classification=[]
    for r in old:
        case=r['case'];gamma=float(r['travel_factor']);rho=float(r['missing_fate_virgin_equivalent_rho']);key=(case,gamma,rho);P=idx[case,'pooled'];D=idx[case,'dedicated_adaptive'];Npv=P['discounted_service_equivalents']
        if key not in coefficients:
            c=P['counts'];s=P['scenario'];N=c['delivered_services'];W=c['wash_inputs']/100
            vector=np.array([(c['virgin_assets']+rho*c['missing_assets'])/100,c['terminal_recovered_assets']/100,c['closed_loop_rejections']/.55,c['closed_loop_rejections']/.55,W])
            travel=(gamma*(c['empty_dirty_crate_km']+s['clean_leg_factor']*c['empty_clean_crate_km'])/139-c['wash_inputs'])/100
            records=[]
            for k,p in source.items():
                opt=linprog(-vector,A_ub=p['A'],b_ub=p['b'],bounds=[(None,None)]*5,method='highs');assert opt.success
                denom=N/100*p['u_lower_ratio'];base=(float(vector@opt.x)+max(travel*p['process']['return_transport'][0],travel*p['process']['return_transport'][1]))/denom
                E=W*p['process']['electricity'][1]/denom;H=W*p['process']['water_heating'][1]/denom
                records.append(dict(indicator=k,a=E,b=H,c=1-base,base_ratio=base))
            coefficients[key]=records
        planes=[(z['a'],z['b'],z['c']) for z in coefficients[key]]
        pe=float(r['electricity_price_context_eur_per_kwh']);ph=float(r['declared_effective_heat_price_eur_per_kwh']);Wgap=P['discounted_wash_equivalents']-D['discounted_wash_equivalents']
        g=(D['included_present_cost_eur']-P['included_present_cost_eur']+(gamma-1)*(D['present_components_eur']['empty_transport']-P['present_components_eur']['empty_transport']))/Npv
        e=.0248*pe*Wgap/Npv;h=(.322/3.6)*ph*Wgap/Npv
        # Flat controls can have roundoff-only coefficients. Apply the existing
        # campaign-neutrality rule uniformly if ALL possible savings are tiny.
        if max(abs(g-e*x-h*y)*Npv for x in [0.,1.] for y in [0.,1.])<1e-6:g=e=h=0.
        env=integrate(planes);joint=integrate([*planes,(e,h,g)]);budget=inverse_budget(planes,g,e,h)
        prior_budget=float(r['joint95_budget_eur_per_discounted_service']) if r['joint95_budget_eur_per_discounted_service'] else None
        pd=env-float(r['environmental_prior_coverage']);jd=joint-float(r['joint_env_and_pool_cost_coverage_at_zero_extra_cost']);probdiff.extend([abs(pd),abs(jd)])
        bd=None if budget is None or prior_budget is None else budget-prior_budget
        if bd is not None:budgetdiff.append(abs(bd))
        changed=(budget is None)!=(prior_budget is None) or (budget is not None and prior_budget is not None and (budget>1e-9)!=(prior_budget>1e-9))
        if changed:classification.append(dict(case=case,gamma=gamma,rho=rho,pe=pe,ph=ph,MC_budget=prior_budget,integrated_budget=budget))
        row={k:r[k] for k in ['case','parent_scenario','annual_discount_rate','terminal_net_eur_per_asset','travel_factor','missing_fate_virgin_equivalent_rho','electricity_price_context_eur_per_kwh','declared_effective_heat_price_eur_per_kwh']}
        row.update(integrated_environmental_prior_coverage=env,integrated_joint_prior_coverage_at_zero_extra_cost=joint,integrated_joint95_budget_eur_per_discounted_service=budget,MC_environmental_coverage=float(r['environmental_prior_coverage']),MC_joint_coverage=float(r['joint_env_and_pool_cost_coverage_at_zero_extra_cost']),MC_joint95_budget_eur_per_discounted_service=prior_budget,integrated_minus_MC_environmental_probability=pd,integrated_minus_MC_joint_probability=jd,integrated_minus_MC_budget_eur_per_discounted_service=bd,cost_base_allowance_eur_per_discounted_service=g,cost_energy_slope_e=e,cost_energy_slope_h=h,classification_changed=changed)
        rows.append(row)
    assert len(rows)==432
    dump_csv(out/'integrated_joint_decisions.csv',rows)
    dump_csv(out/'environmental_halfplanes.csv',[dict(case=key[0],travel_factor=key[1],missing_fate_virgin_equivalent_rho=key[2],**r) for key,records in coefficients.items() for r in records])
    report=dict(controls=controls,decisions=len(rows),unique_environmental_vectors=len(coefficients),maximum_absolute_probability_difference=max(probdiff),maximum_absolute_budget_difference_eur_per_discounted_service=max(budgetdiff),classification_changes=classification,independent_check_pending=True,scope='Deterministic integration of the same declared independent triangular prior; no new empirical probabilities or policy reoptimization',source_network_sha256=hashlib.sha256((a.results_dir/'calendar_network_results.json').read_bytes()).hexdigest(),method='Degree-two density 4(1-u)(1-h) over joint half-plane intersection; Green/shoelace polynomial moments and monotone budget bisection')
    (out/'integration_summary.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
