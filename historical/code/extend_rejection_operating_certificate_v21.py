"""Source-conditional extra POST-WASH rejection and operating-burden screen.

This is a response of published stage vectors, not a new LCA, measured failure
rate, collection-loss model, or empirical EU threshold. All added rejected crates
follow the source's closed-loop recycling/replacement inventories. The original
reconditioning vector is held fixed: added failures occur after that process.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import re
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from couple_pool_budget_and_cohort_v20 import (
    vertex_budget_coefficients, read_csv, COMPARE, FIXED)

INDICATORS=list(COMPARE)
R0=F(55,10000)

def interval(token):
    text=token.upper()
    mantissa,_,exponent=text.partition('E')
    exponent=int(exponent) if exponent else 0
    decimals=len(mantissa.partition('.')[2])
    value=F(mantissa)*F(10)**exponent
    hw=F(1,2)*F(10)**(exponent-decimals)
    return value,hw

def extract_stage_rows(path):
    before=hashlib.sha256(path.read_bytes()).hexdigest()
    block=path.read_text(encoding='utf-8').split('Table S12.',1)[1].split('Table S13.',1)[0]
    found=[]
    for line_number,line in enumerate(block.splitlines(),1):
        pairs=re.findall(r'(-?\d+(?:\.\d+)?(?:E[+-]?\d+)?)\s+(-?\d+)%',line)
        if len(pairs)==6:
            found.append((line_number,line.strip(),[x[0] for x in pairs]))
    assert len(found)==14,len(found)
    result=[]
    for k,(line_number,line,tokens) in zip(INDICATORS,found):
        cells={}
        for field,token in zip(['P20','W20','E20','P125','W125','E125'],tokens):
            value,hw=interval(token)
            cells[field]={'printed':token,'value':float(value),'lower':float(value-hw),'upper':float(value+hw)}
        # End-of-life contains the source recycling accounting. These are NET
        # stage increments, not virgin polymer or generic asset F coefficients.
        scale=100*R0*(125-20)
        lower=(F(str(cells['P125']['lower']))-F(str(cells['P20']['upper']))+
               F(str(cells['E125']['lower']))-F(str(cells['E20']['upper'])))/scale
        upper=(F(str(cells['P125']['upper']))-F(str(cells['P20']['lower']))+
               F(str(cells['E125']['upper']))-F(str(cells['E20']['lower'])))/scale
        assert 0<lower<=upper,(k,lower,upper)
        result.append({'indicator':k,'source_table':'S12','source_doi':'10.3390/resources8020110',
                       'source_line_in_table_block':line_number,'cells':cells,
                       'net_replacement_increment_lower':float(lower),
                       'net_replacement_increment_upper':float(upper),
                       'increment_unit':'source indicator unit per additional rejected/replaced crate',
                       'interpretation':'source age-invariant stage response; all added rejects recycled, not external non-return'})
    assert before==hashlib.sha256(path.read_bytes()).hexdigest()
    return result,before

def feasibility_check(stage,grid,contributions):
    """Independent five-coordinate decomposition within all printed intervals.

    p0,e0 are initial production/terminal terms divided by U; p,e are the
    per-transition production and recycling slopes divided by U; w is the
    per-transition reconditioning vector. e0 may be negative (source credit).
    """
    out=[]
    for row in stage:
        k=row['indicator'];cells=row['cells']
        rs=[v for v in grid if v['source_scenario']=='W1_RE1' and v['indicator']==k]
        u=float(rs[0]['single_use_impact_per_service_100_crates'])
        A=[];rhs=[]
        def add(vector,lower,upper):
            A.extend([vector,[-v for v in vector]])
            rhs.extend([upper/u,-lower/u])
        for v in rs:
            n=int(v['reported_rotations']);y=float(v['reusable_cumulative_impact_100_crates'])
            hw=FIXED[k]/2 if k in FIXED else 10**(np.floor(np.log10(y))-2)/2
            add([1,1,n-1,n-1,n-1],y-hw,y+hw)
        for n in [20,125]:
            for label,vec in [('P',[1,0,n-1,0,0]),('E',[0,1,0,n-1,0]),('W',[0,0,0,0,n-1])]:
                cell=cells[f'{label}{n}'];add(vec,cell['lower'],cell['upper'])
        sol=linprog(np.zeros(5),A_ub=np.array(A),b_ub=np.array(rhs),bounds=[(None,None)]*5,method='highs',
                    options={'primal_feasibility_tolerance':1e-10,'dual_feasibility_tolerance':1e-10})
        if not sol.success:
            raise AssertionError(f'Printed stage-total decomposition infeasible: {k}: {sol.message}')
        # Check tolerances relative to actual printed half-widths, not raw units.
        residual=np.array(A)@sol.x-np.array(rhs)
        out.append({'indicator':k,'feasible':True,'max_dimensionless_residual':float(residual.max())})
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input-dir',type=Path,required=True)
    source=ap.add_mutually_exclusive_group(required=True)
    source.add_argument('--source-layout',type=Path)
    source.add_argument('--stage-facts',type=Path,
                        help='Attributed S12 numeric facts; avoids redistributing full source text')
    ap.add_argument('--output-dir',type=Path,required=True)
    args=ap.parse_args()
    grid=read_csv(args.input_dir/'tua_reported_grid.csv')
    contrib=read_csv(args.input_dir/'tua_process_contributions.csv')
    if args.source_layout:
        stage,sha=extract_stage_rows(args.source_layout)
    else:
        attributed=json.loads(args.stage_facts.read_text(encoding='utf-8'))
        stage=attributed['facts'];sha=attributed['source_layout_sha256']
        assert len(stage)==14
    coefficients=vertex_budget_coefficients(grid,contrib)
    g={row['indicator']:row['net_replacement_increment_upper'] for row in stage}
    validation=feasibility_check(stage,grid,contrib)
    stage_lower={k:0. for k in coefficients}
    for row in contrib:
        value,half=interval(row['displayed_value'])
        stage_lower[row['indicator']]+=float(value-half)
    rows=[];fixed_target_rows=[]
    for B in range(3,126):
        margins={k:B*p['u_lower']-max(p['u']*(v[0]+B*v[1]) for v in p['vertices']) for k,p in coefficients.items()}
        for multiplier in [1.,2.,3.]:
            environmental={k:float(R0)+(margins[k]/(B-1)-(multiplier-1)*p['joint_stage_upper'])/(100*g[k])
                           for k,p in coefficients.items()}
            limiting=min(environmental,key=environmental.get)
            inner_pass=all(margins[k]>(B-1)*(multiplier-1)*p['joint_stage_upper']
                           for k,p in coefficients.items())
            favourable_fail=any(
                min(p['u']*(v[0]+B*v[1]) for v in p['vertices'])+
                (B-1)*(multiplier-1)*stage_lower[k] >= B*(p['u']+COMPARE[k]/2)
                for k,p in coefficients.items())
            env_outcome=('certified_advantage' if inner_pass else
                         'certified_failure' if favourable_fail else 'unresolved')
            for alpha in [F(19,20),F(99,100)]:
                cohort=float((1-alpha)/(B-1))
                ceiling=min(environmental[limiting],cohort)
                common={'service_horizon':B,'joint_process_multiplier':multiplier,'alpha':float(alpha),
                             'environmental_rejection_ceiling':environmental[limiting],
                             'environmental_limiting_indicator':limiting,
                             'environmental_outcome_at_source_rejection':env_outcome,
                             'source_rejection':float(R0)}
                rows.append({**common,'cohort_target_services':B,'target_policy':'K=n',
                             'worst_allocation_cohort_rejection_ceiling':cohort,
                             'joint_rejection_ceiling':ceiling,
                             'nonempty_above_source_rejection':bool(ceiling>float(R0)),
                             'limiting_constraint':'cohort' if cohort<environmental[limiting] else limiting})
                fixed_cohort=float((1-alpha)/2)
                fixed_ceiling=min(environmental[limiting],fixed_cohort)
                fixed_target_rows.append({**common,'cohort_target_services':3,'target_policy':'K=3 fixed',
                             'worst_allocation_cohort_rejection_ceiling':fixed_cohort,
                             'joint_rejection_ceiling':fixed_ceiling,
                             'nonempty_above_source_rejection':bool(fixed_ceiling>float(R0)),
                             'limiting_constraint':'cohort' if fixed_cohort<environmental[limiting] else limiting})
    args.output_dir.mkdir(parents=True,exist_ok=True)
    (args.output_dir/'replacement_stage_facts_v21.json').write_text(json.dumps(
        {'source_layout_sha256':sha,'source_table':'S12','source_doi':'10.3390/resources8020110','facts':stage},indent=2)+'\n',encoding='utf-8')
    with (args.output_dir/'joint_rejection_operating_frontier_v21.csv').open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    with (args.output_dir/'fixed_three_service_target_frontier_v22.csv').open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(fixed_target_rows[0]));writer.writeheader();writer.writerows(fixed_target_rows)
    summary={'source_rejection_exact':str(R0),'source_stage_fact_count':14*6,'feasible_stage_total_checks':validation,
             'evidence_status':'controlled additional post-wash reject sensitivity, not measured return failure or a new LCA',
             'scope':['W1-RE1 only','r>=0.0055; additional post-process rejection','closed-loop source replacement/recycling retained',
                      'unchanged washing/transport per 100 inputs','no external loss or credits for unreturned material',
                      'synchronous assured next service for cohort condition','rounding envelopes, not parameter probability'],
             'selected_rows':[r for r in rows if r['service_horizon'] in [3,5,8,10,20]],
             'row_count':len(rows),'fixed_target_row_count':len(fixed_target_rows),
             'cohort_target_definition':'K distinct from environmental horizon n; diagonal and fixed-K=3 both retained',
             'fixed_K3_source_bounds':[0.989,0.9945],
             'three_outcome_rule':'inner sufficient advantage; favourable endpoint failure; otherwise unresolved',
             'fixed_target_selected_rows':[r for r in fixed_target_rows if r['service_horizon'] in [3,5,8,10,20]]}
    (args.output_dir/'joint_rejection_summary_v21.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
