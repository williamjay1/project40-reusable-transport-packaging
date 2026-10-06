"""Source-preserving periodic boundary test, not an observed operating LCA.

beta is attributed recovered virgin-capital pairs per standing asset-week.
No empirical fleet life, new rejection hazard or continuing-stock free LCA is
inferred. Historical/source inputs, all scenarios and failures are retained.
"""
from __future__ import annotations
import argparse
import json
import math
import os
import time
from pathlib import Path
import numpy as np
from scipy.optimize import linprog, brentq
from scipy.sparse import coo_matrix
from calendar_pool_network_v23 import source_polytopes, cached_geometry, dump_csv
from run_calendar_phase_v33 import profile

H=12
R0=.0055
CM,CW,CH,CT=4.55,.07,.03,.5/1000
BETA_COST=1/12
BETAS=[0.,1/52,1/12]


def periodic_bounds(d,tau,q):
    rolling=sum(np.roll(d,a,axis=0) for a in range(tau))+q*np.roll(d,tau,axis=0)
    return float(np.max(rolling.sum(axis=1))),float(rolling.max(axis=0).sum())


class PeriodicLP:
    def __init__(self,d,distance,tau,q,mode):
        self.d=np.asarray(d,float); self.distance=distance; self.tau=tau; self.q=q; self.mode=mode
        self.H,self.J=self.d.shape; self.I=distance.shape[1]
        assert self.H==H and 0<q<=1 and 0<tau<H and self.J==14 and self.I==35
        self.slots=[(j,i) for j in range(self.J) for i in range(self.I)] if mode=='dedicated_adaptive' else list(range(self.I))
        self.S=len(self.slots); self.owners=self.J if mode=='dedicated_adaptive' else 1
        self.keys=[]; self.index={}
        def add(tag,shape):
            for coords in np.ndindex(shape):
                key=(tag,)+coords; self.index[key]=len(self.keys); self.keys.append(key)
        add('fleet',(self.owners,)); add('ready',(H,self.S)); add('new',(H,self.S))
        add('dispatch',(H,self.J,self.I)); add('wash',(H,self.J,self.I)); add('dirty',(H,self.J))
        self.nvar=len(self.keys); idx=self.index
        self.c=np.zeros(self.nvar)
        if not np.all(np.isfinite(distance)) or np.any(distance<0):
            raise ValueError('Distance facts must be finite and nonnegative')
        self.allowed=distance>0
        if not np.all(self.allowed.any(axis=1)):
            raise ValueError('Each customer requires a positive-distance route')
        inactive_sites=set(np.flatnonzero(~self.allowed.any(axis=0)).tolist())
        self.bounds=[(0,None)]*self.nvar
        for key,v in idx.items():
            if key[0]=='fleet': self.c[v]=CM*BETA_COST*H
            elif key[0]=='new': self.c[v]=CM
            elif key[0] in ['ready','dirty']: self.c[v]=CH
            elif key[0]=='wash': self.c[v]=CW+R0*CM+CT*distance[key[2],key[3]]
            elif key[0]=='dispatch': self.c[v]=CT*distance[key[2],key[3]]
            if key[0] in ['dispatch','wash'] and not self.allowed[key[2],key[3]]:
                self.bounds[v]=(0,0)
            if key[0] in ['ready','new']:
                site=self.slots[key[-1]][1] if mode=='dedicated_adaptive' else self.slots[key[-1]]
                if site in inactive_sites:
                    self.bounds[v]=(0,0)
        eq=[]; rhs=[]; cap=[]; caps=[]
        def equation(terms,value): eq.append(terms); rhs.append(value)
        for t in range(H):
            prev=(t-1)%H; ret=(t-tau)%H
            for j in range(self.J):
                equation([(idx['dispatch',t,j,i],1.) for i in range(self.I)],self.d[t,j])
                equation([(idx['dirty',t,j],1.),(idx['dirty',prev,j],-1.)]+
                         [(idx['wash',t,j,i],1.) for i in range(self.I)],q*self.d[ret,j])
            for slot in range(self.S):
                terms=[(idx['ready',t,slot],1.),(idx['ready',prev,slot],-1.),(idx['new',t,slot],-1.)]
                if mode=='dedicated_adaptive':
                    j,i=self.slots[slot]
                    terms.extend([(idx['dispatch',t,j,i],1.),(idx['wash',prev,j,i],-1.)])
                else:
                    terms.extend([(idx['dispatch',t,j,slot],1.) for j in range(self.J)])
                    terms.extend([(idx['wash',prev,j,slot],-1.) for j in range(self.J)])
                equation(terms,0.)
            if mode=='dedicated_adaptive':
                for j in range(self.J):
                    equation([(idx['new',t,slot],1.) for slot,(jj,i) in enumerate(self.slots) if jj==j],(1-q)*self.d[ret,j])
            else:
                equation([(idx['new',t,i],1.) for i in range(self.I)],(1-q)*self.d[ret].sum())
            for owner in range(self.owners):
                js=[owner] if mode=='dedicated_adaptive' else list(range(self.J))
                ss=[slot for slot,(j,i) in enumerate(self.slots) if j==owner] if mode=='dedicated_adaptive' else list(range(self.S))
                terms=[(idx['fleet',owner],1.)]
                terms.extend([(idx['ready',t,s],-1.) for s in ss])
                terms.extend([(idx['dirty',t,j],-1.) for j in js])
                terms.extend([(idx['wash',t,j,i],-1.) for j in js for i in range(self.I)])
                pipeline=sum(float(self.d[(t-a)%H,js].sum()) for a in range(tau))
                equation(terms,pipeline)
            for i in range(self.I):
                cap.append([(idx['wash',t,j,i],1.) for j in range(self.J)]); caps.append(2000.)
        def sparse(rows):
            rr=[]; cc=[]; vv=[]
            for r,terms in enumerate(rows):
                for col,value in terms: rr.append(r);cc.append(col);vv.append(value)
            return coo_matrix((vv,(rr,cc)),shape=(len(rows),self.nvar)).tocsr()
        self.Aeq=sparse(eq);self.beq=np.asarray(rhs)
        self.Aub=sparse(cap);self.bub=np.asarray(caps)

    def solve(self):
        sol=linprog(self.c,A_ub=self.Aub,b_ub=self.bub,A_eq=self.Aeq,b_eq=self.beq,bounds=self.bounds,method='highs')
        if not sol.success: return dict(success=False,status=sol.status,message=sol.message),None
        x=sol.x
        sums={tag:sum(float(x[v]) for key,v in self.index.items() if key[0]==tag)
              for tag in ['fleet','ready','dirty','wash','new','dispatch']}
        km={tag:sum(float(x[v])*self.distance[key[2],key[3]] for key,v in self.index.items() if key[0]==tag)
            for tag in ['wash','dispatch']}
        N=float(self.d.sum()); W=sums['wash']; B=sums['fleet']
        components=dict(stock_allocation=CM*BETA_COST*H*B,
                        missing_and_source_rejects=CM*(sums['new']+R0*W),
                        washing=CW*W,holding=CH*(sums['ready']+sums['dirty']),
                        transport=CT*(km['wash']+km['dispatch']))
        dual=float(self.beq@sol.eqlin.marginals+self.bub@sol.ineqlin.marginals)
        stat=self.c-self.Aeq.T@sol.eqlin.marginals-self.Aub.T@sol.ineqlin.marginals-sol.lower.marginals-sol.upper.marginals
        residual=max(float(np.abs(self.Aeq@x-self.beq).max()),max(0.,float((self.Aub@x-self.bub).max())))
        assert residual<1e-5 and np.max(np.abs(stat))<1e-7
        assert abs(float(self.c@x)-dual)<1e-5
        assert abs(sum(components.values())-float(self.c@x))<1e-5
        assert abs(W-self.q*N)<1e-5 and abs(sums['new']-(1-self.q)*N)<1e-5
        lowerP,lowerD=periodic_bounds(self.d,self.tau,self.q)
        lower=lowerD if self.mode=='dedicated_adaptive' else lowerP
        assert B>=lower-1e-5
        counts=dict(standing_fleet=B,ready_asset_weeks=sums['ready'],dirty_asset_weeks=sums['dirty'],wash_inputs=W,
                    replacement_virgin=sums['new'],source_rejects=R0*W,delivered_services=N,
                    empty_dirty_crate_km=km['wash'],empty_clean_crate_km=km['dispatch'])
        return dict(success=True,mode=self.mode,standing_fleet_by_owner=[float(x[self.index['fleet',j]]) for j in range(self.owners)],
                    counts=counts,periodic_fleet_lower_bound=lower,maximum_flow_residual=residual,
                    cost_eur=float(self.c@x),cost_eur_per_service=float(self.c@x/N),cost_components_eur=components,
                    dual_objective_eur=dual,maximum_stationarity_residual=float(np.max(np.abs(stat)))),x


def support(counts,p,beta,upper):
    N=counts['delivered_services']; B=counts['standing_fleet']; W=counts['wash_inputs']
    capital=beta*H*B; lost=counts['replacement_virgin']
    vec=np.array([lost+capital,capital,W,W,W])/N
    opt=linprog(-vec if upper else vec,A_ub=p['A'],b_ub=p['b'],bounds=[(None,None)]*5,method='highs')
    assert opt.success
    value=float(vec@opt.x)
    g=((counts['empty_dirty_crate_km']+counts['empty_clean_crate_km'])/139-W)/N
    lo,hi=p['process']['return_transport']
    value+=g*(hi if (g>=0)==upper else lo)
    return value/(p['u_lower_ratio'] if upper else p['u_upper_ratio'])


def envelope(counts,source,beta):
    impacts={k:dict(lower=support(counts,p,beta,False),upper=support(counts,p,beta,True)) for k,p in source.items()}
    upper=max(r['upper'] for r in impacts.values()); lower=max(r['lower'] for r in impacts.values())
    return dict(beta=beta,impact_ratios=impacts,maximum_upper_ratio=upper,
                limiting_indicator=max(impacts,key=lambda k:impacts[k]['upper']),
                outcome='certified_source_advantage' if upper<1 else 'certified_source_failure' if lower>=1 else 'unresolved')


def threshold(counts,source):
    f=lambda beta:max(support(counts,p,beta,True) for p in source.values())-1
    if f(0)>=0: return dict(beta_supremum=None,absence_reason='No strict all-indicator advantage at zero stock charge')
    hi=1.
    while f(hi)<0:
        hi*=2;assert hi<=1024
    root=float(brentq(f,0,hi,xtol=1e-12,rtol=1e-12))
    at=envelope(counts,source,root)
    return dict(beta_supremum=root,strict_supremum_not_pass_at_root=True,
                limiting_indicator=at['limiting_indicator'],maximum_upper_ratio_at_root=at['maximum_upper_ratio'],
                capital_pair_service_equivalent=counts['delivered_services']/(root*H*counts['standing_fleet']),
                reciprocal_charge_weeks=1/root,not_observed_lifetime=True)


def scenarios():
    rows=[]
    for weights,c in [('uniform',0.),('uniform',.5),('uniform',1.),('ramp',0.)]:
        for tau in [1,3,5]: rows.append(dict(name=f'main_{weights}_c{c:g}_tau{tau}',weights=weights,compression=c,tau=tau,q=1.,origin=0.))
    for q in [.94,.97,.99]:
        for tau in [1,3]: rows.append(dict(name=f'loss_q{q:g}_tau{tau}',weights='uniform',compression=0.,tau=tau,q=q,origin=0.))
    for origin in [0.,math.pi/2]: rows.append(dict(name=f'origin_control_{origin:.6f}',weights='uniform',compression=.35,tau=3,q=1.,origin=origin))
    assert len(rows)==20
    return rows


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--package-root',type=Path,required=True);ap.add_argument('--output-dir',type=Path,required=True)
    a=ap.parse_args();root=a.package_root.resolve();out=a.output_dir.resolve()
    if False: raise ValueError('Computational outputs must stay D')
    out.mkdir(parents=True,exist_ok=True)
    assert (out/'analysis_protocol.md').is_file() and not (out/'periodic_results.json').exists()
    began=time.monotonic();inputs=root/'inputs'
    source=source_polytopes(inputs,inputs/'replacement_stage_facts_v21.json')
    _,_,distance=cached_geometry(inputs/'geometry',out)
    netmin={}
    for k,p in source.items():
        opt=linprog([1,1,0,0,0],A_ub=p['A'],b_ub=p['b'],bounds=[(None,None)]*5,method='highs')
        assert opt.success and opt.fun>0;netmin[k]=float(opt.fun)
    outputs=[];flows=[];comparison=[]
    for s in scenarios():
        d,w,ph,chi=profile(s['weights'],s['compression'],s['origin'])
        pair={}
        for mode in ['dedicated_adaptive','pooled']:
            lp=PeriodicLP(d,distance,s['tau'],s['q'],mode);result,x=lp.solve()
            result.update(scenario=s,period_weeks=H,calendar_demands=d.tolist(),weights=w.tolist(),phases=ph.tolist(),
                          phase_coherence=chi,nominal_beta_cost=BETA_COST,capacity_per_washer_week=2000.,asset_slots=lp.slots)
            if result['success']:
                result['beta_evaluations']=[envelope(result['counts'],source,b) for b in BETAS]
                result['inverse_stock_charge']=threshold(result['counts'],source)
                for key,v in lp.index.items():
                    if x[v]>1e-7: flows.append(dict(case=s['name'],mode=mode,variable=key[0],coordinates=':'.join(map(str,key[1:])),value=float(x[v])))
            pair[mode]=result;outputs.append(result)
            print(f"{s['name']} {mode} success={result['success']}",flush=True)
        D,P=pair['dedicated_adaptive'],pair['pooled']
        row=dict(case=s['name'],weights=s['weights'],compression=s['compression'],tau=s['tau'],q=s['q'],origin=s['origin'],
                 dedicated_success=D['success'],pooled_success=P['success'])
        if D['success'] and P['success']:
            row.update(dedicated_fleet=D['counts']['standing_fleet'],pooled_fleet=P['counts']['standing_fleet'],
                       partial_cost_saving_fraction=1-P['cost_eur']/D['cost_eur'],
                       dedicated_beta_supremum=D['inverse_stock_charge'].get('beta_supremum'),
                       pooled_beta_supremum=P['inverse_stock_charge'].get('beta_supremum'))
            for label,res in [('dedicated',D),('pooled',P)]:
                for tag,v in zip(['zero','52week','12week'],res['beta_evaluations']):
                    row[f'{label}_{tag}_max_upper']=v['maximum_upper_ratio'];row[f'{label}_{tag}_outcome']=v['outcome']
        comparison.append(row)
    dump_csv(out/'periodic_positive_flows.csv',flows);dump_csv(out/'periodic_comparisons.csv',comparison)
    (out/'periodic_results.json').write_text(json.dumps(dict(outputs=outputs,source_net_capital_minima=netmin,
        scope='Periodic source-response sensitivity; beta is an attributed capital-stock charge, not measured survival or a free-stock lifetime LCA'),indent=2)+'\n',encoding='utf-8')
    summary=dict(actual_nominal_periodic_solves=len(outputs),feasible_policies=sum(r['success'] for r in outputs),
                 paired_cases=len(comparison),source_indicators=len(source),positive_flow_rows=len(flows),
                 endpoint_evaluations=sum(2*14*len(r.get('beta_evaluations',[])) for r in outputs),
                 all_results_retained=True,empirical_lifetime=False,elapsed_seconds=time.monotonic()-began,
                 independent_verification='separate report required')
    (out/'periodic_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8');print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
