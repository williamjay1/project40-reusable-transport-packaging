"""Public-input, conditional calendar RPC operating network.

Distances are a counterfactual geometry from R3PACK primary-packaging nodes,
NOT observations on RPC demand/flows. Calendar demand, delay, capacity and
collection are declared engineering scenarios. Source LCIA stage coefficients
are generalized under explicitly age-invariant source accounting, not rebuilt
ecoinvent LCIA or a representative EU-27 LCA. Costs are source model priors.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from dataclasses import dataclass, asdict
from pathlib import Path
import numpy as np
from openpyxl import load_workbook
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, vstack, csr_matrix
from couple_pool_budget_and_cohort_v20 import read_csv, COMPARE, FIXED
from extend_rejection_operating_certificate_v21 import interval

R0 = .0055

def dump_csv(path, rows):
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

def extract_geometry(raw, out):
    before = hashlib.sha256(raw.read_bytes()).hexdigest()
    wb = load_workbook(raw, read_only=True, data_only=True)
    sorters = [r[0] for r in list(wb['Sorter'].values)[1:]]
    washers = [r[0] for r in list(wb['Washer'].values)[1:]]
    records = [dict(origin=r[0], destination=r[1], distance_km=float(r[2]))
               for r in list(wb['Distances'].values)[1:]
               if r[0] in sorters and r[1] in washers]
    transport = list(wb['Transport'].values)
    wb.close()
    assert len(records) == len(sorters)*len(washers) == 490
    assert len({(r['origin'], r['destination']) for r in records}) == 490
    assert all(np.isfinite(r['distance_km']) and r['distance_km'] >= 0 for r in records)
    assert before == hashlib.sha256(raw.read_bytes()).hexdigest()
    dump_csv(out/'geometry_numeric_facts.csv', records)
    meta = dict(source_doi='10.5281/zenodo.21392137', source_sha256=before,
                source_sheet='Distances', sorters=sorters, washers=washers,
                directed_edges=490, reverse_edge_status='not reported; symmetric clean-leg distance is an explicit scenario',
                no_rpc_demand_transfer=True, no_capacity_or_LCIA_transfer=True,
                transport_header=list(transport[0]), transport_row=list(transport[1]),
                license='Zenodo item publicly downloadable; license not displayed; workbook not redistributed')
    (out/'geometry_provenance.json').write_text(json.dumps(meta, indent=2)+'\n', encoding='utf-8')
    d = np.array([[next(r['distance_km'] for r in records if r['origin']==j and r['destination']==i)
                   for i in washers] for j in sorters])
    return sorters, washers, d

def cached_geometry(cache, out):
    meta=json.loads((cache/'geometry_provenance.json').read_text(encoding='utf-8'))
    records=read_csv(cache/'geometry_numeric_facts.csv')
    d=np.array([[float(next(r['distance_km'] for r in records if r['origin']==j and r['destination']==i))
                 for i in meta['washers']] for j in meta['sorters']])
    assert d.shape==(14,35) and len(records)==490
    # Derived cache, not new/raw material; preserve original source provenance.
    dump_csv(out/'geometry_numeric_facts.csv',records)
    # Preserve exact provenance bytes: reserialization can change platform newlines
    # and would falsely break the downstream cross-stage fingerprint requirement.
    (out/'geometry_provenance.json').write_bytes((cache/'geometry_provenance.json').read_bytes())
    return meta['sorters'],meta['washers'],d

def source_polytopes(input_dir, stage_facts):
    grid = read_csv(input_dir/'tua_reported_grid.csv')
    contrib = read_csv(input_dir/'tua_process_contributions.csv')
    stage = json.loads(stage_facts.read_text(encoding='utf-8'))['facts']
    result = {}
    for row in stage:
        k = row['indicator']
        gs = [g for g in grid if g['source_scenario']=='W1_RE1' and g['indicator']==k]
        u = float(gs[0]['single_use_impact_per_service_100_crates'])
        vectors=[]; centres=[]; widths=[]
        for g in gs:
            n=int(g['reported_rotations']); y=float(g['reusable_cumulative_impact_100_crates'])
            h=FIXED[k]/2 if k in FIXED else 10**(np.floor(np.log10(y))-2)/2
            vectors.append([1,1,n-1,n-1,n-1]); centres.append(y/u); widths.append(h/u)
        for n in [20,125]:
            for tag, v in [('P',[1,0,n-1,0,0]),('E',[0,1,0,n-1,0]),('W',[0,0,0,0,n-1])]:
                cell=row['cells'][f'{tag}{n}']
                vectors.append(v); centres.append(cell['value']/u)
                widths.append((cell['upper']-cell['lower'])/(2*u))
        a=np.asarray(vectors,float); c=np.asarray(centres); h=np.asarray(widths)
        A=np.vstack([a,-a]); b=np.r_[c+h,-c+h]
        # Fractional interval Chebyshev centre; all source correlations retained.
        mid=linprog([0,0,0,0,0,-1], A_ub=np.c_[A,np.r_[h,h]], b_ub=b,
                    bounds=[(None,None)]*5+[(0,1)], method='highs')
        assert mid.success,(k,mid.message)
        process={}
        for p in [r for r in contrib if r['indicator']==k]:
            v,hw=interval(p['displayed_value'])
            process[p['process']] = (float(v-hw)/u, float(v+hw)/u)
        result[k]=dict(A=A,b=b,centre=mid.x[:5],u=u,
                       u_lower_ratio=1-COMPARE[k]/(2*u),
                       u_upper_ratio=1+COMPARE[k]/(2*u), process=process,
                       maximum_fractional_rounding_slack=float(mid.x[-1]))
    assert len(result)==14
    return result

@dataclass
class Scenario:
    name: str
    weeks: int = 12
    mean_weekly_demand: float = 35000.
    seasonal_amplitude: float = 0.
    seasonal_phase: str = 'aligned'
    seasonal_period_weeks: int = 12
    return_delay_weeks: int = 1
    collection_fraction: float = 1.
    wash_capacity_ratio: float = 2.
    total_rejection_fraction: float = R0  # source combined baseline; any increase is post-process
    transport_burden_factor: float = 1.
    electricity_burden_factor: float = 1.
    heating_burden_factor: float = 1.
    clean_leg_factor: float = 1.
    cost_per_truck_km: float = .5
    declared_empty_crates_per_truck: float = 1000.
    manufacturing_eur_per_crate: float = 4.55
    wash_eur_per_crate: float = .07
    holding_eur_per_crate_week: float = .03

def demand(s, J):
    phases=np.zeros(J) if s.seasonal_phase=='aligned' else np.arange(J)*2*np.pi/J
    d=s.mean_weekly_demand/J*(1+s.seasonal_amplitude*np.sin(
        2*np.pi*np.arange(s.weeks)[:,None]/s.seasonal_period_weeks+phases[None,:]))
    assert np.all(d>=0)
    assert abs(d.sum()-s.weeks*s.mean_weekly_demand)<1e-6
    return d

class CalendarLP:
    def __init__(self,s,distance,mode,*,demands=None):
        self.s=s; self.distance=distance; self.mode=mode
        H=s.weeks; J,I=distance.shape; self.H,self.J,self.I=H,J,I
        self.dem=demand(s,J) if demands is None else np.asarray(demands,dtype=float).copy()
        if self.dem.shape!=(H,J) or not np.all(np.isfinite(self.dem)) or np.any(self.dem<0):
            raise ValueError('Explicit calendar demand must be finite, nonnegative and H by J')
        if abs(float(self.dem.sum())-H*s.mean_weekly_demand)>1e-6:
            raise ValueError('Explicit demand must preserve the declared total service')
        # v36 source-quality restriction: no unverified exact-zero movement.
        # Keep original numeric facts and per-site capacities unchanged.
        if not np.all(np.isfinite(distance)) or np.any(distance<0):
            raise ValueError('Distance facts must be finite and nonnegative')
        self.allowed=distance>0
        if not np.all(self.allowed.any(axis=1)):
            raise ValueError('Each customer requires a positive-distance route')
        self.inactive_sites=set(np.flatnonzero(~self.allowed.any(axis=0)).tolist())
        self.nearest=np.argmin(np.where(self.allowed,distance,np.inf),axis=1)
        self.fractions=None; self.slot_pairs=None
        if mode=='dedicated_adaptive':
            self.slot_pairs=[(j,i) for j in range(J) for i in range(I)]
            self.S=len(self.slot_pairs)
        elif mode=='dedicated_static':
            # Strong static baseline: globally minimum-distance assignment with
            # steady-state wash capacity, NOT a deliberately overloaded nearest
            # heuristic. Each positive assigned route has its own asset stock.
            ae=np.kron(np.eye(J),np.ones((1,I)))
            au=np.tile(np.eye(I),(1,J))
            q=linprog(distance.ravel(),A_ub=au,
                      b_ub=np.repeat(s.mean_weekly_demand*s.wash_capacity_ratio/I,I),
                      A_eq=ae,b_eq=np.repeat(s.mean_weekly_demand/J,J),
                      bounds=[(0,None) if allowed else (0,0) for allowed in self.allowed.ravel()],method='highs')
            assert q.success,q.message
            self.fractions=q.x.reshape(J,I)/(s.mean_weekly_demand/J)
            self.slot_pairs=[(j,i) for j in range(J) for i in range(I) if self.fractions[j,i]>1e-10]
            self.S=len(self.slot_pairs)
        else:
            self.S=J if mode=='dedicated_nearest' else I
        self.names=[]; self.index={}
        def add(tag,shape):
            for p in np.ndindex(shape):
                key=(tag,)+p; self.index[key]=len(self.names); self.names.append(key)
        add('initial',(self.S,)); add('ready',(H,self.S)); add('new',(H,self.S))
        add('dispatch',(H,J,I)); add('wash',(H,J,I)); add('dirty',(H,J))
        N=len(self.names); self.N=N; self.cost=np.zeros(N); self.bounds=[(0,None)]*N
        idx=self.index
        for slot in range(self.S): self.cost[idx['initial',slot]]=s.manufacturing_eur_per_crate
        transport_cost=s.cost_per_truck_km/s.declared_empty_crates_per_truck
        for key,p in idx.items():
            if key[0] in ['ready','dirty']: self.cost[p]=s.holding_eur_per_crate_week
            if key[0]=='new': self.cost[p]=s.manufacturing_eur_per_crate
            if key[0]=='wash': self.cost[p]=s.wash_eur_per_crate+s.total_rejection_fraction*s.manufacturing_eur_per_crate+transport_cost*distance[key[2],key[3]]
            if key[0]=='dispatch': self.cost[p]=transport_cost*s.clean_leg_factor*distance[key[2],key[3]]
            if key[0] in ['dispatch','wash'] and mode=='dedicated_nearest' and key[3]!=self.nearest[key[2]]:
                self.bounds[p]=(0,0)
            if key[0] in ['dispatch','wash'] and mode=='dedicated_static' and self.fractions[key[2],key[3]]<=1e-10:
                self.bounds[p]=(0,0)
            if key[0] in ['dispatch','wash'] and not self.allowed[key[2],key[3]]:
                self.bounds[p]=(0,0)
            if key[0] in ['initial','ready','new'] and mode!='dedicated_nearest':
                slot=key[-1]
                site=self.slot_pairs[slot][1] if self.slot_pairs is not None else slot
                if site in self.inactive_sites:
                    self.bounds[p]=(0,0)
        eq=[]; erhs=[]; ub=[]; urhs=[]
        def equation(terms,rhs): eq.append(terms); erhs.append(rhs)
        tau=s.return_delay_weeks; q=s.collection_fraction
        for t in range(H):
            for j in range(J):
                equation([(idx['dispatch',t,j,i],1) for i in range(I)],self.dem[t,j])
                if mode=='dedicated_static':
                    for i in range(I):
                        if self.fractions[j,i]>1e-10:
                            equation([(idx['dispatch',t,j,i],1)],self.dem[t,j]*self.fractions[j,i])
                terms=[(idx['dirty',t,j],1)]+[(idx['wash',t,j,i],1) for i in range(I)]
                if t: terms.append((idx['dirty',t-1,j],-1))
                equation(terms,q*self.dem[t-tau,j] if t>=tau else 0.)
            for slot in range(self.S):
                terms=[(idx['ready',t,slot],1),(idx['new',t,slot],-1)]
                terms.append((idx['ready',t-1,slot],-1) if t else (idx['initial',slot],-1))
                if mode=='dedicated_nearest':
                    terms += [(idx['dispatch',t,slot,i],1) for i in range(I)]
                    if t: terms += [(idx['wash',t-1,slot,i],-1) for i in range(I)]
                elif mode in ['dedicated_static','dedicated_adaptive']:
                    j,i=self.slot_pairs[slot]
                    terms.append((idx['dispatch',t,j,i],1))
                    if t: terms.append((idx['wash',t-1,j,i],-1))
                else:
                    terms += [(idx['dispatch',t,j,slot],1) for j in range(J)]
                    if t: terms += [(idx['wash',t-1,j,slot],-1) for j in range(J)]
                equation(terms,0.)
            # Replenish only losses detected by the declared return deadline.
            # Newly supplied crates do not create end-of-life credits for losses.
            if mode=='dedicated_nearest':
                for j in range(J): equation([(idx['new',t,j],1)], (1-q)*self.dem[t-tau,j] if t>=tau else 0.)
            elif mode in ['dedicated_static','dedicated_adaptive']:
                for slot,(j,i) in enumerate(self.slot_pairs):
                    if mode=='dedicated_static':
                        equation([(idx['new',t,slot],1)],(1-q)*self.dem[t-tau,j]*self.fractions[j,i] if t>=tau else 0.)
                if mode=='dedicated_adaptive':
                    for j in range(J):
                        equation([(idx['new',t,slot],1) for slot,(jj,i) in enumerate(self.slot_pairs) if jj==j],(1-q)*self.dem[t-tau,j] if t>=tau else 0.)
            else:
                equation([(idx['new',t,i],1) for i in range(I)],(1-q)*self.dem[t-tau].sum() if t>=tau else 0.)
            for i in range(I):
                ub.append([(idx['wash',t,j,i],1) for j in range(J)])
                urhs.append(s.mean_weekly_demand*s.wash_capacity_ratio/I)
        def sparse(rows):
            ri=[];ci=[];v=[]
            for r, terms in enumerate(rows):
                for c,x in terms: ri.append(r);ci.append(c);v.append(x)
            return coo_matrix((v,(ri,ci)),shape=(len(rows),N)).tocsr()
        self.Aeq=sparse(eq); self.beq=np.array(erhs)
        self.Aub=sparse(ub); self.bub=np.array(urhs)
        self.total_services=float(self.dem.sum())
        self.total_missing=(1-q)*self.total_services

    def coefficients(self,theta,process):
        s=self.s; p0,e0,p,e,w=theta
        tt,te,th=process
        a=np.zeros(self.N)
        for key,ix in self.index.items():
            tag=key[0]
            if tag in ['initial','new']: a[ix]=(p0+e0)/100
            if tag=='wash':
                a[ix]=(s.total_rejection_fraction/R0*(p+e)+w-tt+
                       tt*s.transport_burden_factor*self.distance[key[2],key[3]]/139+
                       (s.electricity_burden_factor-1)*te+(s.heating_burden_factor-1)*th)/100
            if tag=='dispatch':
                a[ix]=s.clean_leg_factor*tt*s.transport_burden_factor*self.distance[key[2],key[3]]/(100*139)
        constant=-self.total_missing*e0/100
        return a,constant

    def counts(self,x):
        sums={tag:sum(x[ix] for key,ix in self.index.items() if key[0]==tag)
              for tag in ['initial','new','ready','wash','dispatch','dirty']}
        V=sums['initial']+sums['new']; T=V-self.total_missing
        R=self.s.total_rejection_fraction*sums['wash']; W=sums['wash']
        dirty_distance=sum(x[ix]*self.distance[key[2],key[3]] for key,ix in self.index.items() if key[0]=='wash')
        clean_distance=sum(x[ix]*self.distance[key[2],key[3]] for key,ix in self.index.items() if key[0]=='dispatch')
        assert T>=-1e-7
        return dict(virgin_assets=V,terminal_recovered_assets=T,closed_loop_rejections=R,
                    wash_inputs=W,initial_assets=sums['initial'],replacement_virgin=sums['new'],
                    ready_crate_weeks=sums['ready'],dirty_crate_weeks=sums['dirty'],
                    delivered_services=sums['dispatch'],missing_assets=self.total_missing,
                    empty_dirty_crate_km=dirty_distance,empty_clean_crate_km=clean_distance,
                    mean_dirty_distance_km=dirty_distance/W if W else 0,
                    initial_fleet_service_utilization=self.total_services/sums['initial'],
                    source_mass_recycled_kg=1.49*(T+R), source_washing_water_m3=.00055*W,
                    source_washing_electricity_kwh=.0248*W,
                    source_washing_heat_mj=.322*W)

    def support(self,x,source,k,upper=True):
        p=source[k]; c=self.counts(x)
        vec=np.array([c['virgin_assets']/100,c['terminal_recovered_assets']/100,
                      c['closed_loop_rejections']/(100*R0),c['closed_loop_rejections']/(100*R0),c['wash_inputs']/100])
        opt=linprog(-vec if upper else vec,A_ub=p['A'],b_ub=p['b'],
                    bounds=[(None,None)]*5,method='highs')
        assert opt.success,(k,opt.message)
        s=self.s
        scale=np.array([(s.transport_burden_factor*(c['empty_dirty_crate_km']+s.clean_leg_factor*c['empty_clean_crate_km'])/139-c['wash_inputs'])/100,
                        (s.electricity_burden_factor-1)*c['wash_inputs']/100,
                        (s.heating_burden_factor-1)*c['wash_inputs']/100])
        proc=[]
        for tag,g in zip(['return_transport','electricity','water_heating'],scale):
            lo,hi=p['process'][tag]
            proc.append(hi if (g>=0)==upper else lo)
        a,constant=self.coefficients(opt.x,proc)
        value=float(a@x+constant)
        denominator=self.total_services/100*(p['u_lower_ratio'] if upper else p['u_upper_ratio'])
        return value/denominator,a/denominator,constant/denominator

    def solve(self,source,robust=False,max_cuts=12):
        A=self.Aub; b=self.bub; cuts=[]
        for iteration in range(max_cuts):
            sol=linprog(self.cost,A_ub=A,b_ub=b,A_eq=self.Aeq,b_eq=self.beq,
                        bounds=self.bounds,method='highs')
            if not sol.success:
                return dict(success=False,status=sol.status,message=sol.message,cuts=len(cuts)),None
            extrema={k:(self.support(sol.x,source,k,True)[0],self.support(sol.x,source,k,False)[0]) for k in source}
            bad=[]
            for k in source:
                ratio,a,constant=self.support(sol.x,source,k,True)
                if ratio>1+1e-8: bad.append((k,a,1-constant))
            if not robust or not bad:
                residual=max(float(np.abs(self.Aeq@sol.x-self.beq).max()),float((self.Aub@sol.x-self.bub).max()))
                c=self.counts(sol.x)
                allpass=all(v[0]<1 for v in extrema.values())
                anyfail=any(v[1]>=1 for v in extrema.values())
                outcome='certified_source_advantage' if allpass else 'certified_source_failure' if anyfail else 'unresolved'
                dual=float(b@sol.ineqlin.marginals+self.beq@sol.eqlin.marginals)
                stationarity=self.cost-A.T@sol.ineqlin.marginals-self.Aeq.T@sol.eqlin.marginals-sol.lower.marginals-sol.upper.marginals
                # Finite upper bounds are zero; other upper marginals must be zero.
                assert max(sol.ineqlin.marginals)<=1e-7 and min(sol.lower.marginals)>=-1e-7
                gap=abs(float(self.cost@sol.x)-dual)/max(1,abs(dual))
                assert gap<1e-8 and np.max(np.abs(stationarity))<1e-7
                return dict(success=True,cost_eur=float(self.cost@sol.x),cost_eur_per_service=float(self.cost@sol.x/self.total_services),
                            counts=c,impact_ratios={k:dict(lower=v[1],upper=v[0]) for k,v in extrema.items()},
                            maximum_upper_ratio=max(v[0] for v in extrema.values()),
                            limiting_indicator=max(extrema,key=lambda k:extrema[k][0]),environmental_outcome=outcome,
                            maximum_flow_residual=residual,cut_iterations=iteration,cuts=len(cuts),
                            dual_objective_eur=dual,relative_primal_dual_gap=gap,
                            maximum_dual_stationarity_residual=float(np.max(np.abs(stationarity))),
                            asset_slots=self.slot_pairs if self.slot_pairs is not None else list(range(self.S)),
                            calendar_demands=self.dem.tolist(),
                            cost_components_eur=dict(manufacturing=self.s.manufacturing_eur_per_crate*(c['virgin_assets']+c['closed_loop_rejections']),
                                                     reconditioning=self.s.wash_eur_per_crate*c['wash_inputs'],
                                                     holding=self.s.holding_eur_per_crate_week*(c['ready_crate_weeks']+c['dirty_crate_weeks']),
                                                     empty_transport=self.s.cost_per_truck_km/self.s.declared_empty_crates_per_truck*(c['empty_dirty_crate_km']+self.s.clean_leg_factor*c['empty_clean_crate_km'])),
                            collection_fate_scope='zero extra missing-fate burden; no recycling credit for missing assets' if self.total_missing else 'all assets collected; source terminal treatment retained'),sol.x
            cuts+=bad
            A=vstack([A,csr_matrix(np.array([v[1] for v in bad]))],format='csr')
            b=np.r_[b,[v[2] for v in bad]]
        return dict(success=False,status='cut_limit',message='No convergence declared',cuts=len(cuts)),None

    def schedule(self,x,name):
        rows=[]
        for key,ix in self.index.items():
            if x[ix]<=1e-7: continue
            rows.append(dict(scenario=name,mode=self.mode,variable=key[0],
                             coordinates=':'.join(map(str,key[1:])),value=float(x[ix])))
        return rows

def calibrate(source,input_dir):
    grid=read_csv(input_dir/'tua_reported_grid.csv'); rows=[]
    # Zero time before return, one wash-period turnaround, source return-only leg.
    for n in [1,3,5,8,20,30,40,50,70,90,100,125]:
        s=Scenario(f'calibration_n{n}',weeks=n,mean_weekly_demand=100,
                   return_delay_weeks=0,clean_leg_factor=0,wash_capacity_ratio=2)
        model=CalendarLP(s,np.array([[139.]]),'pooled')
        result,x=model.solve(source)
        assert result['success']
        c=result['counts']
        assert abs(c['initial_assets']-100)<1e-6 and abs(c['wash_inputs']-100*(n-1))<1e-6
        for k,p in source.items():
            vals=[g for g in grid if g['source_scenario']=='W1_RE1' and g['indicator']==k and int(g['reported_rotations'])==n]
            assert len(vals)==1
            raw=float(vals[0]['reusable_cumulative_impact_100_crates'])
            a,constant=model.coefficients(p['centre'],[sum(p['process'][tag])/2 for tag in ['return_transport','electricity','water_heating']])
            v=float(a@x+constant)*p['u']
            half=FIXED[k]/2 if k in FIXED else 10**(np.floor(np.log10(raw))-2)/2
            assert abs(v-raw)<=half+1e-7*max(abs(raw),1e-20),(n,k,v,raw,half)
            rows.append(dict(services=n,indicator=k,reported=raw,reconstructed=v,
                             printed_half_width=half,inside_printed_interval=True))
    return rows

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input-dir',type=Path,required=True)
    ap.add_argument('--stage-facts',type=Path,required=True)
    geom=ap.add_mutually_exclusive_group(required=True)
    geom.add_argument('--raw-network',type=Path)
    geom.add_argument('--geometry-cache',type=Path)
    ap.add_argument('--output-dir',type=Path,required=True)
    ap.add_argument('--quick',action='store_true')
    ap.add_argument('--adaptive-benchmark',action='store_true')
    args=ap.parse_args(); out=args.output_dir; out.mkdir(parents=True,exist_ok=True)
    _,_,dist=extract_geometry(args.raw_network,out) if args.raw_network else cached_geometry(args.geometry_cache,out)
    source=source_polytopes(args.input_dir,args.stage_facts)
    calibration=calibrate(source,args.input_dir); dump_csv(out/'source_calendar_calibration.csv',calibration)
    scenarios=[Scenario('flat_tau1'),Scenario('staggered_tau1',seasonal_amplitude=.5,seasonal_phase='staggered'),
               Scenario('aligned_tau1',seasonal_amplitude=.5),
               Scenario('staggered_tau3',seasonal_amplitude=.5,seasonal_phase='staggered',return_delay_weeks=3),
               Scenario('flat_tau5',return_delay_weeks=5),
               Scenario('staggered_q98',seasonal_amplitude=.5,seasonal_phase='staggered',collection_fraction=.98)]
    if args.quick: scenarios=scenarios[:2]
    if args.adaptive_benchmark:
        scenarios=[Scenario('flat_tau1'),Scenario('staggered_tau1',seasonal_amplitude=.5,seasonal_phase='staggered'),
                   Scenario('staggered_tau3',seasonal_amplitude=.5,seasonal_phase='staggered',return_delay_weeks=3),
                   Scenario('staggered_tau1_H24',weeks=24,seasonal_amplitude=.5,seasonal_phase='staggered'),
                   Scenario('staggered_tau1_H48',weeks=48,seasonal_amplitude=.5,seasonal_phase='staggered')]
    outputs=[]; schedules=[]
    for s in scenarios:
        for mode in (['dedicated_adaptive','pooled'] if args.adaptive_benchmark else ['dedicated_nearest','dedicated_static','pooled']):
            model=CalendarLP(s,dist,mode); result,x=model.solve(source)
            print(s.name,mode,result.get('cost_eur_per_service'),result.get('maximum_upper_ratio'),flush=True)
            result.update(scenario=asdict(s),mode=mode,environmental_constraints=False)
            outputs.append(result)
            if x is not None: schedules+=model.schedule(x,s.name)
            if mode=='pooled' and result['success'] and result['maximum_upper_ratio']>1:
                constrained,cx=model.solve(source,robust=True)
                print(s.name,'robust_constrained',constrained.get('success'),constrained.get('maximum_upper_ratio'),flush=True)
                constrained.update(scenario=asdict(s),mode='pooled_robust',environmental_constraints=True)
                outputs.append(constrained)
                if cx is not None:
                    cr=model.schedule(cx,s.name)
                    for row in cr: row['mode']='pooled_robust'
                    schedules+=cr
    dump_csv(out/'calendar_positive_flows.csv',schedules)
    metadata=dict(scope='conditional public-geometry RPC operating-network scenarios; not observed RPC flows or EU-wide estimates',
                  source_lcia='Tua2019 W1-RE1 source stage polytope; age-invariant, source recycling and terminal allocation. Every recovered asset is liquidated after final delayed returns; no continuing-stock credit. R matches recycled rejects and source-recipe replacements, ready next week.',
                  source_cost='Accorsi et al.2022 SPC31:331-345, accepted manuscript section3 pp13-14; 4.55 EUR manufacture, .07 EUR wash, .03 EUR/crate-week holding. Historical model assumptions, not 2026 prices.',
                  cost_transport='R3PACK .5 EUR/truck-km reference model coefficient; 1000 crates/truck is declared scenario, not observed load',
                  excluded_costs='depot capital/land/installation, taxes, finance, labour outside published unit costs, revenues and product-distribution costs',
                  missing_fate='q<1 results omit additional unobserved lost-fate burdens but omit all lost recycling credits; report conditional screening only',
                  asset_order='return after tau weeks; full reconditioning during one further week; matched source-recipe replacements are available with the reconditioned batch; no opening secondary feed. Fractional flows are expected/scalable units, not discrete tracked assets.',
                  source_recipe='1490g finished crate; nominal .994 molding yield, .93 recycling yield; .91kg finished-equivalent secondary/.58kg virgin recipe. Recipe feed availability is checked; no extra credit assigned to surplus recyclate and no empirical material-yield claim.',
                  holding_location='historical .03 EUR/crate-week rate applied to both ready and dirty queues as a declared uniform proxy; customer queues are not the source operator plant inventory',
                  uncertainty='joint rounding polytopes are epistemic reporting envelopes, not probability distributions',
                  baseline='Primary adaptive comparison: dedicated_adaptive retains customer-owned clean stocks with all arcs, capacities and timing identical to pooling. When included, dedicated_static and dedicated_nearest are secondary diagnostics.',
                  calibration_rows=len(calibration),model_result_count=len(outputs),outputs=outputs)
    (out/'calendar_network_results.json').write_text(json.dumps(metadata,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(calibration_rows=len(calibration),model_result_count=len(outputs),positive_flow_rows=len(schedules)),indent=2))

if __name__=='__main__': main()
