"""Exact finite-asset path pilot, independent subset DP and rational cuts."""
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from fractions import Fraction
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, linprog, milp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/lifecycle_v38/pilot"


@dataclass(frozen=True)
class Job:
    id: int
    owner: int
    start: int
    ready: int
    lost: bool = False


def path_masks(jobs, K, L, ownership, old_owner=None, age=0, previous=0):
    """Every path obeys expiry through final ready time; losses are terminal."""
    paths = []
    order = sorted(range(len(jobs)), key=lambda i:(jobs[i].start, jobs[i].id))
    def extend(chosen, after):
        if chosen:
            paths.append(tuple(chosen))
            if jobs[chosen[-1]].lost:
                return
        if len(chosen) >= K - previous:
            return
        for position in range(after,len(order)):
            index=order[position]; job=jobs[index]
            if chosen and job.start < jobs[chosen[-1]].ready:
                continue
            if ownership == "dedicated":
                required = old_owner if old_owner is not None else (jobs[chosen[0]].owner if chosen else job.owner)
                if job.owner != required:
                    continue
            birth = -age if old_owner is not None else (jobs[chosen[0]].start if chosen else job.start)
            if job.ready > birth + L:
                continue
            extend(chosen+[index], position+1)
    extend([],0)
    return [{"jobs":p,"mask":sum(1<<i for i in p)} for p in paths]


def subset_dp(new, olds, n):
    """Separate exhaustive coverage DP, no optimizer/duals."""
    total=(1<<n)-1
    best=[math.inf]*(total+1); best[0]=0
    for paths in olds:
        later=best.copy()
        for mask,val in enumerate(best):
            if not math.isfinite(val): continue
            for p in paths:
                if not(mask & p["mask"]):
                    merged=mask|p["mask"]
                    later[merged]=min(later[merged],val)
        best=later
    # First-uncovered-job branching avoids duplicate partitions.
    for mask in range(total+1):
        if not math.isfinite(best[mask]) or mask==total: continue
        first=next(i for i in range(n) if not(mask>>i & 1))
        for p in new:
            if p["mask"]>>first & 1 and not(mask & p["mask"]):
                merged=mask|p["mask"]
                best[merged]=min(best[merged],best[mask]+1)
    return int(best[total]) if math.isfinite(best[total]) else None


def weighted_certificate(new, olds, n):
    """Standard path-cover dual, exact rational feasibility after LP proposal."""
    rows=[]; rhs=[]
    for p in new:
        rows.append([int(p["mask"]>>i & 1) for i in range(n)]+[0]*len(olds)); rhs.append(1.)
    for resource, paths in enumerate(olds):
        for p in paths:
            row=[int(p["mask"]>>i & 1) for i in range(n)]+[0]*len(olds)
            row[n+resource]=-1
            rows.append(row);rhs.append(0.)
    objective=np.r_[-np.ones(n),np.ones(len(olds))]
    lp=linprog(objective,A_ub=np.array(rows),b_ub=rhs,bounds=(0,None),method="highs")
    if not lp.success:
        return None
    w=[Fraction(float(max(0,x))).limit_denominator(100000) for x in lp.x[:n]]
    peak=max((sum(w[i] for i in p["jobs"]) for p in new),default=Fraction(0))
    if peak>1: w=[x/peak for x in w]
    allowances=[max((sum(w[i] for i in p["jobs"]) for p in paths),default=Fraction(0)) for paths in olds]
    exact=sum(w)-sum(allowances)
    assert all(sum(w[i] for i in p["jobs"])<=1 for p in new)
    assert all(sum(w[i] for i in p["jobs"])<=allowances[a] for a,paths in enumerate(olds) for p in paths)
    ceil=max(0,-(-exact.numerator//exact.denominator))
    return dict(weights=[str(x) for x in w],old_allowances=[str(x) for x in allowances],
                value_fraction=str(exact),value=float(exact),ceil=ceil,
                all_path_constraints_exactly_verified=True)


def solve_case(name,jobs,K,L,ownership,state):
    n=len(jobs)
    new=path_masks(jobs,K,L,ownership)
    initial=[];olds=[]
    if state!="none":
        for owner in range(3):
            age=0 if state=="young" else L-3
            previous=0 if state=="young" else K-1
            initial.append(dict(id=owner,owner=owner,age=age,previous_services=previous))
            olds.append(path_masks(jobs,K,L,ownership,owner,age,previous))
    columns=[dict(p,old=None) for p in new]
    for a,paths in enumerate(olds):
        columns.extend(dict(p,old=a) for p in paths)
    A=np.zeros((n+len(olds),len(columns)))
    for col,p in enumerate(columns):
        for i in p["jobs"]:A[i,col]=1
        if p["old"] is not None:A[n+p["old"],col]=1
    costs=np.array([int(p["old"] is None) for p in columns],dtype=float)
    lower=np.r_[np.ones(n),np.zeros(len(olds))]
    upper=np.ones(n+len(olds))
    result=milp(costs,integrality=np.ones(len(columns)),bounds=Bounds(0,1),
                constraints=LinearConstraint(A,lower,upper),
                options={"time_limit":30,"mip_rel_gap":0})
    dp=subset_dp(new,olds,n)
    certificate=weighted_certificate(new,olds,n)
    if not result.success:
        assert dp is None
        return dict(name=name,K=K,L=L,ownership=ownership,state=state,
                    feasible=False,dp_optimum=dp,solver_status=int(result.status))
    selected=[p for p,x in zip(columns,result.x) if x>.5]
    coverage=[sum(i in p["jobs"] for p in selected) for i in range(n)]
    old_ids=[p["old"] for p in selected if p["old"] is not None]
    births=sum(p["old"] is None for p in selected)
    assert coverage==[1]*n and len(old_ids)==len(set(old_ids))
    assert births==dp and certificate["ceil"]<=births
    peak=max(sum(job.start<=t<job.ready for job in jobs) for t in range(max(j.ready for j in jobs)))
    totalusebound=max(0,math.ceil((n-sum(K-a["previous_services"] for a in initial))/K))
    return dict(name=name,K=K,L=L,ownership=ownership,state=state,feasible=True,
                jobs=[asdict(j) for j in jobs],initial_assets=initial,
                introductions=births,dp_optimum=dp,
                certificate=certificate,certificate_gap=births-certificate["ceil"],
                selected_paths=[dict(jobs=[jobs[i].id for i in p["jobs"]],old=p["old"]) for p in selected],
                service_coverage=coverage,job_count=n,occupied_peak=peak,
                simple_global_use_bound=totalusebound,
                ready_transitions=sum(len(p["jobs"])-1 for p in selected),
                loss_terminal_jobs=sum(j.lost for j in jobs),
                path_count=len(columns),solver_status=int(result.status),
                mip_gap=float(result.mip_gap),minimum_births_proved=True)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    assert not (OUT/"pilot_results.json").exists()
    protocol=ROOT/"results/lifecycle_v38/pilot_protocol.md"
    started=datetime.now(timezone.utc).isoformat(); tick=time.perf_counter()
    results=[]
    for timing,delay,K,L,loss,state,ownership in itertools.product(
        ("synchronous","staggered"),(2,4),(2,3,5,9),(4,8,12),
        (False,True),("none","young","old"),("dedicated","pooled")):
        jobs=[]
        for cycle in range(3):
            for owner in range(3):
                index=3*cycle+owner
                start=3*cycle if timing=="synchronous" else index
                jobs.append(Job(index,owner,start,start+delay,loss and index in (2,6)))
        name=f"{timing}_delay{delay}_K{K}_L{L}_loss{int(loss)}_{state}_{ownership}"
        row=solve_case(name,jobs,K,L,ownership,state)
        row.update(timing=timing,delay=delay,loss=loss)
        results.append(row)
    assert len(results)==576
    extra=[Job(0,0,0,8)]+[Job(i,1+(i%2),i,i+1) for i in range(1,7)]
    for ownership in ("dedicated","pooled"):
        results.append(solve_case("long_job_counterexample_"+ownership,extra,4,20,ownership,"none"))
    index={row["name"]:row for row in results}
    shared_monotonic=[]
    for row in results:
        if row["ownership"]!="dedicated":continue
        partner=index[row["name"].rsplit("_",1)[0]+"_pooled"]
        ok=(not row["feasible"]) or (partner["feasible"] and partner["introductions"]<=row["introductions"])
        assert ok
        shared_monotonic.append(ok)
    budget_rows=[]
    for row in results:
        if not row["feasible"]:continue
        for env,cash in itertools.product((1.5,2.5,3.5,4.5,5.5,8.5),repeat=2):
            budget_rows.append(dict(case=row["name"],births=row["introductions"],
                environmental_birth_allowance=env,cash_birth_allowance=cash,
                necessary_certificate_passes=row["certificate"]["ceil"]<=min(env,cash),
                constructive_joint_feasible=row["introductions"]<=min(env,cash),
                budgets_are_declared_asset_equivalents=True))
    document=dict(all_passed=True,declared_synthetic_scenarios=True,
                  protocol_sha256=hashlib.sha256(protocol.read_bytes()).hexdigest(),
                  started_at_utc=started,finished_at_utc=datetime.now(timezone.utc).isoformat(),
                  elapsed_seconds=time.perf_counter()-tick,
                  factorial_policies=576,extra_counterexample_policies=2,
                  independent_DP_matches_every_feasible_policy=True,
                  rational_certificate_valid_every_path=True,
                  ownership_inclusion_checked_pairs=len(shared_monotonic),
                  results=results)
    with (OUT/"pilot_results.json").open("x",encoding="utf-8") as h:json.dump(document,h,indent=2);h.write("\n")
    with (OUT/"normalized_joint_budget_grid.csv").open("x",encoding="utf-8",newline="") as h:
        writer=csv.DictWriter(h,fieldnames=list(budget_rows[0]));writer.writeheader();writer.writerows(budget_rows)
    summary={key:value for key,value in document.items() if key!="results"}
    summary.update(feasible_policies=sum(r["feasible"] for r in results),
                   tight_certificate_policies=sum(r.get("certificate_gap")==0 for r in results),
                   positive_certificate_gaps=sum(r.get("certificate_gap",0)>0 for r in results),
                   normalized_budget_rows=len(budget_rows))
    with (OUT/"pilot_summary.json").open("x",encoding="utf-8") as h:json.dump(summary,h,indent=2);h.write("\n")
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    main()
