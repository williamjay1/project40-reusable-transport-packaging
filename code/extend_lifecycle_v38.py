"""Structural and state-envelope checks; all inputs explicitly synthetic."""
import itertools
import json
import math
import hashlib
from pathlib import Path
from collections import Counter

from run_lifecycle_v38 import Job, path_masks, subset_dp

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/lifecycle_v38/extension"


def formula_test():
    rows=[]
    for N,h,K,L,delta,tau,mode in itertools.product(
        (3,6,9),(1,2,3),(2,3,5),(2,4,8),(1,2),(1,2),
        ("pooled","dedicated")):
        if tau>delta or L<tau: continue
        jobs=[Job(i,i%h,i*delta,i*delta+tau) for i in range(N)]
        m=min(K,1+(L-tau)//(delta if mode=="pooled" else h*delta))
        expected=math.ceil(N/m) if mode=="pooled" else sum(
            math.ceil(sum(j.owner==o for j in jobs)/m) for o in range(h))
        paths=path_masks(jobs,K,L,mode)
        exact=subset_dp(paths,[],N)
        assert expected==exact,(N,h,K,L,delta,tau,mode,expected,exact)
        # Construct contiguous shared blocks or contiguous owner blocks.
        groups=[list(range(N))] if mode=="pooled" else [
            [i for i,j in enumerate(jobs) if j.owner==o] for o in range(h)]
        witnesses=[g[i:i+m] for g in groups for i in range(0,len(g),m)]
        feasible={p["mask"] for p in paths}
        assert all(sum(1<<i for i in p) in feasible for p in witnesses)
        assert len(witnesses)==exact
        rows.append(dict(N=N,owners=h,K=K,L=L,spacing=delta,delay=tau,
                         ownership=mode,m=m,exact_births=exact,formula=expected))
    return rows


def old_box_test():
    rows=[]
    for mode in ("pooled","dedicated"):
        endpoint_values=[]
        mode_rows=[]
        for a0,u0,a1,u1,tau in itertools.product((0,2,4),(0,1,2),(0,2,4),(0,1,2),(1,2)):
            jobs=[Job(i,0 if i<2 else 1,t,t+tau) for i,t in enumerate((0,1,4,5))]
            new=path_masks(jobs,3,8,mode)
            old=[path_masks(jobs,3,8,mode,owner,age,uses)
                 for owner,age,uses in ((0,a0,u0),(1,a1,u1))]
            value=subset_dp(new,old,4)
            row=dict(ownership=mode,ages=[a0,a1],previous_uses=[u0,u1],
                     delay=tau,new_births=value)
            mode_rows.append(row)
            if (a0,u0,a1,u1,tau) in ((0,0,0,0,1),(4,2,4,2,2)):
                endpoint_values.append(((a0,u0,a1,u1,tau),value))
        lower=next(v for s,v in endpoint_values if s==(0,0,0,0,1))
        upper=next(v for s,v in endpoint_values if s==(4,2,4,2,2))
        assert min(r["new_births"] for r in mode_rows)==lower
        assert max(r["new_births"] for r in mode_rows)==upper
        assert all(lower<=r["new_births"]<=upper for r in mode_rows)
        rows.append(dict(ownership=mode,lower=lower,upper=upper,
                         endpoint_extremes_sharp=True,states=mode_rows))
    return rows


def summary_pilot():
    pilot=json.loads((ROOT/"results/lifecycle_v38/pilot/pilot_results.json").read_text())
    rows=pilot["results"]; ix={r["name"]:r for r in rows}
    gains=[]; aging=[]
    for r in rows:
        if r["ownership"]!="dedicated":continue
        pooled=ix[r["name"].replace("_dedicated","_pooled")]
        gain=r["introductions"]-pooled["introductions"]
        gains.append(gain)
        if "_young_" in r["name"]:
            old_d=ix[r["name"].replace("_young_","_old_")]
            old_p=ix[old_d["name"].replace("_dedicated","_pooled")]
            aging.append(dict(young_dedicated=r["name"],young_gain=gain,
                              old_gain=old_d["introductions"]-old_p["introductions"],
                              young_births=[r["introductions"],pooled["introductions"]],
                              old_births=[old_d["introductions"],old_p["introductions"]]))
    return dict(pairs=len(gains),positive_sharing_pairs=sum(x>0 for x in gains),
                equality_pairs=sum(x==0 for x in gains),gain_counts=dict(Counter(gains)),
                young_to_old_comparisons=aging,
                gaps=[{k:r[k] for k in ("name","introductions","certificate","selected_paths")}
                      for r in rows if r["certificate_gap"]>0])


def main():
    OUT.mkdir(exist_ok=True,parents=True)
    path=OUT/"structural_results.json"
    assert not path.exists()
    formulas=formula_test(); box=old_box_test(); pilot=summary_pilot()
    result=dict(all_passed=True,declared_synthetic=True,
        protocol_sha256=hashlib.sha256((ROOT/"results/lifecycle_v38/extension_protocol.md").read_bytes()).hexdigest(),
        formula_checks=len(formulas),state_box_checks=sum(len(r["states"]) for r in box),
        formula_rows=formulas,state_boxes=box,pilot_summary=pilot)
    with path.open("x",encoding="utf-8") as stream:
        json.dump(result,stream,indent=2);stream.write("\n")
    print(json.dumps({k:v for k,v in result.items() if k not in ("formula_rows","state_boxes","pilot_summary")},indent=2))
    print("Sharp conditional introduction intervals:",[(r["ownership"],r["lower"],r["upper"]) for r in box])
    print("Ownership pairs:",pilot["pairs"],"strict improvements:",pilot["positive_sharing_pairs"],"equalities:",pilot["equality_pairs"])
    print("Aging removes strict advantage:",[r for r in pilot["young_to_old_comparisons"] if r["young_gain"]>0 and r["old_gain"]==0][:6])


if __name__=="__main__": main()
