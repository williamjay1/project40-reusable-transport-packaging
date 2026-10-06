"""Replacement-preserving original-cohort bounds coupled to source budgets.

Synchronous fixed-capacity service is an explicit model completion, not observed
telemetry. All source initial/replacement material and terminal flows stay fixed.
The source LCA coefficients are not interpreted as per-asset F/v.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import numpy as np

FIXED = dict(zip(["Climate change", "Photochemical ozone formation", "Acidification",
    "Terrestrial eutrophication", "Freshwater ecotoxicity", "Cumulative energy demand", "Water depletion"],
    [.01,.01,.01,.01,1.,1.,.01]))
COMPARE = dict(zip(["Climate change", "Ozone depletion", "Human toxicity, non-cancer effects",
    "Human toxicity, cancer effects", "Particulate matter", "Photochemical ozone formation",
    "Acidification", "Terrestrial eutrophication", "Freshwater eutrophication", "Marine eutrophication",
    "Freshwater ecotoxicity", "Mineral, fossil & renewable resource depletion", "Cumulative energy demand", "Water depletion"],
    [.01,1e-7,1e-7,1e-7,1e-4,1e-3,1e-3,1e-3,1e-4,1e-4,.01,1e-5,1.,.01]))

def read_csv(path):
    with path.open(encoding="utf-8-sig",newline="") as stream:
        return list(csv.DictReader(stream))

def allocation(N, rejected_each_round, B, original_first):
    """Exact integer generational ledger: entry generation -> active count."""
    active = {0:N}
    cumulative = []
    for j in range(1,B):
        remove = rejected_each_round
        order = sorted(active) if original_first else sorted(active,reverse=True)
        for generation in order:
            retired = min(remove,active[generation])
            active[generation] -= retired
            remove -= retired
            if not remove:
                break
        assert remove == 0
        active[j] = rejected_each_round
        assert sum(active.values()) == N
        cumulative.append({"transition":j,"original_survivors":active[0],
                           "replacement_survivors":N-active[0],
                           "new_replacements":rejected_each_round})
    return {"original_survivors_at_B":active[0], "active_by_entry_generation":active,
            "initial_manufacturing":N,
            "replacement_manufacturing":rejected_each_round*(B-1),
            "reconditioning_inputs":N*(B-1),
            "recycling_total_including_terminal_pool":N+rejected_each_round*(B-1),
            "delivered_services":N*B, "trace":cumulative}

def vertex_budget_coefficients(grid,contributions):
    groups,stage=defaultdict(list),defaultdict(float)
    for row in grid:
        if row["source_scenario"]=="W1_RE1":
            groups[row["indicator"]].append(row)
    for row in contributions:
        mantissa,_,exp=row["displayed_value"].upper().partition("E")
        half=.5*10**((int(exp) if exp else 0)-len(mantissa.partition(".")[2]))
        stage[row["indicator"]]+=float(row["impact_per_100_crates_reconditioned"])+half
    assert len(groups)==len(stage)==14
    result={}
    for k,rows in groups.items():
        u=float(rows[0]["single_use_impact_per_service_100_crates"])
        x=np.array([int(row["reported_rotations"]) for row in rows],float)
        raw=np.array([float(row["reusable_cumulative_impact_100_crates"]) for row in rows])
        h=np.array([FIXED[k]/2 if k in FIXED else 10**(math.floor(math.log10(v))-2)/2 for v in raw])/u
        y=raw/u
        A=np.array([p for n in x for p in ([1,n],[-1,-n])])
        rhs=np.array([p for v,hw in zip(y,h) for p in (v+hw,-v+hw)])
        vertices=[]
        for i,j in combinations(range(len(rhs)),2):
            if abs(np.linalg.det(A[[i,j]]))<1e-12:
                continue
            v=np.linalg.solve(A[[i,j]],rhs[[i,j]])
            if max(abs(v[0]+x*v[1]-y)/h)<=1+2e-7:
                vertices.append(v)
        assert vertices
        result[k]={"u":u,"u_lower":u-COMPARE[k]/2,
                   "joint_stage_upper":stage[k],"vertices":vertices}
    return result

def main():
    root=Path(__file__).resolve().parents[1]
    ap=argparse.ArgumentParser()
    ap.add_argument("--input-dir",type=Path,required=True)
    ap.add_argument("--output-dir",type=Path,required=True)
    args=ap.parse_args()
    coefficients=vertex_budget_coefficients(read_csv(args.input_dir/"tua_reported_grid.csv"),
                                           read_csv(args.input_dir/"tua_process_contributions.csv"))
    rows=[]
    examples={}
    N=100000  # integer enlargement of the source's fractional expected flows
    for B in range(2,126):
        ceilings={}
        for k,p in coefficients.items():
            upper=max(p["u"]*(v[0]+B*v[1]) for v in p["vertices"])
            ceilings[k]=1+(B*p["u_lower"]-upper)/((B-1)*p["joint_stage_upper"])
        limiting=min(ceilings,key=ceilings.get)
        for case,r in [("W1",F(55,10000)),("W2",F(46,10000))]:
            lo,hi=max(F(0),1-(B-1)*r),1-r
            retire=int(N*r)
            low=allocation(N,retire,B,True)
            high=allocation(N,retire,B,False)
            assert F(low["original_survivors_at_B"],N)==lo
            assert F(high["original_survivors_at_B"],N)==hi
            for field in ["initial_manufacturing","replacement_manufacturing","reconditioning_inputs",
                          "recycling_total_including_terminal_pool","delivered_services"]:
                assert low[field]==high[field]
            rows.append({"service_horizon":B,"rejection_case":case,
                         "rejection_fraction_exact":str(r),
                         "original_cohort_lower_exact":str(lo),"original_cohort_upper_exact":str(hi),
                         "original_cohort_lower":float(lo),"original_cohort_upper":float(hi),
                         "original_N":N,"replacement_count":low["replacement_manufacturing"],
                         "joint_multiplier_robust_W1_RE1":ceilings[limiting] if case=="W1" else None,
                         "limiting_W1_RE1_indicator":limiting if case=="W1" else None,
                         "interpretation":"synchronous-model identified set, not an observed asset cohort"})
            if case=="W1" and B in [3,8,10,20,125]:
                examples[str(B)]={"original_first":low,"replacement_first":high}
    frontiers=[]
    for alpha in [F(19,20),F(99,100)]:
        eligible=[r for r in rows if r["rejection_case"]=="W1" and
                  F(r["original_cohort_lower_exact"])>=alpha and r["joint_multiplier_robust_W1_RE1"]>1]
        best=max(eligible,key=lambda r:r["joint_multiplier_robust_W1_RE1"]) if eligible else None
        frontiers.append({"alpha_exact":str(alpha),"alpha":float(alpha),
                          "robustly_admissible_horizons":[r["service_horizon"] for r in eligible],
                          "maximum_headroom_horizon":None if best is None else best["service_horizon"],
                          "maximum_joint_multiplier_ceiling":None if best is None else best["joint_multiplier_robust_W1_RE1"],
                          "not_legal_requirement":True})
    args.output_dir.mkdir(parents=True,exist_ok=True)
    with (args.output_dir/"pool_budget_cohort_frontier_v20.csv").open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    out={"row_count":len(rows),"source":"Tua2019 Table1, sections2.3/2.7.1 and3; source vectors S13/S14",
         "assumptions":["N active positions serve synchronously each round", "r*N rejections at every inter-service reconditioning",
                        "immediate one-for-one replacement and guaranteed next dispatch", "no in-use non-return",
                        "source age-invariant process/recycling coefficients and stated replacement material inventory"],
         "bound":"max(0,1-(B-1)*r) <= original-cohort service attainment <= 1-r, B>=2",
         "evidence_status":"sharp mathematical model-completion bounds; not observed survival or changed LCA inventory",
         "coupled_frontiers":frontiers,
         "construction_examples":examples}
    (args.output_dir/"pool_budget_cohort_summary_v20.json").write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"row_count":len(rows),"coupled_frontiers":frontiers,
                      "selected_rows":[r for r in rows if r["rejection_case"]=="W1" and r["service_horizon"] in [3,5,8,10,20]]},indent=2))

if __name__=="__main__":
    main()
