"""Independent Project40 synthetic lifecycle audit; portable synthetic model verification.

State transitions independently generate paths and solve sequential job-to-asset
assignments. The original path_masks is called ONLY for differential comparison;
it is never used by the independent optimizer, generator, or witness simulator.
Run with Python -B. Outputs are confined to this fresh working copy.
"""

import ast
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction
from functools import lru_cache
import hashlib
from itertools import product
import json
from math import inf
from pathlib import Path
import sys
from time import perf_counter
from types import SimpleNamespace

import numpy as np
from scipy.optimize import linprog


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "code/run_lifecycle_v38.py"
PILOT = ROOT / "results/lifecycle_v38/pilot/pilot_results.json"
PROTOCOL = ROOT / "results/lifecycle_v38/pilot_protocol.md"
OUTPUT = ROOT / "results/lifecycle_v38/path_verification.json"
EXTENSION = ROOT / "results/lifecycle_v38/extension/structural_results.json"
EXTENSION_PROTOCOL = ROOT / "results/lifecycle_v38/extension_protocol.md"


def ceiling(value):
    value = Fraction(value)
    return -(-value.numerator // value.denominator)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def job_key(jobs):
    return tuple((j["id"], j["owner"], j["start"], j["ready"], j["lost"])
                 for j in jobs)


def old_key(assets):
    return tuple((a["id"], a["owner"], Fraction(a["age"]), a["previous_services"])
                 for a in assets)


def jobs_from_key(key):
    return [dict(zip(("id", "owner", "start", "ready", "lost"), item))
            for item in key]


def assets_from_key(key):
    return [dict(id=a, owner=o, age=age, previous_services=u) for a, o, age, u in key]


@lru_cache(maxsize=None)
def independent_paths(jkey, K, L, ownership, asset=None, endpoint="ready"):
    """Enumerate all ordered, non-repeating trajectories by actual state steps.

    No sorted-subset extension: every remaining job is tried at each state.
    An initial asset is explicitly identified by asset, even for pooled ownership.
    For the alternative start endpoint, completion=start is a SYNTHETIC proxy.
    """
    jobs = jobs_from_key(jkey)
    paths = set()
    if asset is None:
        birth, free, remaining, owner = None, None, K, None
    else:
        _, original_owner, age, prior = asset
        birth, free, remaining, owner = -age, 0, K - prior, original_owner

    def visit(sequence, born, available, quota, required_owner, absorbed):
        if sequence:
            paths.add(sequence)
        if absorbed or quota <= 0:
            return
        used = set(sequence)
        for i, job in enumerate(jobs):
            if i in used or (available is not None and job["start"] < available):
                continue
            if ownership == "dedicated" and required_owner is not None and job["owner"] != required_owner:
                continue
            actual_birth = job["start"] if born is None else born
            completion = job[endpoint]
            if completion > actual_birth + L:
                continue
            next_owner = job["owner"] if required_owner is None else required_owner
            visit(sequence + (i,), actual_birth, job["ready"], quota - 1,
                  next_owner, bool(job["lost"]))

    visit((), birth, free, remaining, owner, False)
    return frozenset(paths)


@lru_cache(maxsize=None)
def assignment_optimum(jkey, K, L, ownership, initial=(), endpoint="ready"):
    """Independent exact state assignment, not a path-cover DP or MILP.

    A state is (absolute expiry, next available time, remaining uses, owner).
    Each extant asset is updated once when it serves the next calendar job.
    Dead/exhausted assets cannot re-enter. New births occur at first dispatch.
    """
    jobs = sorted(jobs_from_key(jkey), key=lambda j: (j["start"], j["id"]))
    first_states = tuple((L - age, 0, K - prior,
                          owner if ownership == "dedicated" else -1)
                         for _, owner, age, prior in initial if K > prior)

    def normalize(position, states):
        if position == len(jobs):
            return ()
        now = jobs[position]["start"]
        earliest_end = min(j[endpoint] for j in jobs[position:])
        return tuple(sorted((expiry, available if available > now else 0, quota, owner)
                            for expiry, available, quota, owner in states
                            if quota > 0 and expiry >= earliest_end))

    @lru_cache(maxsize=None)
    def search(position, states):
        if position == len(jobs):
            return 0
        job = jobs[position]
        answer = inf
        for state in set(states):  # identical resources are interchangeable
            expiry, available, quota, owner = state
            if (available > job["start"] or expiry < job[endpoint]
                    or (ownership == "dedicated" and owner != job["owner"])):
                continue
            later = list(states)
            later.remove(state)
            if quota > 1 and not job["lost"]:
                later.append((expiry, job["ready"], quota - 1, owner))
            answer = min(answer, search(position + 1, normalize(position + 1, later)))
        if K > 0 and job[endpoint] <= job["start"] + L:
            later = list(states)
            if K > 1 and not job["lost"]:
                later.append((job["start"] + L, job["ready"], K - 1,
                              job["owner"] if ownership == "dedicated" else -1))
            answer = min(answer, 1 + search(position + 1, normalize(position + 1, later)))
        return answer

    answer = search(0, normalize(0, first_states))
    return (None if answer == inf else int(answer), search.cache_info().currsize)


def simulate_witness(row, selected=None, ownership=None):
    """Check actual listed asset trajectories directly, independent of any library."""
    jobs = {j["id"]: j for j in row["jobs"]}
    initial = {a["id"]: a for a in row["initial_assets"]}
    ownership = ownership or row["ownership"]
    selected = row["selected_paths"] if selected is None else selected
    cover = Counter()
    used_old = set()
    traces = []
    for number, path in enumerate(selected):
        ids = path["jobs"]
        assert ids and len(ids) == len(set(ids))
        old_id = path["old"]
        if old_id is None:
            birth = jobs[ids[0]]["start"]
            free, prior = birth, 0
            owner = jobs[ids[0]]["owner"]
        else:
            assert old_id in initial and old_id not in used_old
            used_old.add(old_id)
            asset = initial[old_id]
            birth, free, prior, owner = -asset["age"], 0, asset["previous_services"], asset["owner"]
        absorbed = False
        events = []
        for step, job_id in enumerate(ids):
            job = jobs[job_id]
            assert not absorbed, (row["name"], "service after loss", ids)
            assert job["start"] >= free
            assert job["ready"] <= birth + row["L"]
            assert prior + step + 1 <= row["K"]
            if ownership == "dedicated":
                assert job["owner"] == owner
            events.append(dict(job=job_id, waiting_before_dispatch=job["start"] - free,
                               calendar_age_at_dispatch=job["start"] - birth,
                               calendar_age_at_ready=job["ready"] - birth,
                               cumulative_services=prior + step + 1,
                               lost_after=job["lost"]))
            free, absorbed = job["ready"], bool(job["lost"])
            cover[job_id] += 1
        traces.append(dict(path_number=number, old=old_id, birth=birth,
                           expiry=birth + row["L"], jobs=ids, events=events))
    assert cover == Counter({j: 1 for j in jobs}), row["name"]
    births = sum(p["old"] is None for p in selected)
    assert births == row["introductions"]
    return traces


def primal_gap_witness(jobs, new, olds, claimed_value):
    """Exact rational primal proof that a gap is an LP/integer gap, not bad duals."""
    columns = [(p, None) for p in sorted(new)]
    columns += [(p, old) for old, family in enumerate(olds) for p in sorted(family)]
    coverage = np.array([[int(j in p) for p, _ in columns] for j in range(len(jobs))], dtype=float)
    quotas = np.array([[int(old == a) for _, old in columns] for a in range(len(olds))], dtype=float)
    cost = np.array([int(old is None) for _, old in columns], dtype=float)
    solved = linprog(cost, A_eq=coverage, b_eq=np.ones(len(jobs)),
                     A_ub=quotas if olds else None, b_ub=np.ones(len(olds)) if olds else None,
                     bounds=(0, None), method="highs")
    assert solved.success
    x = [Fraction(float(v)).limit_denominator(100000) for v in solved.x]
    assert all(v >= 0 for v in x)
    assert all(sum(v for v, (p, _) in zip(x, columns) if j in p) == 1 for j in range(len(jobs)))
    assert all(sum(v for v, (_, old) in zip(x, columns) if old == a) <= 1 for a in range(len(olds)))
    value = sum(v for v, (_, old) in zip(x, columns) if old is None)
    assert value == claimed_value, (value, claimed_value)
    return [dict(jobs=[jobs[i]["id"] for i in p], old=old, multiplicity=str(v))
            for v, (p, old) in zip(x, columns) if v]


def verify_certificate(row, new, olds):
    stored = row["certificate"]
    weights = [Fraction(w) for w in stored["weights"]]
    assert len(weights) == len(row["jobs"]) and all(w >= 0 for w in weights)
    weight = lambda p: sum((weights[i] for i in p), Fraction(0))
    assert all(weight(p) <= 1 for p in new)
    allowances = [max((weight(p) for p in family), default=Fraction(0)) for family in olds]
    assert [str(a) for a in allowances] == stored["old_allowances"]
    value = sum(weights) - sum(allowances)
    assert value == Fraction(stored["value_fraction"])
    assert max(0, ceiling(value)) == stored["ceil"]
    assert abs(float(value) - stored["value"]) < 1e-10
    return weights, allowances, value


def semantic_probes(source_function):
    def source_paths(jobs, K, L, ownership, old_owner=None, age=0, prior=0):
        return {tuple(p["jobs"]) for p in source_function(
            [SimpleNamespace(**j) for j in jobs], K, L, ownership, old_owner, age, prior)}

    loss = [dict(id=0, owner=0, start=0, ready=2, lost=True),
            dict(id=1, owner=0, start=2, ready=4, lost=False)]
    assert (0, 1) not in source_paths(loss, 3, 8, "pooled")
    assert (0, 1) not in independent_paths(job_key(loss), 3, 8, "pooled")
    waits = [dict(id=0, owner=0, start=4, ready=6, lost=False)]
    old = (0, 0, Fraction(3), 0)
    assert not independent_paths(job_key(waits), 3, 8, "pooled", old)
    assert not source_paths(waits, 3, 8, "pooled", 0, 3)
    assert independent_paths(job_key(waits), 3, 8, "pooled") == {(0,)}
    # Same helper called with owner=None silently chooses fresh-birth semantics.
    assert source_paths(waits, 3, 8, "pooled", None, 3) == {(0,)}
    final = [dict(id=0, owner=0, start=0, ready=4, lost=False)]
    assert not independent_paths(job_key(final), 1, 3, "pooled")
    assert independent_paths(job_key(final), 1, 3, "pooled", endpoint="start") == {(0,)}
    boundary = [dict(id=0, owner=0, start=0, ready=3, lost=False)]
    assert independent_paths(job_key(boundary), 1, 3, "pooled") == {(0,)}
    return dict(loss_is_absorbing=True, initial_waiting_consumes_age=True,
                inclusive_ready_expiry_boundary=True,
                final_ready_rule_stronger_than_instantaneous_delivery=True,
                owner_none_helper_probe=dict(
                    old_age=3, L=8, job_start=4, job_ready=6,
                    explicit_old_paths=[], original_helper_old_owner_none_paths=[[0]],
                    affects_any_of_578=False,
                    implication="old_owner is an origin sentinel in the helper; ownerless pooled old assets need an explicit origin flag"))


def scan_age_pairs(rows):
    index = {r["name"]: r for r in rows}
    pilot_contrasts, age_profiles = [], []
    for dedicated in rows:
        if dedicated["ownership"] != "dedicated" or dedicated["state"] != "young":
            continue
        base = dedicated["name"].rsplit("_young_", 1)[0]
        pooled = index[base + "_young_pooled"]
        old_d, old_p = index[base + "_old_dedicated"], index[base + "_old_pooled"]
        initial_gap = dedicated["introductions"] - pooled["introductions"]
        old_gap = old_d["introductions"] - old_p["introductions"]
        if initial_gap > 0 and old_gap == 0:
            pilot_contrasts.append(dict(base=base,
                young_names=[dedicated["name"], pooled["name"]],
                young_births=[dedicated["introductions"], pooled["introductions"]],
                old_names=[old_d["name"], old_p["name"]],
                old_births=[old_d["introductions"], old_p["introductions"]],
                old_age=dedicated["L"] - 3, old_prior=dedicated["K"] - 1,
                confounded_age_and_prior=True))
        if initial_gap <= 0:
            continue
        profile = []
        for age in range(dedicated["L"] + 1):
            initial = tuple((o, o, Fraction(age), 0) for o in range(3))
            d = assignment_optimum(job_key(dedicated["jobs"]), dedicated["K"], dedicated["L"], "dedicated", initial)[0]
            p = assignment_optimum(job_key(dedicated["jobs"]), dedicated["K"], dedicated["L"], "pooled", initial)[0]
            assert d >= p
            profile.append(dict(age=age, dedicated=d, pooled=p, gain=d-p))
        erased = [x["age"] for x in profile if x["gain"] == 0]
        prior_only = []
        for u in range(dedicated["K"]):
            initial = tuple((o, o, Fraction(0), u) for o in range(3))
            d = assignment_optimum(job_key(dedicated["jobs"]), dedicated["K"], dedicated["L"], "dedicated", initial)[0]
            p = assignment_optimum(job_key(dedicated["jobs"]), dedicated["K"], dedicated["L"], "pooled", initial)[0]
            prior_only.append(dict(prior=u, dedicated=d, pooled=p, gain=d-p))
        age_profiles.append(dict(base=base, fixed_prior=0, integer_age_profile=profile,
                                 first_integer_erasure=min(erased) if erased else None,
                                 prior_only_at_age0=prior_only,
                                 age_gain_is_nonincreasing=all(b["gain"] <= a["gain"]
                                     for a,b in zip(profile, profile[1:]))))
    erased_profiles = [p for p in age_profiles if p["first_integer_erasure"] is not None]
    minimal_age = min(p["first_integer_erasure"] for p in erased_profiles)
    minimal_profiles = [p for p in erased_profiles if p["first_integer_erasure"] == minimal_age]
    preferred = next((p for p in minimal_profiles if p["base"] == "synchronous_delay2_K2_L8_loss0"),
                     minimal_profiles[0])
    reference = index[preferred["base"] + "_young_dedicated"]
    # The family changes at a=L-ready_j. Equality remains feasible, so real-age
    # erasure may have an infimum but no smallest positive age.
    critical = sorted({Fraction(0), Fraction(reference["L"])} |
                      {Fraction(reference["L"] - j["ready"]) for j in reference["jobs"]
                       if 0 <= reference["L"] - j["ready"] <= reference["L"]})
    cells = []
    for i, left in enumerate(critical):
        initial = tuple((o, o, left, 0) for o in range(3))
        d = assignment_optimum(job_key(reference["jobs"]), reference["K"], reference["L"], "dedicated", initial)[0]
        p = assignment_optimum(job_key(reference["jobs"]), reference["K"], reference["L"], "pooled", initial)[0]
        cells.append(dict(domain=f"a={left}", lower=str(left), upper=str(left), point=True,
                          dedicated=d, pooled=p, gain=d-p))
        if i+1 < len(critical):
            right = critical[i+1]
            sample = (left+right)/2
            initial = tuple((o, o, sample, 0) for o in range(3))
            d = assignment_optimum(job_key(reference["jobs"]), reference["K"], reference["L"], "dedicated", initial)[0]
            p = assignment_optimum(job_key(reference["jobs"]), reference["K"], reference["L"], "pooled", initial)[0]
            cells.append(dict(domain=f"{left}<a<{right}", lower=str(left), upper=str(right), point=False,
                              dedicated=d, pooled=p, gain=d-p))
    first_zero = next(cell for cell in cells if cell["gain"] == 0)
    minimal_pair = dict(base=preferred["base"], K=reference["K"], L=reference["L"],
        delay=reference["delay"], loss=reference["loss"], common_prior=0,
        old_stock_count=3, old_assets_initially_ready=True,
        reference_names=[preferred["base"] + "_age0_prior0_dedicated", preferred["base"] + "_age0_prior0_pooled"],
        aged_names=[preferred["base"] + f"_age{minimal_age}_prior0_dedicated",
                    preferred["base"] + f"_age{minimal_age}_prior0_pooled"],
        reference_births=[preferred["integer_age_profile"][0]["dedicated"], preferred["integer_age_profile"][0]["pooled"]],
        aged_births=[preferred["integer_age_profile"][minimal_age]["dedicated"], preferred["integer_age_profile"][minimal_age]["pooled"]],
        minimum_tested_integer_age=minimal_age, real_age_infimum=first_zero["lower"],
        real_infimum_attained=first_zero["point"], piecewise_real_age_cells=cells,
        scope="minimum over the pilot's initially advantageous young contexts, uniform integer initial age, holding prior=0 fixed")
    return dict(pilot_young_old_erasure_pairs=pilot_contrasts,
                smallest_pilot_old_age_among_erasure_pairs=min(p["old_age"] for p in pilot_contrasts),
                age_only_profiles=age_profiles, named_minimum_integer_age_pair=minimal_pair)


def regular_structure_checks():
    checked = 0
    examples = []
    for N, h, K, spacing, delay, life_steps in product(range(1,7),range(1,4),range(1,5),
                                                    (1,2),(1,2),(0,1,2,4)):
        if delay > spacing:
            continue
        L = delay + life_steps * spacing
        jobs = [dict(id=j, owner=j%h, start=j*spacing, ready=j*spacing+delay, lost=False)
                for j in range(N)]
        mp = min(K, 1 + (L-delay)//spacing)
        md = min(K, 1 + (L-delay)//(h*spacing))
        counts = Counter(j["owner"] for j in jobs)
        expected_p = ceiling(Fraction(N, mp))
        expected_d = sum(ceiling(Fraction(n, md)) for n in counts.values())
        actual_p = assignment_optimum(job_key(jobs),K,L,"pooled")[0]
        actual_d = assignment_optimum(job_key(jobs),K,L,"dedicated")[0]
        assert (actual_p,actual_d)==(expected_p,expected_d)
        q = len(counts)
        lower_gain = ceiling(Fraction(N, md))-ceiling(Fraction(N, mp))
        assert lower_gain <= actual_d-actual_p <= lower_gain+q-1
        checked += 1
        if (N,h,K,spacing,delay,life_steps)==(6,3,2,1,1,4):
            examples.append(dict(N=N,h=h,K=K,delta=spacing,tau=delay,L=L,
                                 mP=mp,mD=md,BP=actual_p,BD=actual_d))
    def run(times, K, L, tau=1):
        jobs=[dict(id=i,owner=0,start=t,ready=t+tau,lost=False) for i,t in enumerate(times)]
        return assignment_optimum(job_key(jobs),K,L,"pooled")[0]
    isolated=run((0,1,2,10,11,12),2,3)
    assert isolated==4
    boundary=run((0,2),2,3)
    assert boundary==1
    irregular=run((0,1,2,4,6,8,12),4,3)
    assert irregular==4
    assert run((0,1,2,4,6,8),4,3)==3
    return dict(equal_spacing_cases_checked=checked, both_modes_exact=True,
        regular_examples=examples,
        isolated_regular_clusters=dict(times=[0,1,2,10,11,12],tau=1,L=3,K=2,
            cluster_births=[2,2],actual_total=4,global_ceil_bound=3),
        gap_equality_counterexample=dict(times=[0,2],tau=1,L=3,K=2,gap=2,L_minus_tau=2,
                                        actual_births=1,strict_greater_than_required=True),
        arbitrary_cluster_count_counterexample=dict(times=[0,1,2,4,6,8,12],tau=1,L=3,K=4,
            clusters=[[0,1,2,4,6,8],[12]],max_path_sizes=[3,1],naive_sum_of_ceils=3,
            actual_cluster_births=[3,1],actual_total=4,
            obstruction_subset=[0,4,8],message="add exact cluster optima; ceil(n/max_path_size) needs regular-block assumptions"))


def identification_examples():
    jobs=[dict(id=3*c+o,owner=o,start=3*c,ready=3*c+1,lost=False) for c in range(3) for o in range(3)]
    K,L=3,12
    ranges=[dict(id=0,owner=0,age=[3,3],prior=[1,1],observation="known birth and prior use count"),
            dict(id=1,owner=1,age=[4,7],prior=[0,1],observation="interval known age and partial use log"),
            dict(id=2,owner=2,age=[6,10],prior=[1,2],
                 observation="lower-censored initial age; upper=10 is a declared external bound, not inferred from censoring")]
    scenarios=[]
    for a1,a2,u1,u2,tau in product(range(4,8),range(6,11),range(2),(1,2),(1,2)):
        actual_jobs=[dict(j,ready=j["start"]+tau) for j in jobs]
        initial=((0,0,Fraction(3),1),(1,1,Fraction(a1),u1),(2,2,Fraction(a2),u2))
        assert all(age>=u*tau for _,_,age,u in initial)  # this example supports stationary prior delay
        births={mode:assignment_optimum(job_key(actual_jobs),K,L,mode,initial)[0]
                for mode in ("dedicated","pooled")}
        scenarios.append(dict(ages=[3,a1,a2],prior=[1,u1,u2],tau=tau,**births))
    optimistic=next(s for s in scenarios if s["ages"]==[3,4,6] and s["prior"]==[1,0,1] and s["tau"]==1)
    pessimistic=next(s for s in scenarios if s["ages"]==[3,7,10] and s["prior"]==[1,1,2] and s["tau"]==2)
    extrema={mode:[min(s[mode] for s in scenarios),max(s[mode] for s in scenarios)]
             for mode in ("dedicated","pooled")}
    assert all(extrema[mode]==[optimistic[mode],pessimistic[mode]] for mode in extrema)
    # A sharp interval hull need not contain every achievable integer value.
    simple=[dict(id=o,owner=o,start=0,ready=1,lost=False) for o in range(3)]
    initial=tuple((o,o,Fraction(2),0) for o in range(3))
    discrete=[]
    for tau in (Fraction(1),Fraction(3,2),Fraction(2)):
        actual=[dict(j,ready=tau) for j in simple]
        B=assignment_optimum(job_key(actual),1,3,"pooled",initial)[0]
        discrete.append(dict(tau=str(tau),births=B))
    assert {s["births"] for s in discrete}=={0,3}
    return dict(declared_synthetic=True,K=K,L=L,initial_observation_bounds=ranges,tau_box=[1,2],
        finite_grid_states=len(scenarios),optimistic_corner=optimistic,pessimistic_corner=pessimistic,
        sharp_interval_hulls=extrema,
        realizable_integer_values={mode:sorted({s[mode] for s in scenarios}) for mode in extrema},
        interval_hull_not_full_identified_set=dict(K=1,L=3,old_ages=[2,2,2],prior=[0,0,0],
            simultaneous_job_count=3,tau_box=[1,2],checks=discrete,
            exact_identified_set=[0,3],interval_hull=[0,3],unattainable_integers=[1,2]))


def propositions():
    return [
        dict(name="Equal-spacing exact birth count",status="proved under supplied restrictions",
            assumptions="N>=1, integer K>=1,h>=1,delta>0,0<tau<=delta,L>=tau; new birth at first dispatch; ready endpoint; no old/loss/coupled capacities",
            statement="mP=min(K,1+floor((L-tau)/delta)); BP=ceil(N/mP). For cyclic h owners, mD=min(K,1+floor((L-tau)/(h*delta))); BD=sum_o ceil(n_o/mD).",
            proof="q ordered shared jobs span at least (q-1)*delta; final ready needs that span+tau<=L, and q<=K. Thus at most mP per asset. Consecutive blocks of mP meet both expiry and readiness because tau<=delta, and attain ceil(N/mP). Dedicated paths have successive owner times at least h*delta; apply the same argument separately to each owner and construct consecutive owner blocks. No coupling prevents simultaneous execution.",
            gain="Let q=#nonempty owners. Delta=BD-BP>=0 and ceil(N/mD)-ceil(N/mP)<=Delta<=ceil(N/mD)-ceil(N/mP)+q-1.",
            thresholds="mP>=q iff K>=q and L>=tau+(q-1)*delta; replace delta by h*delta for mD. For a birth ceiling b>=1, BP<=b iff K>=ceil(N/b) and L>=tau+(ceil(N/b)-1)*delta.",
            limitation="Even mP=mD=K can leave an owner-rounding gain. This is a model-level structural result, not a novelty or physical-lifetime claim."),
        dict(name="Isolated cluster additivity",status="proved with qualification",
            statement="If ordered cluster separation min_start(next)-max_start(previous)>L-tau, no new path crosses clusters. Therefore B*=sum_b B_b*. If each cluster satisfies equal-spacing assumptions, B_b*=ceil(n_b/m_b).",
            proof="A path born at its first job in the earlier cluster would reach the later cluster's ready endpoint more than L after birth: gap+tau>L. Hence feasible columns partition by cluster. Restrict any cover to each cluster for the lower bound, and concatenate independent optimal covers for equality.",
            limitation="Strict > is sufficient, not necessary; equality may allow crossing. Isolation alone does not make ceil(n_b/global maximum path length) exact for an irregular cluster; an explicit counterexample is recorded."),
        dict(name="Sharp monotone box bounds",status="proved conditionally",
            statement="Fix starts, owners,K,L,initial identities and the same deterministic loss/rules. For admissible closed box coordinates a,u,tau, P(worst) subset P(theta) subset P(best). Thus B*(best)<=B*(theta)<=B*(worst). If both simultaneous corners belong to the allowed parameter set, these are attained sharp extrema (possibly infeasible/infinite).",
            proof="Increasing a lowers absolute old expiry; increasing integer prior u lowers its residual use count; increasing tau raises ready endpoints and next-availability times. Each move can only remove paths. A cover feasible in the smallest family remains feasible in every larger family. Corner membership supplies both attaining parameter configurations.",
            common_policy="For a pure birth ceiling with fixed starts and no extra nonmonotone rows, a worst-corner assignment is one assignment feasible throughout this nested box. Hence worst-case adaptive feasibility and a common robust assignment coincide here; this implication does not extend automatically to correlated nonrectangular states or additional budgets.",
            identification="Use I={B*(theta):theta in C}. [min I,max I] is its sharp interval hull, not necessarily every achievable integer. A censoring lower bound does not supply an upper age bound or a future-life guarantee. Marginal corners excluded by history/selection constraints give only conservative outer bounds; optimize over the joint admissible C for sharp bounds.",
            limitation="Do not transfer box extrema to the shared gain difference: two separately monotone optimal birth counts need not have a monotone difference. No coverage for unknown losses, undefined prehistory, or real probabilities is claimed."),
        dict(name="Endpoint bracketing",status="proved conditionally",
            statement="With actual delivery completions e_j in [start_j,ready_j], fixed availability/loss rules and all else equal, B_start<=B_actual_delivery<=B_ready. Ready-endpoint certificates apply to the protocol's stronger service boundary, not automatically to the weaker delivery-only boundary.",
            proof="Every ready-expiry-feasible trajectory meets all delivery endpoints, and every delivery-expiry-feasible trajectory meets dispatch endpoints. Family inclusion gives the inequalities.",
            limitation="Start-endpoint tests explicitly assume instantaneous synthetic delivery; the pilot does not provide actual e_j. A lost terminal job's ready label is its declared end/unavailability deadline, not proof of successful return."),
    ]


def main():
    if OUTPUT.exists():
        raise FileExistsError(f"Refusing to overwrite another edit: {OUTPUT}")
    tick=perf_counter()
    before={p.relative_to(ROOT).as_posix():digest(p) for p in (SOURCE,PILOT,PROTOCOL)}
    source_text=SOURCE.read_text(encoding="utf-8")
    tree=ast.parse(source_text)
    original=next(node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name=="path_masks")
    namespace={}
    exec(compile(ast.Module(body=[original],type_ignores=[]),str(SOURCE),"exec"),namespace)
    original_function=namespace["path_masks"]
    doc=json.loads(PILOT.read_text(encoding="utf-8"))
    rows=doc["results"]
    assert len(rows)==578 and len({r["name"] for r in rows})==578
    assert doc["declared_synthetic_scenarios"] and all(r["feasible"] for r in rows)
    assert doc["protocol_sha256"]==before[PROTOCOL.relative_to(ROOT).as_posix()]
    verified=[]; gaps=[]; changed_endpoint=[]; history_incompatible=[]
    seen_paths={}
    for number,row in enumerate(rows,1):
        jobs=row["jobs"]; key=job_key(jobs); old=old_key(row["initial_assets"])
        new=independent_paths(key,row["K"],row["L"],row["ownership"])
        olds=[independent_paths(key,row["K"],row["L"],row["ownership"],a) for a in old]
        observed=[SimpleNamespace(**j) for j in jobs]
        assert new=={tuple(p["jobs"]) for p in original_function(observed,row["K"],row["L"],row["ownership"])}
        for a,family in zip(old,olds):
            _,owner,age,prior=a
            assert family=={tuple(p["jobs"]) for p in original_function(observed,row["K"],row["L"],row["ownership"],owner,age,prior)}
        assert len(new)+sum(map(len,olds))==row["path_count"]
        traces=simulate_witness(row)
        optimum,states=assignment_optimum(key,row["K"],row["L"],row["ownership"],old)
        assert optimum==row["introductions"]==row["dp_optimum"]
        weights,allowances,value=verify_certificate(row,new,olds)
        relaxed,_=assignment_optimum(key,row["K"],row["L"],row["ownership"],old,"start")
        assert relaxed<=optimum
        peak=max(sum(j["start"]<=t<j["ready"] for j in jobs)
                 for t in sorted({j["start"] for j in jobs}|{j["ready"] for j in jobs}))
        assert peak==row["occupied_peak"]
        assert sum(len(p["jobs"])-1 for p in row["selected_paths"])==row["ready_transitions"]
        assert sum(j["lost"] for j in jobs)==row["loss_terminal_jobs"]
        assert row["certificate_gap"]==optimum-row["certificate"]["ceil"]
        assert max(0,ceiling(Fraction(len(jobs)-sum(row["K"]-a[3] for a in old),row["K"])))==row["simple_global_use_bound"]
        if relaxed!=optimum:
            changed_endpoint.append(dict(name=row["name"],ready_births=optimum,
                                         instantaneous_delivery_births=relaxed))
        if "delay" in row and any(age<prior*row["delay"] for _,_,age,prior in old):
            history_incompatible.append(row["name"])
        entry=dict(name=row["name"],independent_state_optimum=optimum,states_explored=states,
                   independent_new_path_count=len(new),independent_old_path_counts=list(map(len,olds)),
                   witness_valid=True,exact_certificate_valid=True,certificate_ceiling=row["certificate"]["ceil"],
                   instantaneous_delivery_optimum=relaxed,selected_asset_traces=traces,
                   unused_old_ids=[a[0] for a in old if a[0] not in {p["old"] for p in row["selected_paths"]}])
        verified.append(entry)
        seen_paths[row["name"]]=(new,olds)
        if row["certificate_gap"]>0:
            assert row["ownership"]=="dedicated"
            owner_components=[]
            for owner in sorted({j["owner"] for j in jobs}):
                weight_total=sum(w for j,w in zip(jobs,weights) if j["owner"]==owner)
                old_allowance=sum(v for a,v in zip(row["initial_assets"],allowances) if a["owner"]==owner)
                component=weight_total-old_allowance
                owner_components.append(dict(owner=owner,dual_component=str(component),
                                              rounded_bound=max(0,ceiling(component))))
            rounded=sum(p["rounded_bound"] for p in owner_components)
            assert rounded==optimum
            primal=primal_gap_witness(jobs,new,olds,value)
            gaps.append(dict(name=row["name"],integer_optimum=optimum,
                exact_LP_optimum=str(value),global_certificate_ceiling=row["certificate"]["ceil"],
                gap=row["certificate_gap"],owner_components=owner_components,
                ownerwise_rounded_certificate=rounded,integer_constructive_paths=row["selected_paths"],
                exact_fractional_cover=primal,
                false_pass_budget_example=dict(environmental_birth_allowance=row["certificate"]["ceil"]+0.5,
                                               cash_birth_allowance=row["certificate"]["ceil"]+0.5,
                                               necessary_passes=True,integer_joint_feasible=False,
                                               normalized_asset_equivalents_only=True)))
        if number%100==0:
            print(f"SYNTHETIC independently checked {number}/578 rows",flush=True)
    assert len(gaps)==14
    pairs=0
    index={r["name"]:r for r in rows}
    for row in rows:
        if row["ownership"]!="dedicated": continue
        shared=index[row["name"].rsplit("_",1)[0]+"_pooled"]
        dn,do=seen_paths[row["name"]];pn,po=seen_paths[shared["name"]]
        assert dn<=pn and all(d<=p for d,p in zip(do,po))
        assert shared["introductions"]<=row["introductions"]
        simulate_witness(row,ownership="pooled")  # its birth count need not equal the shared minimum
        pairs+=1
    print("SYNTHETIC all 578 trajectories/optima and 14 exact LP gaps verified; scanning age-only pairs",flush=True)
    age=scan_age_pairs(rows)
    print("SYNTHETIC age-only pairs verified; checking structural formulas and partial-age bounds",flush=True)
    structures=regular_structure_checks()
    identification=identification_examples()
    after={p.relative_to(ROOT).as_posix():digest(p) for p in (SOURCE,PILOT,PROTOCOL)}
    assert before==after,"Input changed concurrently; no output written"
    result=dict(status="PASS_FOR_DECLARED_SYNTHETIC_READY_ENDPOINT_MODEL",
        completed_at_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=perf_counter()-tick,
        declared_synthetic=True,current_EU_physical_validation=False,
        inputs_sha256=before,input_hashes_unchanged=True,
        methodology=dict(generator="unsorted all-remaining-job state/permutation enumeration",
            optimizer="direct sequential job-to-extant-asset assignment with explicit expiry/remaining-use state",
            witness_check="per-asset event simulator over original selected paths, not membership alone",
            original_path_library_use="differential comparison only; not used for lifecycle feasibility or optimum",
            certificate_check="exact Fraction inequalities on independently generated paths; 14 gaps also have exact fractional primal witnesses"),
        summary=dict(rows_checked=len(rows),feasible_witnesses_verified=len(verified),
            independent_minimum_births_matches=578,path_family_mismatches=0,invalid_selected_paths=0,
            ownership_inclusion_pairs=pairs,weighted_certificate_gaps=len(gaps),
            all_gaps_closed_by_ownerwise_rounding=True,
            cases_changed_under_instantaneous_delivery_endpoint=len(changed_endpoint),
            conditional_stationary_prehistory_incompatibilities=len(history_incompatible),
            regular_formula_cases=structures["equal_spacing_cases_checked"]),
        semantic_probes=semantic_probes(original_function),
        endpoint_review=dict(protocol_ready_end_consistent=True,
            stronger_than_delivery_only=True,start_endpoint_is_declared_instantaneous_proxy=True,
            no_actual_delivery_completion_times_provided=True,
            changed_cases=changed_endpoint),
        old_state_review=dict(initial_assets_are_available_at_time0=True,
            age_continues_during_waiting=True,unused_old_resources_not_new_births=True,
            occupied_peak_is_job_occupancy_not_total_fleet=True,
            history_condition="Only if every previous service used the same delay tau, birth was first prior dispatch, and all prior services reached ready before time0: age>=prior_services*tau is necessary.",
            incompatible_case_names_under_that_extra_condition=history_incompatible,
            statement="These do not invalidate forward-state pilot arithmetic: the protocol does not identify stationary prehistory. They forbid treating all initial states as verified histories."),
        necessary_certificate_gaps=gaps,age_pair_analysis=age,
        structural_results=structures,identification_design=identification,
        conservative_mathematical_propositions=propositions(),verified_cases=verified,
        limitations=["All demand, delay, lifetime, prior states and loss scenarios are synthetic.",
            "DP/MILP sharing a path library is not independent lifecycle validation; this audit generates and simulates states separately.",
            "No contemporary LCI/cash map or physical service history was audited.",
            "Age-only integer minima are grid/scope-specific; real-age infima need not be attained.",
            "Rectangular uncertainty endpoints must satisfy joint observation/history restrictions before sharpness can be claimed.",
            "Known job losses are absorbing scenario events, not predictions, probability coverage or online recourse."])
    with OUTPUT.open("x",encoding="utf-8") as handle:
        json.dump(result,handle,ensure_ascii=False,indent=2,default=lambda v:str(v) if isinstance(v,Fraction) else TypeError(type(v)))
        handle.write("\n")
    print(json.dumps(dict(summary=result["summary"],named_age_pair=age["named_minimum_integer_age_pair"],
                          identification_hulls=identification["sharp_interval_hulls"],output=str(OUTPUT)),
                     ensure_ascii=False,indent=2),flush=True)


def canonical_review():
    """Append a focused current-TeX audit using independent state trajectories.

    This reads but never imports extension or original lifecycle modules. The
    previously completed 578 audit is preserved, not replayed through its paths.
    """
    tick = perf_counter()
    output_hash = digest(OUTPUT)
    result = json.loads(OUTPUT.read_text(encoding="utf-8"))
    inputs = (EXTENSION, EXTENSION_PROTOCOL, SOURCE, PILOT, PROTOCOL)
    before = {p.relative_to(ROOT).as_posix(): digest(p) for p in inputs}
    extension = json.loads(EXTENSION.read_text(encoding="utf-8"))
    assert extension["protocol_sha256"] == before[EXTENSION_PROTOCOL.relative_to(ROOT).as_posix()]
    assert all(before[p] == value for p, value in result["inputs_sha256"].items())
    formulas = []
    expected_keys = set()
    for N,h,K,L,delta,tau,mode in product((3,6,9),(1,2,3),(2,3,5),
                                         (2,4,8),(1,2),(1,2),("pooled","dedicated")):
        if tau <= delta and L >= tau:
            expected_keys.add((N,h,K,L,delta,tau,mode))
    observed_keys = set()
    for r in extension["formula_rows"]:
        N,h,K,L,delta,tau,mode = (r[k] for k in ("N","owners","K","L","spacing","delay","ownership"))
        parameter_key = (N,h,K,L,delta,tau,mode)
        assert parameter_key not in observed_keys
        observed_keys.add(parameter_key)
        jobs = [dict(id=j, owner=j%h, start=j*delta, ready=j*delta+tau, lost=False) for j in range(N)]
        m = min(K, 1+(L-tau)//(delta if mode == "pooled" else h*delta))
        counts = Counter(j["owner"] for j in jobs)
        formula = ceiling(Fraction(N,m)) if mode == "pooled" else sum(ceiling(Fraction(n,m)) for n in counts.values())
        value, states = assignment_optimum(job_key(jobs), K, L, mode)
        assert value == formula == r["formula"] == r["exact_births"] and m == r["m"]
        groups = [list(range(N))] if mode == "pooled" else [[j for j in range(N) if j%h == o] for o in range(h)]
        paths = [dict(old=None,jobs=g[i:i+m]) for g in groups for i in range(0,len(g),m)]
        witness = dict(name=f"canonical_formula_{parameter_key}", jobs=jobs, K=K, L=L,
                       ownership=mode, initial_assets=[], selected_paths=paths, introductions=value)
        simulate_witness(witness)
        formulas.append(dict(**r, independent_state_births=value, states_explored=states,
                             independently_simulated_block_paths=paths))
    assert observed_keys == expected_keys and len(formulas) == extension["formula_checks"] == 486
    print("SYNTHETIC canonical 486 formula values and consecutive-block witnesses independently verified", flush=True)

    boxes = []
    history_exclusions = 0
    for box in extension["state_boxes"]:
        mode = box["ownership"]
        values = []
        checked_keys = set()
        def jobs_at(tau):
            return [dict(id=i,owner=0 if i<2 else 1,start=t,ready=t+tau,lost=False)
                    for i,t in enumerate((0,1,4,5))]
        best_jobs, worst_jobs = jobs_at(1), jobs_at(2)
        best_assets = ((0,0,Fraction(0),0),(1,1,Fraction(0),0))
        worst_assets = ((0,0,Fraction(4),2),(1,1,Fraction(4),2))
        best_new = independent_paths(job_key(best_jobs),3,8,mode)
        worst_new = independent_paths(job_key(worst_jobs),3,8,mode)
        best_old = [independent_paths(job_key(best_jobs),3,8,mode,a) for a in best_assets]
        worst_old = [independent_paths(job_key(worst_jobs),3,8,mode,a) for a in worst_assets]
        pessimistic_paths = ([dict(old=0,jobs=[0]),dict(old=1,jobs=[1]),
                              dict(old=None,jobs=[2]),dict(old=None,jobs=[3])]
                             if mode == "pooled" else
                             [dict(old=0,jobs=[0]),dict(old=None,jobs=[1]),
                              dict(old=None,jobs=[2]),dict(old=None,jobs=[3])])
        optimistic_paths = [dict(old=0,jobs=[0,1]),dict(old=1,jobs=[2,3])]
        for r in box["states"]:
            a0,a1 = r["ages"]; u0,u1 = r["previous_uses"]; tau = r["delay"]
            setting = (a0,u0,a1,u1,tau)
            assert setting not in checked_keys
            checked_keys.add(setting)
            jobs = jobs_at(tau)
            assets = ((0,0,Fraction(a0),u0),(1,1,Fraction(a1),u1))
            value, states = assignment_optimum(job_key(jobs),3,8,mode,assets)
            assert value == r["new_births"]
            new = independent_paths(job_key(jobs),3,8,mode)
            olds = [independent_paths(job_key(jobs),3,8,mode,a) for a in assets]
            assert worst_new <= new <= best_new
            assert all(lo <= mid <= hi for lo,mid,hi in zip(worst_old,olds,best_old))
            common = dict(name=f"canonical_box_robust_{mode}_{setting}",jobs=jobs,K=3,L=8,
                          ownership=mode, initial_assets=assets_from_key(assets),
                          selected_paths=pessimistic_paths,introductions=box["upper"])
            simulate_witness(common)
            history_exclusions += any(a < u*tau for a,u in zip(r["ages"],r["previous_uses"]))
            values.append(dict(**r, independent_state_births=value, states_explored=states,
                               endpoint_path_nesting_verified=True,common_pessimistic_witness_valid=True))
        assert checked_keys == set(product((0,2,4),(0,1,2),(0,2,4),(0,1,2),(1,2)))
        best = assignment_optimum(job_key(best_jobs),3,8,mode,best_assets)[0]
        worst = assignment_optimum(job_key(worst_jobs),3,8,mode,worst_assets)[0]
        assert [best,worst] == [box["lower"],box["upper"]]
        assert [best,worst] == [min(r["new_births"] for r in values),max(r["new_births"] for r in values)]
        assert [best,worst] == ([0,2] if mode=="pooled" else [0,3])
        simulate_witness(dict(name=f"canonical_box_optimistic_{mode}",jobs=best_jobs,K=3,L=8,
                             ownership=mode,initial_assets=assets_from_key(best_assets),
                             selected_paths=optimistic_paths,introductions=0))
        boxes.append(dict(ownership=mode,sharp_min_max=[best,worst],state_checks=values,
                          optimistic_paths=optimistic_paths,common_pessimistic_paths=pessimistic_paths))
    assert sum(len(b["state_checks"]) for b in boxes) == extension["state_box_checks"] == 324
    print("SYNTHETIC canonical 324 box comparisons, corner extrema, nesting and common worst-state witnesses verified",flush=True)

    analytic = []
    for times,K,L in (((0,1,2,3,10),4,3),((0,1,2,10,11,12),2,3)):
        jobs = [dict(id=i,owner=0,start=t,ready=t+1,lost=False) for i,t in enumerate(times)]
        paths = independent_paths(job_key(jobs),K,L,"pooled")
        integer = assignment_optimum(job_key(jobs),K,L,"pooled")[0]
        cuts = []
        for mask in range(1,1<<len(jobs)):
            subset = {i for i in range(len(jobs)) if mask&(1<<i)}
            rank = max(len(subset.intersection(p)) for p in paths)
            cuts.append(dict(jobs=sorted(subset),rank=rank,ceiling=ceiling(Fraction(len(subset),rank))))
        if len(times)==5:
            assert integer==3 and max(len(p) for p in paths)==3
            assert next(c for c in cuts if c["jobs"]==[0,3,4])["ceiling"]==3
            selected = [dict(old=None,jobs=[0,1,2]),dict(old=None,jobs=[3]),dict(old=None,jobs=[4])]
            simulate_witness(dict(name="appendix_five_jobs",jobs=jobs,K=K,L=L,ownership="pooled",
                                  initial_assets=[],selected_paths=selected,introductions=3))
            intervals = [(jobs[p["jobs"][0]]["start"],jobs[p["jobs"][-1]]["ready"]) for p in selected]
            peak = max(sum(a<=t<b for a,b in intervals) for t in times)
            assert peak==1
            analytic.append(dict(name="five_job_subset_and_birth_vs_peak",times=times,K=K,L=L,
                                 integer_births=integer,global_rank=3,global_ceiling=2,
                                 peak_physical_presence_for_declared_retirement_policy=peak,all_subset_cuts=cuts))
        else:
            assert integer==4 and len(cuts)==63 and max(c["ceiling"] for c in cuts)==3
            # Three pair columns per isolated triple, multiplicity 1/2 each.
            fractional = {p:Fraction(1,2) for p in paths if len(p)==2}
            assert len(fractional)==6
            assert all(sum(x for p,x in fractional.items() if i in p)==1 for i in range(6))
            assert sum(fractional.values())==3
            assert all(Fraction(len(p),2)<=1 for p in paths)  # w_j=1/2 supplies LP lower bound 3
            analytic.append(dict(name="two_cluster_all_scalar_cuts_insufficient",times=times,K=K,L=L,
                                 cluster_counts=[2,2],integer_births=4,exact_fractional_LP_births="3",
                                 all_63_subset_cuts=cuts,max_single_cut=3,
                                 fractional_pairs=[dict(jobs=p,multiplicity=str(x)) for p,x in fractional.items()]))
    closed_examples=[]
    for L,wanted in ((1,(9,9)),(4,(3,6)),(7,(2,3))):
        jobs=[dict(id=i,owner=i%3,start=i,ready=i+1,lost=False) for i in range(9)]
        actual=tuple(assignment_optimum(job_key(jobs),5,L,mode)[0] for mode in ("pooled","dedicated"))
        assert actual==wanted
        closed_examples.append(dict(N=9,h=3,K=5,delta=1,tau=1,L=L,pooled=actual[0],dedicated=actual[1]))
    life_example=[]
    jobs=[dict(id=i,owner=0 if i<2 else 1,start=t,ready=t+1,lost=False) for i,t in enumerate((0,1,4,5))]
    for L,wanted in ((101,(1,2)),(2,(2,2))):
        actual=tuple(assignment_optimum(job_key(jobs),4,L,mode)[0] for mode in ("pooled","dedicated"))
        assert actual==wanted
        life_example.append(dict(L=L,pooled=actual[0],dedicated=actual[1]))
    joint=[dict(X=n,Y=3-n,cash=n+3*(3-n),environment=3*n+(3-n)) for n in range(4)]
    assert min(r["cash"] for r in joint)==min(r["environment"] for r in joint)==3
    assert not any(r["cash"]<=5 and r["environment"]<=5 for r in joint)
    budget_policies=[]
    def adjacent_covers(next_job,groups):
        if next_job==4:
            births=len(groups);switches=sum(len(g)==2 for g in groups)
            budget_policies.append(dict(paths=groups,cash=3*births+switches,environment=births+4*switches))
            return
        adjacent_covers(next_job+1,groups+[[next_job]])
        if next_job+1<4:
            adjacent_covers(next_job+2,groups+[[next_job,next_job+1]])
    adjacent_covers(0,[])
    assert min(p["cash"] for p in budget_policies)==8
    assert min(p["cash"] for p in budget_policies if p["environment"]<=4)==12
    assert any((p["cash"],p["environment"])==(12,4) for p in budget_policies)
    pilot = json.loads(PILOT.read_text(encoding="utf-8"))["results"]
    ix = {r["name"]:r for r in pilot}
    gains=[]
    for r in pilot:
        if r["ownership"]=="dedicated":
            shared=ix[r["name"].rsplit("_",1)[0]+"_pooled"]
            assert r["jobs"]==shared["jobs"] and r["initial_assets"]==shared["initial_assets"]
            simulate_witness(r,ownership="pooled")
            gains.append(r["introductions"]-shared["introductions"])
    assert (len(gains),sum(g>0 for g in gains),sum(g==0 for g in gains),min(gains),max(gains))==(289,94,195,0,4)

    review=dict(status="NO_FATAL_CONTRADICTION_IN_DECLARED_READY_ENDPOINT_SYNTHETIC_MODEL",
        completed_at_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=perf_counter()-tick,
        inputs_sha256=before,manuscript_inspection_performed=False,independent_original_path_library_use=False,
        inspected_scope=["sec:assetpaths through small tests","sec:lifecycle_results before modern-sources","app:lifecycleproofs"],
        summary=dict(closedform_rows=486,state_box_comparisons=324,distinct_box_settings=162,
                     pooled_sharp_extrema=[0,2],dedicated_sharp_extrema=[0,3],
                     direct_dedicated_witness_embedding_checks=289,
                     all_two_cluster_subset_cuts_checked=63,
                     fatal_theory_or_numeric_findings=0),
        formula_checks=formulas,state_box_checks=boxes,analytic_falsifications=analytic,
        nine_job_closed_examples=closed_examples,lifetime_erasure_example=life_example,
        simultaneous_type_budget_enumeration=joint,adjacent_cover_budget_enumeration=budget_policies,
        fatal_findings=[],
        scope_and_wording_notes=[
            dict(location="sec:assetpaths and app:lifecycleproofs",severity="scope restriction, not a fatal contradiction",
                 issue="Lost jobs use a declared detection deadline r, and r<=birth+L remains required even for terminal loss.",
                 action="Keep that explicit scenario requirement; do not infer delivery-only infeasibility, completed returns or physical detection-before-expiry necessity. A delivery-end conversion needs explicit e_j; subtracting one interval is exact only if e_j=t_j for these examples."),
            dict(location="sec:lifecycleclosed",severity="proof wording only",
                 issue="The main-text necessity sentence uses mP dispatches rather than arbitrary q; the appendix already gives the correct arbitrary-q proof.",
                 action="If editing later, say q jobs require q<=K and (q-1)*delta+tau<=L; state N,h,K integer. Cluster addition is sum of exact optima; ceil(n_b/m_b) is tight only under each cluster's regular-block premises."),
            dict(location="sec:lifecycle_results ageing comparison and appendix table",severity="attribution/label only",
                 issue="The pilot young-to-old comparison changes age and prior uses jointly; the appendix Age erases row changes lifespan, not old initial age.",
                 action="Call the former old-state deterioration, the latter lifetime erasure. To isolate age, use the verified synchronous_delay2_K2_L8_loss0 age0/age1 pair with prior0 fixed, (D,P)=(3,2)/(3,3); no smallest positive real-age erasure exists in that example."),
            dict(location="sec:lifecycleidentified and old-state support box",severity="conditional support, not physical-history validation",
                 issue="324 is 162 parameter settings evaluated in two ownership modes. Sharpness is attained extrema/interval hull, not probabilities or every intermediate integer. Some initial states fail stationary positive-delay prehistory under an extra historical rule.",
                 stationary_prehistory_exclusions_in_324_comparisons=history_exclusions,
                 action="Retain the declared synthetic forward-state support. If real history imposes age>=prior*tau or other correlated restrictions, optimize the joint admissible support before calling its bounds sharp; marginal box endpoints may be excluded."),
        ])
    assert before == {p.relative_to(ROOT).as_posix():digest(p) for p in inputs}, "An input changed concurrently; audit not saved"
    assert output_hash == digest(OUTPUT), "Side report changed concurrently; audit not saved"
    result["canonical_review"] = review
    result["summary"]["canonical_closedform_rows_independently_verified"] = 486
    result["summary"]["canonical_box_comparisons_independently_verified"] = 324
    result["summary"]["all_dedicated_witnesses_embedded_in_shared"] = 289
    OUTPUT.write_text(json.dumps(result,ensure_ascii=False,indent=2,default=str)+"\n",encoding="utf-8")
    print(json.dumps(dict(status=review["status"],summary=review["summary"],
                          input_hashes_unchanged=True,manuscript_inspection_performed=False,
                          output=str(OUTPUT)),ensure_ascii=False,indent=2),flush=True)


if __name__=="__main__":
    canonical_review() if "--canonical-only" in sys.argv else main()
