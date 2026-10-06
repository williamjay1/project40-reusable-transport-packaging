"""Project40 sidecar: explicitly SYNTHETIC finite-job reachability checks.

No empirical inputs; no raw/manuscript/release access; stdout only.
Run from D: with Python -B. Exact enumeration and SciPy MILP independently
solve the same tiny set-partition instances. Not an EU physical validation.
"""

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import ceil

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp


@dataclass(frozen=True)
class Job:
    time: int
    owner: str


@dataclass(frozen=True)
class Asset:
    name: str
    birth: int
    life: int
    remaining: int
    old_id: str | None = None
    permitted_owner: str | None = None
    cash: int = 0
    env: int = 0


@dataclass(frozen=True)
class Path:
    asset: Asset
    jobs: tuple[int, ...]
    cash: int
    env: int

    @property
    def mask(self):
        return sum(1 << j for j in self.jobs)

    @property
    def births(self):
        return int(self.asset.old_id is None)


def enumerate_paths(jobs, assets, dedicated=False, lost_after=(),
                    cross_cash=0, cross_env=0):
    """Zero lead/travel, ready at delivery+1; age tested at delivery endpoint.

    Explicit birth options are part of each declared synthetic instance.
    Loss after job j terminates that asset's path (scenario known in advance).
    Old identity capacity is imposed when selecting paths, not here.
    """
    paths = []
    for asset in assets:
        for size in range(1, min(asset.remaining, len(jobs)) + 1):
            for ids in combinations(range(len(jobs)), size):
                if any(not (asset.birth <= jobs[j].time <= asset.birth + asset.life)
                       for j in ids):
                    continue
                if any(jobs[j].time < jobs[i].time + 1
                       for i, j in zip(ids, ids[1:])):
                    continue
                if any(i in lost_after for i in ids[:-1]):
                    continue
                if asset.permitted_owner is not None and any(
                    jobs[j].owner != asset.permitted_owner for j in ids
                ):
                    continue
                if dedicated and len({jobs[j].owner for j in ids}) > 1:
                    continue
                switches = sum(jobs[i].owner != jobs[j].owner
                               for i, j in zip(ids, ids[1:]))
                paths.append(Path(asset, ids, asset.cash + switches * cross_cash,
                                  asset.env + switches * cross_env))
    return paths


def new_options(jobs, life, remaining, name="new", cash=0, env=0):
    return [Asset(name, birth, life, remaining, cash=cash, env=env)
            for birth in sorted({j.time for j in jobs})]


def partitions(jobs, paths):
    """Enumerate exact disjoint covers, respecting each old asset once."""
    full = (1 << len(jobs)) - 1
    containing = {j: [p for p in paths if j in p.jobs] for j in range(len(jobs))}
    result = []

    def walk(mask, used_old, selected):
        if mask == full:
            result.append(tuple(selected))
            return
        first = next(j for j in range(len(jobs)) if not mask & (1 << j))
        for path in containing[first]:
            if mask & path.mask:
                continue
            old = path.asset.old_id
            if old is not None and old in used_old:
                continue
            next_old = used_old | {old} if old is not None else used_old
            walk(mask | path.mask, next_old, selected + [path])

    walk(0, frozenset(), [])
    return result


def score(selected, objective):
    return sum(p.births if objective == "births" else getattr(p, objective)
               for p in selected)


def admissible(selected, caps):
    return all(score(selected, name) <= limit for name, limit in caps.items())


def solve_milp(jobs, paths, objective="births", caps=None, relaxed=False):
    caps = caps or {}
    rows = [[int(j in p.jobs) for p in paths] for j in range(len(jobs))]
    lower, upper = [1.0] * len(jobs), [1.0] * len(jobs)
    for old in sorted({p.asset.old_id for p in paths if p.asset.old_id is not None}):
        rows.append([int(p.asset.old_id == old) for p in paths])
        lower.append(-np.inf)
        upper.append(1.0)
    for name, limit in caps.items():
        rows.append([p.births if name == "births" else getattr(p, name) for p in paths])
        lower.append(-np.inf)
        upper.append(float(limit))
    objective_vector = np.array([p.births if objective == "births"
                                 else getattr(p, objective) for p in paths], dtype=float)
    return milp(objective_vector,
                integrality=np.zeros(len(paths)) if relaxed else np.ones(len(paths)),
                bounds=Bounds(np.zeros(len(paths)), np.ones(len(paths))),
                constraints=LinearConstraint(np.array(rows, dtype=float), lower, upper),
                options={"presolve": True})


def crosscheck(name, jobs, paths, expected, objective="births", caps=None):
    feasible = [p for p in partitions(jobs, paths) if admissible(p, caps or {})]
    exact = min((score(p, objective) for p in feasible), default=None)
    result = solve_milp(jobs, paths, objective, caps)
    if exact is None:
        assert result.status == 2, (name, result.message)
    else:
        assert result.status == 0 and abs(result.fun - exact) < 1e-8, (name, result)
    assert exact == expected, (name, exact, expected)
    witness = next((p for p in feasible if score(p, objective) == exact), ())
    print(f"SYNTHETIC {name}: exact={exact}; MILP_status={result.status}; "
          f"columns={len(paths)}; witness="
          f"{[(p.asset.name, p.asset.birth, p.jobs) for p in witness]}")
    return witness


def rank(paths, subset, old_id=None, only_new=False):
    allowed = [p for p in paths if (not only_new or p.births)
               and (old_id is None or p.asset.old_id == old_id)]
    return max((len(set(p.jobs) & set(subset)) for p in allowed), default=0)


def scalar_bound(paths, subset):
    capacity = rank(paths, subset, only_new=True)
    return ceil(Fraction(len(subset), capacity)) if capacity else float("inf")


def peak(jobs, selected):
    # Retain an asset from birth until the final delivery's collection/loss
    # endpoint t+1, even if it has no remaining usable life. Half-open intervals.
    points = sorted({p.asset.birth for p in selected}
                    | {jobs[p.jobs[-1]].time + 1 for p in selected})
    return max((sum(p.asset.birth <= t < jobs[p.jobs[-1]].time + 1
                    for p in selected) for t in points), default=0)


def main():
    print("ALL INPUTS BELOW ARE DECLARED SYNTHETIC; not empirical EU parameters.")

    # A: three-node lifetime fragmentation, birth count differs from peak fleet.
    jobs = [Job(t, "A") for t in (0, 1, 2)]
    long = enumerate_paths(jobs, new_options(jobs, 100, 3))
    short = enumerate_paths(jobs, new_options(jobs, 1, 3))
    long_chosen = crosscheck("A_long_life", jobs, long, 1)
    chosen = crosscheck("A_short_life", jobs, short, 2)
    assert scalar_bound(short, range(3)) == 2 and peak(jobs, chosen) == 1
    stated_witness = [next(p for p in short if p.jobs == (0, 1) and p.asset.birth == 0),
                      next(p for p in short if p.jobs == (2,) and p.asset.birth == 2)]
    assert peak(jobs, stated_witness) == 1 and peak(jobs, long_chosen) == 1
    print("SYNTHETIC A: tight subset bound=2 births; constructed peak fleet=1.")

    # B: a five-node strict subset obstruction is stronger than the global rank.
    jobs = [Job(t, "A") for t in (0, 1, 2, 3, 10)]
    paths = enumerate_paths(jobs, new_options(jobs, 2, 4))
    crosscheck("B_subset_fragmentation", jobs, paths, 3)
    assert scalar_bound(paths, range(5)) == 2
    assert rank(paths, (0, 3, 4)) == 1
    assert scalar_bound(paths, (0, 3, 4)) == 3
    print("SYNTHETIC B: global rank bound=2; subset {0,3,4} bound=3 (tight).")

    # C: four-node equal service, owner cuts; age removes a strict sharing gain.
    jobs = [Job(t, owner) for t, owner in zip((0, 1, 4, 5), "AABB")]
    for life, expected_shared in ((100, 1), (1, 2)):
        assets = new_options(jobs, life, 4)
        shared = enumerate_paths(jobs, assets)
        dedicated = enumerate_paths(jobs, assets, dedicated=True)
        crosscheck(f"C_life{life}_shared", jobs, shared, expected_shared)
        crosscheck(f"C_life{life}_dedicated", jobs, dedicated, 2)
        owner_bound = sum(scalar_bound(dedicated,
                                      [j for j, job in enumerate(jobs) if job.owner == owner])
                          for owner in ("A", "B"))
        assert owner_bound == 2
        assert {(p.asset, p.jobs, p.cash, p.env) for p in dedicated} <= {
            (p.asset, p.jobs, p.cash, p.env) for p in shared}
    print("SYNTHETIC C: sharing gap falls from 1 to 0 under ageing; no reversal.")

    # D: three-node old-stock capacity, then a separate known-loss scenario.
    jobs = [Job(t, "A") for t in (0, 1, 2)]
    old = Asset("old", -2, 3, 1, old_id="old1", permitted_owner="A")
    assets = [old] + new_options(jobs, 3, 3)
    paths = enumerate_paths(jobs, assets)
    crosscheck("D_old_stock", jobs, paths, 1)
    assert rank(paths, range(3), old_id="old1") == 1
    assert rank(paths, range(3), only_new=True) == 3
    crosscheck("D_old_only", jobs, [p for p in paths if not p.births], None)
    loss_paths = enumerate_paths(jobs, new_options(jobs, 100, 3), lost_after=(0,))
    crosscheck("D_loss_after_first_job", jobs, loss_paths, 2)
    assert scalar_bound(loss_paths, range(3)) == 2

    # E: four-node cost-selected environmental nondominance, common technology.
    jobs = [Job(t, owner) for t, owner in zip(range(4), "ABAB")]
    assets = new_options(jobs, 1, 2, cash=3, env=1)
    shared = enumerate_paths(jobs, assets, cross_cash=1, cross_env=4)
    dedicated = enumerate_paths(jobs, assets, dedicated=True, cross_cash=1, cross_env=4)
    selected_shared = crosscheck("E_shared_cash", jobs, shared, 8, "cash")
    selected_dedicated = crosscheck("E_dedicated_cash", jobs, dedicated, 12, "cash")
    assert score(selected_shared, "env") == 10
    assert score(selected_dedicated, "env") == 4
    crosscheck("E_shared_env", jobs, shared, 4, "env")
    crosscheck("E_shared_joint_env4_cash12", jobs, shared, 12, "cash", {"env": 4})
    print("SYNTHETIC E: cash-selected shared (8,10), dedicated (12,4); "
          "shared still contains the dedicated environmental witness.")

    # F: three concurrent nodes, separate minima cannot certify a joint budget.
    jobs = [Job(0, owner) for owner in "ABC"]
    assets = [Asset("X", 0, 1, 1, cash=1, env=3),
              Asset("Y", 0, 1, 1, cash=3, env=1)]
    paths = enumerate_paths(jobs, assets)
    crosscheck("F_cash_only", jobs, paths, 3, "cash")
    crosscheck("F_env_only", jobs, paths, 3, "env")
    crosscheck("F_joint_budget5_5", jobs, paths, None,
               caps={"cash": 5, "env": 5})
    assert all(4 * len(p.jobs) <= p.cash + p.env for p in paths)
    assert 4 * len(jobs) > 5 + 5
    print("SYNTHETIC F: dual y=(4,4,4), lambda=(1,1); 12>10 certifies failure.")

    # G: six-node limitation (outside the 3--5-node constructive examples).
    # Two isolated odd clusters, K=2: every single-subset scalar cut admits B=3,
    # but each cluster requires two births, so the integer optimum is four.
    jobs = [Job(t, "A") for t in (0, 1, 2, 10, 11, 12)]
    paths = enumerate_paths(jobs, new_options(jobs, 2, 2))
    crosscheck("G_integer_gap", jobs, paths, 4)
    cuts = [scalar_bound(paths, ids) for size in range(1, 7)
            for ids in combinations(range(6), size)]
    assert max(cuts) == 3
    fractional = solve_milp(jobs, paths, relaxed=True)
    assert fractional.status == 0 and abs(fractional.fun - 3) < 1e-8
    crosscheck("G_birth_budget3", jobs, paths, None, caps={"births": 3})
    print("SYNTHETIC G: all 63 scalar subset cuts max=3; LP=3; integer=4.")
    print("PASS: exact covers/MILP, tight certificates, inclusion, loss, "
          "old stock, joint dual and the stated integrality limitation.")


if __name__ == "__main__":
    main()
