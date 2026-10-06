# Project40 lifecycle reachability pilot

Recorded before pilot execution on 2026-10-02. This is a deterministic,
declared small-network proof and computational test, not observed RPC demand.
The existing 471-arc calculation and v37 release are preserved unchanged.

Primary question: can an ownership system cover every specified service job
within both a birth-material budget and a cash-equivalent birth budget when
asset use-capacity, calendar expiry, ready delay and loss stop future reuse?
Minimum cumulative asset introductions are not maximum concurrent fleet.

Three owners each request three identical-capacity services. Synchronous
starts are 0,3,6; staggered starts are 3*m+owner, m=0,1,2. One time unit is
a declared planning interval, not a fitted return distribution. Normal
ready delay is 2 or 4 intervals. Each job ends at its specified ready time.
An asset cannot next serve before that time. No routing-cost coefficient is
inferred: this pilot's transition delays represent fixed admissible mobility.
Both ownership systems have the same transition opportunities; dedicated
paths cannot change owner.

Asset maximum services K={2,3,5,9}; calendar lifespan L={4,8,12} intervals;
loss stops={none,jobs2and6}; initial stocks={none,young,old}. Young stocks
are three known assets, one per original owner, with zero previous services
and age0. Old stocks have K-1 previous services and ageL-3: their remaining
capacity is1service and calendar remaining time3. Initial assets are optional
resources, not compulsory purchases; prehistory impacts/cash are not charged
as new births. Their original ownership constrains dedicated reuse only.

The full matrix is2calendars*2delays*4K*3L*2loss patterns*3initial states
*2ownership systems=576exact integer policies. All adverse results are retained.
Life/age/loss inputs are scenario limits, not lifetime estimates.

Each admissible complete single-asset path is enumerated. A binary set
partition MILP covers every job exactly once and uses each initial asset at
most once. Introductions are charged at first dispatch (zero procurement
lead time). This tests existence with known jobs/scenario; it is not an
online or nonanticipative uncertainty policy.

A separately solved nonnegative weighted-job dual provides a necessary
introduction certificate. Float weights are converted to rational values,
scaled to ensure every new path has weight<=1, and each old-resource
allowance is recomputed as its exact maximum path weight. The resulting
inequality is checked exactly with fractions. The rounded certificate is
compared to the integer optimum and an independently coded subset DP.
No claim that every certificate is tight is made in advance.

Predeclared extra falsification instance: one long unavailable asset job
[0,8) and six one-interval jobs starting1,...,6, distributed over three
owners; K=4,L=20. Test whether simple maximum occupancy and total-use
ceil(N/K) can both underestimate required births. The weighted certificate
is evaluated, not tailored by changing demand after results.

Primary records: integer cover, exact certificate, independent DP optimum,
per-asset path, actual introductions, reusable transitions, peak occupied
stock, per-owner service completion and loss-terminal count. Secondary:
conditional joint feasibility for explicitly declared introduction ceilings
{1.5,2.5,3.5,4.5,5.5,8.5}; these are normalized resource budgets, not EUR,
climate effects or statutory thresholds. Modern process vectors can later
map introductions and washing counters within a separate matched functional
unit; old and modern sources must never be blended without provenance.

Identification extension: unobserved remaining capacities are bounded,
not fitted. Pessimistic/optimistic resource-path sets bracket feasibility;
existence for each scenario is distinguished from one policy for all scenarios.
Right-censored assets provide lower attained histories, not known future life.

CPU, <=9normal jobs, <=7extra jobs, no commercial solver. Outputs and
caches stay on D; new raw source metadata stay on F. No empirical survival
claim, global-first claim or high-readiness score follows from pilot PASS.
