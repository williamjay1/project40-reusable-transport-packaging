# Calendar phase/size mechanism extension — before execution

2026-10-01, Project 40. Post-result extension of the existing ownership mechanism, not preregistration and not observed European demand.

Question: does the owner-specific mandatory-stock mechanism survive partial demand synchronisation and unequal customer means, rather than only the equal-weight perfectly staggered example?

Fixed inputs: existing 14 service locations / 35 wash sites / 490 public distance facts; source W1-RE1, full collection, entirely new opening stock, 12 weeks, period 12, mean 35,000 services/week, amplitude 0.5, one preparation week, historical reference partial-cost coefficients, terminal source liquidation. No source coefficients or raw files change.

Demand: d[t,j] = 35000*w[j]*(1 + 0.5*sin(2*pi*t/12 + phase_origin + (1-c)*2*pi*j/14)). c=0 is the previously staggered case; c=1 aligned. Weights: uniform; ramp proportional to 1..14; reverse ramp. Origins: 0, pi/2, 3*pi/2. Bounds-only grid: c=0..1 in steps .05, tau=1,3,5, all weights/origins (1701 rows). All rows retained.

Exact finite-campaign stock bounds use non-circular, zero-pre-horizon rolling demand. R=abs(sum_j w[j]*exp(i*phase[j])) is only an aggregate harmonic coherence descriptor; it does not uniquely identify finite-campaign dedicated stock or capacity/routing feasibility. No R-only universal threshold will be asserted.

Constructive network matrix: uniform and ramp weights; origin 0; tau=1,3; c=0,.5,.9,1; fully adaptive dedicated and pooled ownership (32 new nominal-cost LP policies). Report every policy, including failures and non-advantage. Bounds are necessary, not sufficient; constructive policies are separately verified. Ozone and all 14 indicators are retained. No inference of global minimum collection, operator LCC, or empirical EU probability.

Validation: independently rebuild demand/rolling bounds from explicit arithmetic, verify every flow and source endpoint with the independent existing event ledger, and reconstruct source floor minima separately. Endpoint c=0/uniform and c=1 controls compare against historical schedules. Demand totals and all 1701 rows, all 32 feasible schedules, source conservation and numerical residuals must pass. Any mismatch is resolved before manuscript use. No selected-policy environmental failure is promoted to an all-policy failure unless the mandatory-stock climate floor and nonnegative source remaining terms support it.

Scientific cutoff remains 2026-08-29; this is new computation with existing disclosed inputs, not a later scientific observation.

## Arithmetic-count erratum after first execution

The coordinate matrix written above contains 3 weights × 3 origins × 3 delays × 21 compression values = **567**, not 1701, bounds-only rows. The original prose multiplication was wrong. All declared coordinates were executed; no coordinate, policy, parameter or retention rule changed. The 32 new LP policies and 16 paired comparisons are unchanged. Verification must check all 567 rows rather than inventing 1134 additional observations.
