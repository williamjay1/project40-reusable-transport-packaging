# Lifecycle boundary extension (v38)

Recorded after the 578-policy pilot and before this extension is run.
This is a disclosed, post-pilot structural extension, not preregistration.
All demand and asset-state values are declared mathematical examples.

1. Independently test the closed-form equal-spacing boundary for N=3,6,9,
   owners=1,2,3, K=2,3,5, L=2,4,8, spacing=1,2, ready delay=1,2 when
   delay<=spacing and L>=delay. Compare the formula to exact subset-DP
   covers. No old assets, losses, routing cost or shared capacity coupling.
2. Two old assets, one per owner, jobs at 0,1,4,5, owners A,A,B,B,
   new-service capacity K=3 and calendar life L=8. Enumerate all age values
   0,2,4 and previous-use values 0,1,2 independently for each old asset,
   and ready delay=1 or2. Check all 162 state realizations in each ownership
   system against their optimistic/pessimistic endpoint path families.
   The box includes all endpoint states; this tests sharp min/max support
   bounds, not probabilities and not survival estimated from censoring.
3. Demand-preserving synthetic four-job examples test whether joint
   budgets reject a cash-minimum shared policy while retaining the common
   dedicated witness. Three simultaneous jobs with two birth types test
   failure despite separately attractive cash and environmental minima.
4. No modern LCIA, legal conformity or empirical cohort validation is
   inferred. A separate source compatibility audit must precede physical
   projections. Finite observed use histories constrain attained past
   service; they do not identify these future-life support endpoints.

The ready endpoint is deliberate: calendar expiry covers final return and
reconditioning readiness. New assets enter at first dispatch. Old assets
are ready at time0. Births are cumulative entries, not concurrent stock.
Both ownership systems receive the same initial resources and options.
