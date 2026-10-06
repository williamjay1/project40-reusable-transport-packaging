# Continuing-pool boundary test, v34

Declared 2026-10-01 before the new periodic solves. Post-result extension to v33, not preregistration. Scientific cutoff remains 2026-08-29. Earlier results and complete adverse cases remain preserved.

## Data feasibility first

| Quantity | Actually available source | Role/quality/boundary |
|---|---|---|
| 490 directed collection-to-wash distances, 14 customer positions, 35 washer positions | Existing attributed `inputs/geometry/geometry_numeric_facts.csv` and `geometry_provenance.json`; original public R3PACK D3.5 DOI10.5281/zenodo.21392137 retained read-only F archive | Same regional geometry; symmetric reverse-leg scenario, no inferred RPC traffic |
| Correlated 5-coordinate source LCIA constraints and 14 indicators | Existing public Tua2019 tables, DOI10.3390/resources8020110; attributed `tua_reported_grid.csv`, `replacement_stage_facts_v21.json`, `tua_process_contributions.csv` | Historical source accounting, not a rebuilt/current background inventory |
| Baseline reject rate .0055; manufacturing4.55, washing.07, holding.03 and truck.5/1000 | Existing source/cost provenance and historical-input audit | Historical partial cost and disclosed source recipe; no new observed operator costs |
| Demand, capacities, collection/deadlines | Existing explicitly declared engineering scenarios | Not observations: total420000 services per12week period; 2000 wash/site/week; amplitude.5 |
| Asset calendar life/continuing-stock charge | No compatible observed public fleet-age series established | Do NOT supply an empirical life. Solve an inverse charge limit instead |

All numerical core inputs are already locally available through the authorised public-data route. No new original download or private application is required. Original F/E sources and prior release files are not changed. Outputs stay D. Computation has no paid database/software dependency.

## Engineering question and allocation parameter

Does the stock/attainability mechanism survive a continuing, periodic network without empty prehistory or final liquidation? A fully adaptive periodic dedicated fleet is compared to a periodic pool, retaining the same capacities, arcs, return delay and source-recipe rejection replenishment.

Ready/dirty state equations close modulo H=12. Delivered assets remain unavailable tau weeks and washing supplies clean assets one further week later. Missing assets are replaced at the deadline; q is collected by that deadline and is not an observed eventual recovery rate. Both ownership modes see all public arcs. Fleet is conserved every periodic week, including ready, dirty, washing and in-use/not-yet-resolved pipeline states.

Let beta be the attributed number of additional recovered virgin-stock manufacture/end-of-life pairs per standing asset-week. It is a stock-footprint allocation/renewal scenario, not a measured lifetime or rejection hazard. beta=0 is a limiting inherited-stock incremental boundary that excludes its embodied footprint; beta=1/52 and1/12 are declared allocation stresses, not observed annual/quarterly fleet retirement. The added source response is beta*H*B*(p0+e0). Process rejects use the existing source recycled-replacement recipe, not this virgin-stock pair. Missing assets get virgin manufacture and no recycling credit. No double retirement of process rejects is imposed.

Inverse beta thresholds and service-equivalent N/(beta*H*B) apply to attributed capital-stock pairs, not individual physical survival/PPWR rotations. Full stock-accounting comparisons need justified real allocation/renewal data; beta=0 alone cannot establish a lifetime LCA advantage. Nominal policy selection uses beta_cost=1/12 and the retained partial costs; the environmental inverse holds each selected policy fixed. No new discounted joint-prior matrix is asserted.

## Complete matrix and fixed controls

- Main: weights uniform at phase-compressions0,.5,1 and ramp at0; tau1,3,5; q1; both ownership modes:4*3*2=24 solves.
- Collection: uniform compression0; tau1,3; q.94,.97,.99; both modes:3*2*2=12 solves.
- Start-origin invariance control: uniform compression.35; tau3; q1; origin0 versuspi/2 (exact3week shift); both modes:2*2=4 solves.
- Total **40 nominal periodic LP solves**, all retained. Evaluate three beta coordinates0,1/52,1/12 for each policy and all14 indicators; compute selected fixed-policy all-indicator inverse beta thresholds. No favorable-only filtering.

The primary expected test is algebraic: periodic owner-stock lower bounds should be invariant to an integer cyclic shift, whereas startup zero-prehistory bounds need not be. We expect partial alignment to reduce shareable reserve; we do not require pooling to dominate every indicator or calendar. At beta0, stock reduction alone contributes no embodied-footprint advantage. Any environmental claim must use the actual remaining source response and route coefficients.

## Equations, checks and execution acceptance

Periodic ready balance: R[t,s]-R[t-1,s]-new[t,s]+dispatch[t,s]-wash[t-1,s]=0. Dirty: D[t,j]-D[t-1,j]+sum_i wash[t,j,i]=q*d[t-tau,j]. Source recipe restores accepted/rejected batches together, as in v33. New-supply totals equal (1-q)*d[t-tau]. All indexes moduloH.

At end of weekt, B=sum ready+sum dirty+sum wash[t]+sum_{a=0}^{tau-1}d[t-a]; equality by owner for dedicated. Necessary bound: pooled max_t sum_j(sum_{a=0}^{tau-1}d[t-a,j]+q*d[t-tau,j]); dedicated sum_j max_t of its owner window. Queues/capacities may make the lower bound unattainable, so feasible solutions must be checked separately.

For source coordinates theta=(p0,e0,p,e,w), impact numerator/100 has coefficients((1-q)N+betaHB,betaHB,W,W,W), with source transport subtraction and actual dirty/clean movement addition; divide by the same N-service single-use reference. Reject r=.0055 is already included in source p/e slopes per wash. The independent checker must reconstruct conservation, nominal costs, source constraints and all endpoint/root comparisons without importing the producer.

Acceptance: actual feasibility/duality tests; independent modulo balances and fleet counts; rotation controls; all source endpoint/root checks; no writes to original F/E. No empirical external validation claim. Retain failures, root absences and the distinction between necessary floors and selected passing policies. Anticipated CPU-only scale<10minutes and<100MB work output; this is an estimate, not a completed-run claim. No ML training required.
