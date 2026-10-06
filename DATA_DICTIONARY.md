# Data and estimand dictionary

The JSON/CSV files contain explicit names, units, source keys and/or locators.
This document defines their evidence status and the principal fields needed to
interpret reproduction. It does not label simulations as raw operator data.

## Synthetic lifecycle records

`benchmarks/lifecycle_pilot_rows.json` contains578 author-computed policy rows.
It is a verification target, not an input to the optimizer. The scenario generator is
`code/run_lifecycle_v38.py`; its declaration is `results/lifecycle_v38/pilot_protocol.md`.

| Field | Definition/unit |
|---|---|
| `jobs[].id`, `owner` | Synthetic job identity and ownership group |
| `jobs[].start`, `ready` | Declared planning-interval endpoints; not observed timestamps |
| `jobs[].lost` | Known absorbing scenario stop; ready marks the declared end/detection deadline |
| `K`, `L` | Integer maximum services; ready-end calendar life in intervals |
| `initial_assets[].age` | Initial calendar age in intervals, including waiting |
| `initial_assets[].previous_services` | Prior services; not charged again as new manufacture |
| `introductions` | Minimum cumulative new entries at first dispatch (count) |
| `occupied_peak` | Maximum simultaneously occupied jobs, not total purchased fleet |
| `selected_paths` | A feasible exact-cover witness; old IDs can be used at most once |
| `certificate.value_fraction`, `ceil` | Exact rational necessary lower bound and rounded bound |
| `certificate_gap` | Integer introductions minus rounded certificate;14 gaps retained |

`benchmarks/lifecycle_structural.json` stores486 formula rows and162 initial-state
settings evaluated in each of two ownership systems (=324 comparisons). The interval
is a conditional hull of attainable integer counts, not a confidence interval and
not necessarily a contiguous identified set. The independent audit regenerates state
trajectories and checks the actual MILP witness; subset DP alone is not an independent
lifecycle semantics check when it shares a path library.

## Historical source-response model

`historical/inputs/tua_reported_grid.csv`: printed environmental response grid;
multiple source indicators including CED, not14 EF3.1 categories. `tua_process_contributions.csv`
and `replacement_stage_facts_v21.json`: published process components and rounding
intervals. Preserve native indicator units and source functional unit/normalisation;
do not substitute a modern factor for all categories based on one carbon ratio.

`historical/inputs/geometry/geometry_numeric_facts.csv`:14 origins×35 indexed candidate
wash positions,490 reported distances(km). Nineteen exact zeros are excluded in the
declared conservative calculation;471 positive pairs reach34 positions. Reverse clean
transport symmetry,68,000 units/week accessible capacity and crate demand are explicit
scenarios. This is not geolocated, observed RPC movement or proof that zeros are missing.

`historical/inputs/energy/electricity_nrg_pc_205_frozen.json`: JSON-stat dimensions
`freq,siec,nrg_cons,unit,tax,currency,geo,time`; monetary use selects the specified EUR/kWh
industrial band/tax basis. Gas marked`EXCLUDED_after_cutoff` must not enter costs.

JRC workbook: source Table3_LC_GHG, explicit EU27 whitelist, exact year2024; generation
LC-GHG units tonnes CO2eq/MWh are numerically equal to kgCO2eq/kWh. The calculation changes
only the generation-related climate component; transmission/delivery gap remains unknown.

`benchmarks/historical/*/*.csv` are author-computed comparison targets. The runner
regenerates the network policies from input files, executes independent event/source
and cost checks, and compares all retained benchmark rows, including missing roots,
adverse classifications and finite/infeasible distinctions. Positive-flow ledgers are
created during execution rather than treated as empirical flows.

## Foreground and figures

`modern_lci_numeric_facts.json`: four attributed published scalar facts, with unit,
source DOI, page/table locator and condition. `foreground_bridge.json`: four matched
10kg-fish RPC services with8kg common ice; K4, starts0,1,4,5; calendar life101 or2.
The avoided1.2kg and0.31kWh in the long-life example are gross new-crate mass and
manufacture electricity only, not net environmental impacts, neat-HDPE avoidance,
full energy or superiority to single use. Wash counts are4 in each ownership system.

`figure_data/` stores all numerical series, plotted comparisons and scope warnings.
`figure_semantics.tex` is a minimal author-created semantic input, not a submitted
paper. `figures_reference/` contains vector outputs checked for the v39 manuscript.
Original scenario q/q²/q³ values are design assumptions, not fitted probabilities.

Seeds, factorial sets, grids, terminal conditions and source rounding are in the
retained code/protocols. Fixed priors and probability quadrature must not be presented
as likelihood-based inference from observed cohorts.
