# Post-result joint-dependence robustness protocol

Recorded 2026-10-01 before calculating this extension. This is a disclosed post-result sensitivity extension, not preregistration or estimation of operating correlations.

## Question and fixed evidence

Does the matched, terminal-aware joint environmental/financial transport budget depend materially on assuming independent electricity and heat washing-input multipliers? The source-accounting envelopes, all 56 financial policies, transport loading/technology, delays, demand, lost-asset fate coordinates, prices and terminal assumptions remain unchanged. We retain all 216 distinct contexts at additional net cost kappa=0 and 0.005 EUR/discounted service (432 rows), including flat equality and adverse loss cases. No new optimization is required or claimed.

## Single-factor matrix

Only the joint energy distribution changes. Marginal u=lambda_E-1 and h=lambda_H-1 each have triangular density 2(1-x), x in [0,1]. Compare independent marginals (already independently verified v31 result), perfect positive ranks u=h=1-sqrt(1-z), and perfect negative ranks u=1-sqrt(1-z), h=1-sqrt(z), z uniform [0,1]. These rank pairings are stress cases, not worst/best copulas among all admissible dependence structures. In addition derive the full support guarantee: every (u,h) in [0,1]^2 must satisfy all environmental and strict financial inequalities. This last statement is distribution-free for the declared support, not a 95% empirical probability or a marginal-distribution copula bound.

## Estimands and signs

For each gamma>=0 use the same 14 upper source-impact planes and cash-saving plane as v31. Preserve the sign-dependent source-transport endpoint when replacing the 139-km source return. Compute probability of the intersection directly, not products of marginal probabilities. For rank-negative pairings the feasible set can be a union of intervals: do not assume monotonicity in z. Largest gamma with coverage>=0.95 is a supremum; strict constant equality is failure and may make a supremum unattained. Failure at gamma=0 means no nonnegative budget, not zero. Check non-increasing coverage in gamma from all checked positive transport endpoints and nonpositive cost slope. For full support check every corner, retaining signed cash slopes; failure at gamma=0 again means no budget.

## Before-compute expectations and interpretation

No direction of rank-dependence effects is asserted: environmental and saving planes can have different slope signs. A full-support guarantee, when it exists, cannot exceed a 95% threshold under either declared pairing or independence. The useful main result is a matched decision boundary and its dependence sensitivity, not universal superiority or a distribution inferred from European operations. Existence, below-reference thresholds, flat equality and any ordering violations will be reported across the complete matrix, without selecting only favorable cases.

## Numerical and independent verification

Producer: exact one-dimensional interval lengths after resolving all constraint crossings, with root inversion for gamma; the full-support bound uses all four energy corners and bisection of the minimum margin. Independent implementation: source-constraint LP support and raw cash-ledger reconstruction; analytic angular crossings for negative ranks and direct marginal CDF for positive ranks; analytically solved sign-aware corner gamma ceilings. It must not import the producer or its rank-integral routine. Analytical controls include positive/negative constant ties, diagonal and mixed-sign planes, disconnected rank-negative events, and gamma-root existence/supremum distinctions. Verify every row, classifications, baseline independent probability, units and raw/frozen input hashes. Tolerance 2e-8 in roots, 1e-10 in probabilities; stricter errors retained in report.

## Resources and execution

CPU-only; no GPU, external API cost, paid data or network acquisition. Estimated seconds to a few minutes and under 20 MB added work, on D drive (200 GB available at start). No new raw download; F/E originals remain read-only. Scripts take explicit --travel-dir, --output-dir and --results-dir paths. Fresh output directories; never overwrite frozen results. Independent verifier adds only its own report. If verified and material, update the main methods/results/discussion, numeric artifact checks and reproducibility runner; otherwise retain as negative robustness evidence without manufacturing a favorable conclusion.
