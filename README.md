# Reusable transport packaging: lifetime, return calendar and ownership

Code and numerical evidence accompanying **Effects of Asset Lifetime and Return
Delays on the Environmental Performance of Reusable Transport Packaging**
(EP&SE manuscript content version v39).

**GitHub source repository, not an archival deposit or formal software release.**
Repository: <https://github.com/williamjay1/project40-reusable-transport-packaging>.
The sole human project creator is **junjiezhang**, explicitly confirmed by the user
on 2026-10-06. The authenticated hosting account `williamjay1` is not an additional
project contributor. Author-granted licences for original materials are still pending;
public visibility does not itself grant MIT or CC BY rights.
No DOI or journal acceptance is claimed. Read `RELEASE_CHECKLIST.md` before archiving.
This packaging update changes publication/authorship documentation only; the
scientific content remains manuscript v39 and the unreleased software version remains
provisionally 1.0.0. See `AUTHORSHIP.md` for the boundaries of this confirmation.

## Scientific scope

This is a conditional engineering/theory study. The main asset model uses explicitly
synthetic service jobs, initial states, calendar lives, service quotas and known-loss
scenarios. The historical network uses published Italian source-response coefficients
and R3PACK geometry as a **primary-packaging proxy**, with declared crate demand and
capacity, not observed EU-27 transport-crate operations. A separate Spanish foreground
example multiplies four attributed scalar facts by synthetic lifecycle counts; it is
not a contemporary full LCA. No asset-event cohort, paid LCI database or confidential
operator log is supplied or inferred.

The service boundary requires availability through the final return/reconditioning
endpoint. Minimum cumulative new asset entries, concurrent stock and legal reuse
targets are different quantities. Both ownership systems receive matched service,
technology and initial options; adverse, equality and infeasible cases are retained.

Scientific literature cutoff: 2026-08-29. Full reproduction/rights audit: 2026-10-05.
User-confirmed creator metadata update: 2026-10-06; not a new literature/data search.
Source retrieval dates may be later; these are not silently treated as publication dates.

## Reproduce

Use Python3.12 and a fresh environment. Install `requirements.txt` using your normal
environment manager. The package uses SciPy/HiGHS, not a commercial solver.

```sh
python -B -m unittest discover -s tests -v
python -B tools/audit_repository.py
python -B reproduce.py --output /absolute/path/to/NEW_WORK_DIRECTORY --scope lifecycle
python -B reproduce.py --output /absolute/path/to/ANOTHER_NEW_DIRECTORY --scope all --figures
```

On the author's Windows workstation choose new output directories on the computational
drive; readers elsewhere can choose any lawful local working location. The runner copies
the repository into a new work directory, directs temporary files/caches there, hashes
and preserves repository inputs, executes the actual producers/verifiers and compares
scientific counters/tables with the frozen benchmarks. It never logs in, deposits,
uploads, assigns a DOI or downloads missing data. Assertions must remain enabled.

Scopes:

- `lifecycle`: 578 policies; independent job-to-asset state assignments, trajectory
  simulation and rational certificate checks; 486 constructive formula rows and 324
  admitted old-state comparisons; separate matched foreground example.
- `historical`: all24 producer/verifier stages, 208 nominal network policies and retained
  ancillary robust/source-bound solves; phase, periodic, financial, collection,
  electricity-component and joint-prior/travel/dependence analyses.
- `all`: both. `--figures` also regenerates four main figures, graphical abstract and
  five supplementary figures. Vector PDF/SVG plus native1200dpi PNG/TIFF are exported.
  Vectors have no DPI. Raster export is a new render, not interpolated enlargement.

An optimal schedule or font rendering can vary with the solver/platform. Tests compare
the matched scientific quantities and validate the actual trajectories; they do not
require one degenerate optimizer witness or bit-identical newly generated artwork.
Automatic collision checks do not replace visual examination at the final print size.

The separate `lifecycle_v38_independent_check.py` A–G examples use an instantaneous
delivery endpoint; they do **not** replace the pilot's final-ready endpoint validation.
The path verifier's historical `canonical_*` field names now denote the independent
structural formula/state checks. It no longer opens or claims to inspect a manuscript.
Source-PDF extraction and the retired large seafood-table appendix are not run by this
package; frozen attributed inputs remain distinguishable from newly computed outputs.

## Directory guide

| Location | Contents/data status |
|---|---|
| `code/` | Lifecycle producers, independently coded state checks and figure scripts |
| `results/lifecycle_v38/*protocol.md` | Original scenario declarations, not observations |
| `results/lifecycle_v38/modern_lci_numeric_facts.json` | Four attributed source facts and limitations |
| `historical/code/` | Source-response/calendar/LCC/uncertainty producers and verifiers |
| `historical/inputs/` | Frozen official/source extracts; provenance and original rights retained |
| `historical/reference/` | Four original scenario protocols, not an old manuscript dependency |
| `benchmarks/` | Author-computed frozen targets; not new empirical measurements |
| `figure_data/` | Full plotted numeric records, definitions and source attribution |
| `figures_reference/` | Ten checked figure sets: vector PDF/SVG masters |
| `provenance/` | Component hashes and documented portability-only changes |
| `metadata/` | Confirmed sole-human-creator record and inactive release template; no DOI |

See `DATA_DICTIONARY.md`, `DATA_SOURCES.csv` and `THIRD_PARTY_NOTICES.md`.
No full third-party article/report, supplementary artwork, ecoinvent process inventory,
personal asset record, password, internal review report or historical ZIP is included.

## Interpretation and integrity

Regular-path rules are exact only under their stated timing/ownership assumptions.
The declared state-box extrema are sharp only if the joint admissible support contains
the relevant corners; censoring alone supplies neither an upper age nor a future-life
distribution. Loss events are known absorbing scenarios, not estimated risks.
Budget-prior quantiles are conditional design-prior outputs, not confidence intervals
or empirically calibrated probabilities. Historical cost is a stated partial/cash
account, not a complete commercial contract. Generation-only electricity changes do
not update every impact category or identify a modern country-level full LCA.

Source ambiguities and adverse scenarios have not been erased to create a favourable
result. Source documents, extracted factual values and original results have different
rights; a top-level code licence must not override the source-specific notices.

## Citation and Zenodo

The creator is confirmed as `junjiezhang` and the source repository URL is recorded;
licence grants and any formal name formatting remain to be approved before a formal
licenced archival release. `CITATION.md`
explains the pending state; do not cite an invented DOI. The prepared
metadata finalizer can create `CITATION.cff` from verified author-supplied information.
It currently requires a confirmed given/family-name form; the exact unsplit name in
`metadata/confirmed_contributors.json` must not be arbitrarily split to satisfy it.
After verified metadata/licence edits, `tools/refresh_manifest.py --confirm-metadata-only`
preserves the old manifest and updates document hashes; it refuses scientific changes.
`ZENODO_GUIDE_zh.md` describes the user-operated archive step. No active `.zenodo.json`
is supplied, avoiding conflicting creator/version metadata.

AI assistance contributed to source auditing, model/verification code, figure production
and manuscript structure. Critical sources and numerical outputs require human-author
approval; reproducibility is not evidence of editorial acceptance or physical validation.
