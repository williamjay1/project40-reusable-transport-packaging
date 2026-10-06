# Local release validation

The complete public runner was executed in a new working directory on 2026-10-05.
See `provenance/reproduction_validation.json` for actual run times, versions and
checks. Protected repository inputs were unchanged. No external publication occurred.

- Eight semantic smoke tests passed.
- All 578 lifecycle policies matched their declared inputs/scientific counts; actual
  witnesses and minima were independently checked. All 14 rounded-certificate gaps
  remain recorded, rather than reporting every dual as tight.
- All 486 formula rows and 324 old-state comparisons matched the frozen structural
  benchmark; the four-job foreground bridge matched without a full LCIA claim.
- All 24 historical producer/verifier stages passed. Every retained CSV benchmark
  row was compared under the predeclared tolerances, including adverse/equality cases.
- Ten main/graphical-abstract/supplementary figure sets were regenerated as 40 PDF,
  SVG, PNG and TIFF files. Actual raster DPI metadata was checked at 1200; raster
  pixel dimensions were checked. Automatic rendered collision checks passed.

The first complete attempt failed a **metadata byte-identity** check: cached JSON
reserialization changed newline bytes. This was repaired by preserving the exact
source provenance bytes. Checks were not deleted, tolerances were not widened and
scenario parameters/results were not changed. The final complete rerun passed.
The intermediate attempt completed all scientific/figure subprocesses but its report
writer failed on a TIFF DPI rational object; conversion to a JSON float repaired only
report serialization. The failure is retained in the local development audit trail.

The complete reproduction used Windows and the pinned Python 3.12 environment.
The Ubuntu GitHub Actions smoke-test/audit workflow subsequently passed as documented
below; it does not execute or establish a complete Linux scientific reproduction.
Ten regenerated standalone figure images were actually inspected for readability,
text/line overlap and clipping; the final-run versions matched the reviewed image
bytes. This is not a new manuscript page-by-page visual review, empirical validation, legal
permission, journal-index check or promise of acceptance.

## Metadata-only successor

The user confirmed **junjiezhang** as the sole human project creator on 2026-10-06.
Packaging revision 1.1 records that confirmation without changing the scientific
code, frozen inputs, benchmarks, figure data or reference artwork. Its metadata
manifest preserves the previous manifest; byte comparisons and semantic smoke tests
are reported in the separate local delivery check. This does not relabel the full
2026-10-05 computation as a new scientific rerun.

The user subsequently explicitly authorized GitHub upload. The authenticated hosting
account was verified as `williamjay1` and the source repository URL is
<https://github.com/williamjay1/project40-reusable-transport-packaging>.
This hosting update changes documentation only, not any scientific component.
Licence grants, any requested formal name formatting and
submission declarations (affiliation, funding, interests and CRediT roles) remain
to be confirmed. Authorship confirmation alone is not approval of those declarations.
DOI creation/confirmation is reserved to the user-operated Zenodo archive step.
GitHub source upload is not a formal release, archival deposit or full Linux reproduction.

## Verified GitHub source upload and smoke CI

The public source repository was verified on 2026-10-06 at commit
`0472565f126f34e8758e26f83b885abeaa5144c1`. Its author and committer names are both
`junjiezhang`. A fresh anonymous clone passed the complete repository hash/scope audit
and all eight semantic smoke tests. The remote audit covered 143 files, 106 component
lineage hashes and the declared 578-policy benchmark.

The Ubuntu/Python 3.12 GitHub Actions run
<https://github.com/williamjay1/project40-reusable-transport-packaging/actions/runs/37423668451>
completed successfully at 2026-10-06T06:26:03Z. Dependency installation, all eight
semantic tests and the repository audit passed. This is a real cross-platform smoke
check, not a new full simulation run, empirical validation, licence grant or DOI deposit.
Later documentation-only commits do not change the scientific code tested at that ref.
