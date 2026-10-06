# Archival and licenced release gate

Local scientific reproduction and public archival publication are separate gates.
The user explicitly authorized GitHub source upload after confirming `junjiezhang`
as the sole human creator on 2026-10-06. The verified hosting account is `williamjay1`.
The repository is <https://github.com/williamjay1/project40-reusable-transport-packaging>.
Original-material licence grants and a formal archival release remain separate gates.

Before a formal licenced/archival release:

- Retain exactly one human project creator, `junjiezhang`. Supply the chosen GitHub
  versioned repository reference and any required formal name formatting; optional ORCID/affiliation must
  be real. Do not infer account ownership, surname/given names or contribution roles.
- Approve licensing only for original code/scenarios/results; retain source notices.
- Inspect `DATA_SOURCES.csv` and rights exceptions. Resolve any unclear redistribution
  permission or exclude that material with an auditable download/extraction replacement.
- Run tests, repository audit and an isolated`--scope all --figures` reproduction.
- Examine regenerated figures at print size and verify their values against the paper.
- Confirm funding/competing interests, CRediT roles and AI disclosure with the sole
  human author. Retain upstream source attribution without converting it to authorship.
- Make a versioned tag from the reviewed commit. Do not point manuscript claims at a
  mutable branch alone; record a commit and release version.
- Let the user enable the selected repository in Zenodo before the intended release
  ingestion, or use the official manual archive path.
- Use only the actually assigned version DOI in the manuscript; do not invent a DOI
  badge. No manuscript-journal DOI exists just because code has been archived.

`metadata/release_metadata.example.json` is inactive. Fill a separate JSON file with
verified facts; `tools/finalize_metadata.py` refuses empty identities or unapproved
licensing and never contacts Zenodo. It creates only CFF, after those facts and a
separately prepared author-approved root MIT `LICENSE` are supplied. No creator has
been inferred from a Windows username. The user-confirmed unsplit creator record
in `metadata/confirmed_contributors.json` is not directly accepted by the current
given/family-name finalizer; do not invent missing fields. After creator/licence/doc edits, run
`python -B tools/refresh_manifest.py --confirm-metadata-only` and the repository audit.
This preserves the old manifest and refuses to rehash scientific-code/input changes
as though they were metadata; such changes need a new version and complete run.

This package excludes the full manuscript, reviewer notes, historical ZIPs and
third-party PDFs. Keep the authoritative main DOCX/PDF and supplement in the private
submission package. A journal deposit/acceptance or verified current SCIE index is
not implied by GitHub or Zenodo publication.
