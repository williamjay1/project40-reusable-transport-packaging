"""Read-only syntax, provenance, privacy and scope screen (not a legal opinion)."""
import argparse
import ast
import csv
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
EXCLUDED = {".git", ".venv", "__pycache__"}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def audit(root=ROOT):
    files = sorted(p for p in root.rglob("*") if p.is_file() and not EXCLUDED.intersection(p.parts))
    assert files
    syntax, scanned = [], []
    for path in files:
        rel = path.relative_to(root).as_posix()
        assert path.stat().st_size < 100_000_000, (rel, "GitHub single-file limit")
        assert path.suffix.lower() not in {".zip", ".docx", ".pem", ".key", ".p12"}, rel
        assert path.name not in {".env", "credentials.json", "id_rsa", "id_ed25519"}, rel
        if path.suffix == ".pdf":
            assert rel.startswith("figures_reference/"), (rel, "Only original figure masters allowed")
        if path.suffix.lower() in {".py", ".md", ".json", ".csv", ".tex", ".cff", ".yml", ".txt"}:
            content = path.read_text(encoding="utf-8-sig")
            if path.name != Path(__file__).name:
                assert not re.search(r"(?i)[CDEF]:[\\/](?:Users|MLWork|AcademicData|Europe environment)", content), (rel, "local absolute path")
                assert not re.search(r"\b(?:ghp_|github_pat_)[A-Za-z0-9_]{16,}|-----BEGIN (?:RSA |OPENSSH )?PRIVATE KEY-----", content), (rel, "credential-like string")
            scanned.append(rel)
            if path.suffix == ".py":
                ast.parse(content, filename=rel)
                syntax.append(rel)
    lineage = json.loads((root / "provenance/component_lineage.json").read_text(encoding="utf-8"))
    for record in lineage["components"]:
        assert sha(root / record["file"]) == record["released_sha256"], record["file"]
        assert record["byte_identical"] == (record["original_sha256"] == record["released_sha256"]), record["file"]
    with (root / "historical/inputs/geometry/geometry_numeric_facts.csv").open(encoding="utf-8-sig", newline="") as f:
        geometry = list(csv.DictReader(f))
    positive = sum(float(r["distance_km"]) > 0 for r in geometry)
    zero = sum(float(r["distance_km"]) == 0 for r in geometry)
    assert (len(geometry), positive, zero) == (490, 471, 19)
    baseline = json.loads((root / "benchmarks/lifecycle_pilot_rows.json").read_text(encoding="utf-8"))
    assert len(baseline["results"]) == 578 and baseline["declared_synthetic"]
    manifest = root / "provenance/file_manifest.json"
    if manifest.exists():
        records = json.loads(manifest.read_text(encoding="utf-8"))["files"]
        assert {r["path"] for r in records} == {p.relative_to(root).as_posix() for p in files if p != manifest}
        for r in records:
            p = root / r["path"]
            assert sha(p) == r["sha256"] and p.stat().st_size == r["bytes"], r["path"]
    return dict(all_passed=True, scope="static_release_screen_not_scientific_or_legal_certification",
                files=len(files), python_files_syntax_checked=len(syntax), text_files_privacy_scanned=len(scanned),
                component_hashes_checked=len(lineage["components"]), geometry_rows=490,
                positive_arcs=471, excluded_zero_rows=19, synthetic_policies=578,
                external_publication_performed=False, doi_assigned=False,
                file_manifest_checked=manifest.exists())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    print(json.dumps(audit(), indent=2))
