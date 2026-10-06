"""Rehash approved metadata edits while preserving the old manifest.

Scientific code, scenario inputs and benchmarks cannot be silently changed here.
This utility neither uploads nor deposits nor claims a DOI/licence approval.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def metadata_path(name):
    p = Path(name)
    return (p.suffix in {".md", ".cff"} or p.name in {"LICENSE", "DATA_LICENSE.txt"}
            or p.parts[0] == "metadata")


def refresh(root):
    current = root / "provenance/file_manifest.json"
    original_bytes = current.read_bytes()
    original = json.loads(original_bytes)
    old = {r["path"]: r["sha256"] for r in original["files"]}
    existing = {p.relative_to(root).as_posix(): p for p in sorted(root.rglob("*"))
                if p.is_file() and not {".git", "__pycache__", ".venv"}.intersection(p.parts) and p != current}
    changed = [name for name, digest in old.items()
               if name not in existing or sha(existing[name]) != digest]
    added = sorted(set(existing) - set(old))
    assert all(metadata_path(name) for name in changed + added), (
        "Only metadata/licence/doc changes allowed; scientific changes require a new complete run", changed, added)
    archive_dir = root / "provenance/manifest_history"
    archive_dir.mkdir(exist_ok=True)
    now = datetime.now(timezone.utc)
    backup = archive_dir / ("file_manifest_" + now.strftime("%Y%m%dT%H%M%S%fZ") + ".json")
    with backup.open("xb") as f:
        f.write(original_bytes)
    assert sha(backup) == hashlib.sha256(original_bytes).hexdigest()
    records = []
    for p in sorted(root.rglob("*")):
        if p.is_file() and p != current and not {".git", "__pycache__", ".venv"}.intersection(p.parts):
            records.append(dict(path=p.relative_to(root).as_posix(), bytes=p.stat().st_size, sha256=sha(p)))
    data = dict(created_at_utc=now.isoformat(), operation="author-confirmed metadata-only refresh",
                previous_manifest=backup.relative_to(root).as_posix(), previous_manifest_sha256=sha(backup),
                changed_metadata=changed, added_metadata=added,
                external_publication_performed=False, doi_assigned_by_this_tool=False, files=records)
    temporary = current.with_suffix(".tmp")
    with temporary.open("x", encoding="utf-8") as f:
        json.dump(data, f, indent=2); f.write("\n")
    os.replace(temporary, current)
    return dict(refreshed=True, previous_manifest_preserved=backup.relative_to(root).as_posix(),
                scientific_files_unchanged=True, metadata_changes=changed, metadata_additions=added)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--confirm-metadata-only", action="store_true", required=True,
                        help="Human confirms that creator/licence/repository facts are real and only metadata changed")
    parser.parse_args()
    print(json.dumps(refresh(ROOT), indent=2))
