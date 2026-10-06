"""Create software CFF only from confirmed real author/repository/licence facts.

No upload, login, release, archival or DOI creation. Run after human-author review.
"""
import argparse
from datetime import date
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]


def render(data):
    assert data.get("author_approved_original_code_MIT") is True, "Original code MIT approval required"
    assert data.get("author_approved_original_results_CC_BY_4_0") is True, "Original results approval required"
    assert data.get("third_party_rights_review_approved") is True, "Mixed-rights review required"
    url = data.get("repository_url") or ""
    assert re.fullmatch(r"https://github\.com/[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+/?", url), "Actual GitHub repository URL required"
    authors = data.get("authors") or []
    assert authors and data.get("copyright_holders"), "Real creators and copyright holders required"
    assert data.get("title") and re.fullmatch(r"\d+\.\d+\.\d+", data.get("version", "")), "Title and semantic release version required"
    date.fromisoformat(data["release_date"])
    # No DOI is accepted here: obtain it by actual user-operated archive, then add
    # the verified identifier to a later metadata update.
    assert not data.get("version_doi") and not data.get("concept_doi"), "This pre-deposit finalizer does not assert a DOI"
    q = lambda text: json.dumps(text, ensure_ascii=False)
    lines = ['cff-version: "1.2.0"', 'type: software',
             'message: "Please cite this version of the reproducible research software."',
             "title: " + q(data["title"]), "version: " + q(data["version"]),
             "date-released: " + q(data["release_date"]), "repository-code: " + q(url),
             "license: MIT", "authors:"]
    for author in authors:
        assert author.get("given-names") and author.get("family-names"), "Real given/family names required"
        lines.extend(["  - family-names: " + q(author["family-names"]),
                      "    given-names: " + q(author["given-names"])])
        if author.get("orcid"):
            assert re.fullmatch(r"https://orcid\.org/\d{4}-\d{4}-\d{4}-\d{3}[\dX]", author["orcid"]), "ORCID syntax"
            lines.append("    orcid: " + q(author["orcid"]))
        if author.get("affiliation"):
            lines.append("    affiliation: " + q(author["affiliation"]))
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--metadata", type=Path, required=True)
    args = p.parse_args()
    data = json.loads(args.metadata.read_text(encoding="utf-8"))
    text = render(data)
    # Publication remains blocked until author creates actual LICENSE text and
    # source-specific data notices. Never fabricate a grant from a template.
    assert (ROOT / "LICENSE").is_file(), "Install the author-approved MIT LICENSE before activating citation metadata"
    with (ROOT / "CITATION.cff").open("x", encoding="utf-8") as f:
        f.write(text)
    print("CITATION.cff created from supplied facts; no deposit or DOI performed.")
