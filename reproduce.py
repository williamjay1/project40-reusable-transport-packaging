"""Reproduce in a NEW working directory without altering repository inputs.

No download, GitHub upload, Zenodo login or DOI creation is performed.
Run with assertions enabled: python -B reproduce.py --output NEW_ABSOLUTE_PATH.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from PIL import Image

PACKAGE = Path(__file__).resolve().parent
LIFECYCLE = ("run_lifecycle_v38.py", "extend_lifecycle_v38.py",
             "lifecycle_foreground_v38.py", "lifecycle_v38_independent_check.py",
             "lifecycle_v38_path_verify.py")
REQUIRED_HISTORICAL = (
    "barometer_crate_bins.csv", "barometer_crate_bins_manifest.json",
    "energy/electricity_nrg_pc_205_frozen.json", "energy/gas_nrg_pc_203_EXCLUDED_after_cutoff.json",
    "financial_reference/financial_input_audit.json", "financial_reference/industrial_energy_contexts.csv",
    "geometry/geometry_numeric_facts.csv", "geometry/geometry_provenance.json",
    "JRC_generation_2026_acquired_20261001_v1.xlsx", "process_contributions_provenance.json",
    "r3pack_schema_reference.json", "replacement_stage_facts_v21.json", "source_parameters.json",
    "source_provenance.json", "tua_process_contributions.csv", "tua_reported_grid.csv")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def protected_files(root):
    return {p.relative_to(root).as_posix(): sha(p) for p in sorted(root.rglob("*"))
            if p.is_file() and not any(x in p.parts for x in (".git", ".venv", "__pycache__"))}


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def step(work, env, report, script, *args):
    name = f"{len(report['runs']) + 1:02d}_{Path(script).stem}"
    log = work / "logs" / (name + ".log")
    command = [sys.executable, "-B", str(work / script), *map(str, args)]
    tick = time.monotonic()
    print("START " + name, flush=True)
    with log.open("x", encoding="utf-8") as stream:
        process = subprocess.run(command, cwd=work, env=env, stdout=stream,
                                 stderr=subprocess.STDOUT, check=False)
    report["runs"].append(dict(script=script, arguments=list(map(str, args)),
                              exit_code=process.returncode, seconds=time.monotonic() - tick,
                              log=log.relative_to(work).as_posix()))
    if process.returncode:
        print(log.read_text(encoding="utf-8", errors="replace")[-5000:], file=sys.stderr)
        raise RuntimeError("Actual command failed: " + script)
    print("PASS " + name, flush=True)


def check_lifecycle(work):
    actual = read_json(work / "results/lifecycle_v38/pilot/pilot_results.json")
    old = read_json(work / "benchmarks/lifecycle_pilot_rows.json")
    expected = {r["name"]: r for r in old["results"]}
    observed = {r["name"]: r for r in actual["results"]}
    assert set(expected) == set(observed) and len(observed) == 578
    fields = ("K", "L", "ownership", "state", "jobs", "initial_assets", "feasible",
              "introductions", "dp_optimum", "certificate_gap", "occupied_peak",
              "simple_global_use_bound", "ready_transitions", "loss_terminal_jobs", "path_count")
    for name in expected:
        for field in fields:
            assert expected[name][field] == observed[name][field], (name, field)
    # Optimal witnesses may differ across solver versions. Independently validate
    # the actual witness rather than impose one arbitrarily selected schedule.
    independent = read_json(work / "results/lifecycle_v38/path_verification.json")
    assert independent["summary"]["independent_minimum_births_matches"] == 578
    assert independent["summary"]["weighted_certificate_gaps"] == 14
    assert independent["summary"]["canonical_closedform_rows_independently_verified"] == 486
    assert independent["summary"]["canonical_box_comparisons_independently_verified"] == 324
    extension = read_json(work / "results/lifecycle_v38/extension/structural_results.json")
    expected_ext = read_json(work / "benchmarks/lifecycle_structural.json")
    for field in ("formula_checks", "state_box_checks", "formula_rows", "state_boxes"):
        assert extension[field] == expected_ext[field], field
    gains = extension["pilot_summary"]
    assert (gains["pairs"], gains["positive_sharing_pairs"], gains["equality_pairs"]) == (289, 94, 195)
    foreground = read_json(work / "results/lifecycle_v38/foreground_v2/foreground_bridge.json")
    assert foreground == read_json(work / "benchmarks/foreground_bridge.json")
    return dict(policies=578, matched_declared_inputs_and_scientific_counters=True,
                independently_validated_witnesses=578, certificate_integer_gaps=14,
                formula_rows=486, box_comparisons=324, matched_foreground=True)


def check_historical(work):
    summary = read_json(work / "historical/logs/all_summary.json")
    assert summary["all_passed"] and len(summary["steps"]) == 24
    assert all(x["exit_code"] == 0 for x in summary["steps"])
    expected = work / "benchmarks/historical"
    comparisons = []
    # Predeclared tolerances: source endpoints 1e-8 relative, aggregate counts
    # 1e-4 absolute, EUR totals 1e-5. Preserve absence/status, never turn None to0.
    for path in sorted(expected.rglob("*.csv")):
        relative = path.relative_to(expected)
        target = work / "historical/results" / relative
        with path.open(encoding="utf-8-sig", newline="") as a, target.open(encoding="utf-8-sig", newline="") as b:
            ra, rb = csv.DictReader(a), csv.DictReader(b)
            assert ra.fieldnames == rb.fieldnames, str(relative)
            left, right = list(ra), list(rb)
        assert len(left) == len(right), str(relative)
        max_error = 0.0
        for i, (a, b) in enumerate(zip(left, right)):
            for key in a:
                if a[key] == b[key]:
                    continue
                try:
                    x, y = float(a[key]), float(b[key])
                except (ValueError, TypeError):
                    raise AssertionError((str(relative), i, key, a[key], b[key]))
                assert math.isfinite(x) and math.isfinite(y), (relative, i, key)
                atol = 1e-4 if any(s in key for s in ("assets", "fleet", "stock", "services", "inputs")) else 1e-5 if any(s in key for s in ("cost", "eur", "objective")) else 1e-8
                error = abs(x - y)
                assert error <= atol + 1e-8 * max(abs(x), abs(y)), (relative, i, key, x, y, error)
                max_error = max(max_error, error)
        comparisons.append(dict(table=relative.as_posix(), rows=len(left), maximum_absolute_error=max_error))
    return dict(all24_producer_verifier_steps_pass=True, tables=comparisons,
                no_empirical_EU_fleet_validation=True)


def check_figures(work):
    records = []
    roots = [work / "figures"]
    # Both redraws emit five named sets into this directory. Check real files,
    # not declared export DPI.
    for root in roots:
        assert root.is_dir(), root
        for pdf in sorted(root.glob("*.pdf")):
            for suffix in (".pdf", ".svg", ".png", ".tiff"):
                path = pdf.with_suffix(suffix)
                assert path.is_file() and path.stat().st_size > 0, path
                row = dict(file=path.relative_to(work).as_posix(), sha256=sha(path))
                if suffix in (".png", ".tiff"):
                    with Image.open(path) as im:
                        dpi = im.info["dpi"]
                        assert all(abs(v - 1200) < .02 for v in dpi), (path, dpi)
                        assert im.width > 1000 and im.height > 1000, path
                        row.update(width_px=im.width, height_px=im.height, actual_dpi=[float(v) for v in dpi])
                records.append(row)
    assert len(records) == 40, len(records)
    return dict(exports=len(records), native_raster_dpi=1200, files=records,
                visual_review_not_automatic=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--scope", choices=("all", "lifecycle", "historical"), default="all")
    parser.add_argument("--figures", action="store_true", help="Requires lifecycle scope; export native1200dpi and vectors")
    args = parser.parse_args()
    work = args.output.resolve()
    if not args.output.is_absolute() or work.exists() or work == PACKAGE or work.is_relative_to(PACKAGE):
        parser.error("Choose a NEW absolute working directory outside the repository.")
    if args.figures and args.scope == "historical":
        parser.error("Figure regeneration requires lifecycle results.")
    if not __debug__:
        parser.error("Assertions must be enabled; do not use Python -O.")
    if args.scope in ("all", "historical"):
        missing = [name for name in REQUIRED_HISTORICAL
                   if not (PACKAGE / "historical/inputs" / name).is_file()]
        if missing:
            parser.error("Missing frozen inputs: " + ", ".join(missing) + "; no replacement download is allowed.")
    before = protected_files(PACKAGE)
    shutil.copytree(PACKAGE, work, ignore=shutil.ignore_patterns(".git", ".venv", "__pycache__"))
    for name in ("logs", "temp", "cache", "historical/logs", "historical/temp", "historical/cache"):
        (work / name).mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env.update(PYTHONUTF8="1", PYTHONDONTWRITEBYTECODE="1", TMP=str(work / "temp"), TEMP=str(work / "temp"),
               MPLCONFIGDIR=str(work / "cache/matplotlib"), OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
    report = dict(all_passed=False, scope=args.scope, started_at_utc=datetime.now(timezone.utc).isoformat(),
                  package_files_sha256=before, runs=[], declared_synthetic=True,
                  empirical_asset_data=False, paid_LCA_database_required=False,
                  external_deposit_performed=False, doi_assigned=False,
                  runtime=dict(python=sys.version.split()[0], **{name: importlib.metadata.version(name)
                               for name in ("numpy", "scipy", "matplotlib", "openpyxl", "Pillow")}))
    try:
        if args.scope in ("all", "lifecycle"):
            for name in LIFECYCLE:
                step(work, env, report, "code/" + name)
            step(work, env, report, "code/lifecycle_v38_path_verify.py", "--canonical-only")
            report["lifecycle"] = check_lifecycle(work)
        if args.scope in ("all", "historical"):
            step(work, env, report, "historical/run_geometry_refit_v36.py", "--family", "all")
            report["historical"] = check_historical(work)
        if args.figures:
            step(work, env, report, "code/plot_epse_nature_v39.py")
            step(work, env, report, "code/redraw_supplement.py")
            report["figures"] = check_figures(work)
        report["repository_inputs_unchanged"] = before == protected_files(PACKAGE)
        assert report["repository_inputs_unchanged"]
        report["all_passed"] = True
    except Exception as exc:
        report["error"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        report["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
        report["repository_inputs_unchanged"] = before == protected_files(PACKAGE)
        with (work / "reproduction_report.json").open("x", encoding="utf-8") as stream:
            json.dump(report, stream, indent=2, allow_nan=False); stream.write("\n")
        print(json.dumps({k: v for k, v in report.items() if k != "package_files_sha256"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
