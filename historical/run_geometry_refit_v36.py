"""Full positive-arc primary refit; old reference thresholds are not overwritten."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import traceback
from datetime import datetime, timezone

ROOT=Path(__file__).resolve().parent
ORIGINAL=ROOT/'reference'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--family',choices=['all','adaptive','phase','periodic','financial','collection','joint'],default='all')
    a=ap.parse_args()
    if not __debug__:
        raise ValueError('Assertions must be enabled')
    label=a.family
    summary=ROOT/'logs'/f'{label}_summary.json'
    if summary.exists():
        raise FileExistsError('Preserve prior run summary; use a fresh family/run workspace')
    inputs=ROOT/'inputs'; code=ROOT/'code'; stage=inputs/'replacement_stage_facts_v21.json'; results=ROOT/'results'
    env=os.environ.copy();env.pop('PYTHONPATH',None)
    env.update(PYTHONUTF8='1',PYTHONDONTWRITEBYTECODE='1',TEMP=str(ROOT/'temp'),TMP=str(ROOT/'temp'),
               MPLCONFIGDIR=str(ROOT/'temp'/'matplotlib'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    before={p.relative_to(inputs).as_posix():sha(p) for p in inputs.rglob('*') if p.is_file()}
    original_before={p.relative_to(ORIGINAL).as_posix():sha(p) for p in ORIGINAL.rglob('*') if p.is_file() and 'tests' not in p.relative_to(ORIGINAL).parts}
    start=time.monotonic()
    report=dict(all_passed=False,family=label,steps=[],started_at_utc=datetime.now(timezone.utc).isoformat(),
                restriction='All exact-zero arcs forbidden, no capacity redistribution, immutable numeric inputs',
                declared_protocol='analysis_protocol.md',new_raw_acquisition=False)
    def save():
        report['elapsed_seconds']=time.monotonic()-start
        summary.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    def step(name,args):
        command=[sys.executable,'-B',*map(str,args)]
        item=dict(name=name,command=command,status='running');report['steps'].append(item);save()
        print('START '+name,flush=True)
        log=ROOT/'logs'/f'{label}_{len(report["steps"]):02d}_{name}.log';began=time.monotonic()
        with log.open('x',encoding='utf-8') as stream:
            process=subprocess.run(command,cwd=ROOT/'temp',env=env,stdout=stream,stderr=subprocess.STDOUT,check=False)
        item.update(exit_code=process.returncode,elapsed_seconds=time.monotonic()-began,log=log.relative_to(ROOT).as_posix(),
                    status='passed' if process.returncode==0 else 'failed');save()
        if process.returncode:
            raise RuntimeError(f'{name} failed; see {log}')
        print('PASSED '+name,flush=True)
    def fresh(name,protocol=None):
        out=results/name;out.mkdir(parents=True,exist_ok=False)
        if protocol:
            text=(ORIGINAL/'results'/protocol/'analysis_protocol.md').read_text(encoding='utf-8')
            (out/'analysis_protocol.md').write_text('# Conservative v36 geometry refit\n\nThe root analysis_protocol.md fixes the new positive-arc mask and all unchanged scenario factors. The prior scenario declaration is retained below; its original-geometry outcomes are not current controls.\n\n'+text,encoding='utf-8')
        return out
    base=['--input-dir',inputs,'--stage-facts',stage]
    flowverify=['--source-grid',inputs/'tua_reported_grid.csv','--stage-facts',stage,'--process-facts',inputs/'tua_process_contributions.csv']
    try:
        if label in ['all','adaptive']:
            out=fresh('calendar_network_v23_adaptive')
            step('adaptive_10_policies',[code/'calendar_pool_network_v23.py',*base,'--geometry-cache',inputs/'geometry','--output-dir',out,'--adaptive-benchmark'])
            step('adaptive_independent_event_source',[code/'verify_calendar_pool_v23.py','--results-dir',out,*flowverify])
            thresholds=fresh('calendar_thresholds_v23')
            step('adaptive_operating_roots',[code/'calendar_operating_thresholds_v23.py',*base,'--network-results',out/'calendar_network_results.json','--output-dir',thresholds])
            step('independent_adaptive_roots',[code/'verify_operating_thresholds_v24.py',*flowverify,'--network-results',out/'calendar_network_results.json','--thresholds',thresholds/'fixed_policy_probabilistic_travel_thresholds.csv','--output',thresholds/'independent_threshold_verification.json'])
        if label in ['all','phase']:
            out=fresh('calendar_phase_v33','calendar_phase_v33')
            step('phase_32_policies_567_bounds',[code/'run_calendar_phase_v33.py','--package-root',ROOT,'--output-dir',out])
            step('phase_independent',[code/'verify_calendar_phase_v33.py','--package-root',ROOT,'--results-dir',out])
        if label in ['all','periodic']:
            out=fresh('continuing_pool_v34','continuing_pool_v34')
            step('periodic_40_policies',[code/'run_continuing_pool_v34.py','--package-root',ROOT,'--output-dir',out])
            step('periodic_stock_certificates',[code/'derive_periodic_stock_certificates_v34.py','--package-root',ROOT,'--results-dir',out])
            step('periodic_independent',[code/'verify_continuing_pool_v34.py','--package-root',ROOT,'--results-dir',out,'--conservative-positive-arcs'])
        if label in ['all','financial']:
            prices=fresh('financial_inputs_v29')
            step('financial_public_price_decode',[code/'prepare_financial_inputs_v29.py','--electricity-raw',inputs/'energy/electricity_nrg_pc_205_frozen.json','--gas-raw',inputs/'energy/gas_nrg_pc_203_EXCLUDED_after_cutoff.json','--output-dir',prices])
            out=fresh('calendar_financial_v29')
            step('financial_56_policies',[code/'calendar_financial_completion_v29.py',*base,'--geometry-cache',inputs/'geometry','--price-dir',prices,'--output-dir',out])
            step('financial_independent_events',[code/'verify_calendar_pool_v23.py','--results-dir',out,*flowverify])
            step('financial_independent_cash_source',[code/'verify_financial_completion_v29.py','--results-dir',out,*base,'--electricity-raw',inputs/'energy/electricity_nrg_pc_205_frozen.json'])
            step('financial_joint_432_contexts',[code/'reevaluate_financial_joint_v29.py','--results-dir',out,*base])
        if label in ['all','collection']:
            out=fresh('calendar_collection_v25')
            step('collection_70_policies',[code/'calendar_collection_frontier_v25.py',*base,'--geometry-cache',inputs/'geometry','--output-dir',out])
            step('collection_independent_events',[code/'verify_calendar_pool_v23.py','--results-dir',out,*flowverify])
            grid=fresh('calendar_grid_component_v25')
            raw=inputs/'JRC_generation_2026_acquired_20261001_v1.xlsx'
            step('generation_component_scenarios',[code/'calendar_grid_component_screen_v25.py',*base,'--raw-file',raw,'--network-results',results/'calendar_network_v23_adaptive/calendar_network_results.json','--output-dir',grid])
            step('collection_generation_independent',[code/'verify_collection_grid_v25.py',*base,'--collection-dir',out,'--grid-dir',grid,'--grid-raw',raw,'--adaptive-results',results/'calendar_network_v23_adaptive/calendar_network_results.json','--output',out/'independent_collection_grid_verification.json'])
        if label in ['all','joint']:
            financial=results/'calendar_financial_v29';out=fresh('joint_integration_v30')
            step('joint_integration_432',[code/'integrate_joint_prior_v30.py','--results-dir',financial,*base,'--output-dir',out])
            step('joint_independent_integration',[code/'verify_joint_integration_v30.py','--results-dir',out,'--financial-dir',financial,*base])
            travel=fresh('joint_travel_v31','joint_travel_v31')
            step('joint_transport_432',[code/'run_joint_travel_v31.py','--financial-dir',financial,*base,'--integration-dir',out,'--output-dir',travel])
            step('joint_independent_transport',[code/'verify_joint_travel_v31.py','--financial-dir',financial,*base,'--integration-dir',out,'--results-dir',travel])
            dependence=fresh('joint_dependence_v32','joint_dependence_v32')
            step('joint_dependence_432',[code/'run_joint_dependence_v32.py','--travel-dir',travel,'--output-dir',dependence])
            step('joint_independent_dependence',[code/'verify_joint_dependence_v32.py','--financial-dir',financial,*base,'--travel-dir',travel,'--results-dir',dependence])
        after={p.relative_to(inputs).as_posix():sha(p) for p in inputs.rglob('*') if p.is_file()}
        assert before==after,'Numeric input mutated'
        assert all((ORIGINAL/p).is_file() and sha(ORIGINAL/p)==h for p,h in original_before.items()),'Historical release mutated'
        report.update(all_passed=True,unchanged_input_files=len(before),unchanged_historical_release_files=len(original_before))
    except Exception:
        report['failure']=traceback.format_exc();print(report['failure'],file=sys.stderr,flush=True)
    finally:
        report['finished_at_utc']=datetime.now(timezone.utc).isoformat();save()
    print(json.dumps(dict(all_passed=report['all_passed'],summary=str(summary)),indent=2),flush=True)
    return 0 if report['all_passed'] else 1

if __name__=='__main__':
    raise SystemExit(main())
