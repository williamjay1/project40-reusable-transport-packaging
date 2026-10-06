"""Acquire/decode public non-household energy-price contexts; never overwrite raw.

Prices are marginal-energy stress contexts, not a rebase of the historical
Accorsi bundled washing tariff and not observed RPC contract prices.
"""
import argparse
import csv
import hashlib
import json
import os
import stat
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen
import numpy as np

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--electricity-raw',type=Path,required=True)
    ap.add_argument('--gas-raw',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    a=ap.parse_args(); a.output_dir.mkdir(parents=True,exist_ok=True)
    gas_url=('https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/'
             'nrg_pc_203?format=JSON&lang=en&geo=EU27_2020&sinceTimePeriod=2024')
    if not a.gas_raw.exists():
        with urlopen(Request(gas_url,headers={'User-Agent':'Project40-public-research/1.0'}),timeout=60) as r:
            data=r.read()
        obj=json.loads(data)
        assert obj.get('class')=='dataset' and obj.get('source')=='ESTAT'
        a.gas_raw.parent.mkdir(parents=True,exist_ok=True)
        with a.gas_raw.open('xb') as f: f.write(data)
        os.chmod(a.gas_raw,stat.S_IREAD)
    cells=[]; sources=[]; excluded=[]
    for dataset,path,url in [('nrg_pc_205',a.electricity_raw,
        'https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/nrg_pc_205?format=JSON&lang=en&geo=EU27_2020&sinceTimePeriod=2024'),
        ('nrg_pc_203',a.gas_raw,gas_url)]:
        raw=path.read_bytes(); obj=json.loads(raw)
        # Later retrieval is allowed; source snapshot itself must predate cutoff.
        if obj['updated'][:10]>'2026-08-29':
            excluded.append(dict(dataset=dataset,source_updated=obj['updated'],source_url=url,
                raw_file=str(path),raw_sha256=hashlib.sha256(raw).hexdigest(),
                reason='Later revision cannot enter frozen scientific evidence; preserved read-only, not used'))
            continue
        codes=[list(obj['dimension'][d]['category']['index']) for d in obj['id']]
        for index,value in obj['value'].items():
            if value is None: continue
            coordinates=np.unravel_index(int(index),obj['size'])
            row=dict(zip(obj['id'],[codes[k][v] for k,v in enumerate(coordinates)]))
            if (row.get('geo')!='EU27_2020' or row.get('currency')!='EUR'
                or row.get('tax')!='X_VAT' or row.get('unit')!='KWH'
                or row.get('time') not in ['2024-S1','2024-S2','2025-S1','2025-S2']
                or row.get('nrg_cons','').startswith('TOT')): continue
            label=obj['dimension']['nrg_cons']['category']['label'][row['nrg_cons']]
            cells.append(dict(dataset=dataset,band=row['nrg_cons'],band_label=label,
                semester=row['time'],tax='X_VAT',unit='EUR/kWh',value=float(value)))
        sources.append(dict(dataset=dataset,source_url=url,source_updated=obj['updated'],
            raw_sha256=hashlib.sha256(raw).hexdigest(),raw_file=str(path),
            access='public Eurostat API, no account',licence='Eurostat/European Commission reuse policy; attribution required',
            spatial_coverage='EU27_2020 aggregate; no non-EU countries',
            interpretation='Public marginal-input stress context, not an observed wash-site tariff'))
    assert len([x for x in cells if x['dataset']=='nrg_pc_205'])==28
    with (a.output_dir/'industrial_energy_contexts.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(cells[0]));w.writeheader();w.writerows(cells)
    ranges={d:dict(lower=min(c['value'] for c in cells if c['dataset']==d),
                   upper=max(c['value'] for c in cells if c['dataset']==d),
                   cells=sum(c['dataset']==d for c in cells)) for d in sorted({c['dataset'] for c in cells})}
    meta=dict(acquired_at_utc=datetime.now(timezone.utc).isoformat(),sources=sources,excluded_sources=excluded,ranges=ranges,
        tax_convention='X_VAT excludes VAT and other recoverable taxes/levies; retains non-recoverable levies',
        sampling='No probability distribution assigned to price cells; extrema are context sensitivities',
        accounting='Only extra energy quantities above source cleaning are charged at these prices; historical 0.07 EUR baseline washing remains bundled and is not double charged',
        heat_mapping='Source 0.322 MJ/crate is delivered conventional-boiler heat; effective purchased-heat price p_H is a declared coordinate, not an admitted gas observation. Gas-to-heat conversion is not silently assumed.',
        cutoff='All admitted snapshots updated no later than 2026-08-29; no current API revision silently replaces archived electricity',
        selected_rows=len(cells))
    (a.output_dir/'financial_input_audit.json').write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(meta,indent=2))

if __name__=='__main__': main()
