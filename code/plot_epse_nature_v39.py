"""Rebuild Project40 figures from audited outputs, not manuscript numbers.

Nature-oriented design uses the MIT-licensed SciencePlots nature/science
styles at pinned GitHub commit b9b16959570bd2fbc9ff5118bacc423c3bddd592.
Vector PDF/SVG retain editable text. Raster exports are natively 1200 dpi.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.text import Text
from matplotlib.transforms import Bbox
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'figures'
RESULT = ROOT/'results/epse_v39'
STYLE = ROOT/'code/vendor/scienceplots_v39'
BLUE, ORANGE, GREY = '#0072B2', '#D55E00', '#6B7075'
INPUTS = [ROOT/'historical/results/calendar_thresholds_v23/strong_baseline_comparisons.csv',
          ROOT/'results/lifecycle_v38/path_verification.json',
          ROOT/'results/lifecycle_v38/extension/structural_results.json',
          ROOT/'historical/results/continuing_pool_v34/periodic_comparisons.csv']
INPUTS = [p if p.exists() else ROOT/'benchmarks/historical'/p.relative_to(ROOT/'historical/results')
          if p.is_relative_to(ROOT/'historical/results') else p for p in INPUTS]
MANIFEST = []

def load_csv(p):
    with p.open(encoding='utf-8',newline='') as f:
        return list(csv.DictReader(f))

def panel(ax, letter, title):
    ax.text(-.02,1.14,letter,transform=ax.transAxes,fontweight='bold',fontsize=8)
    ax.text(.045,1.14,title,transform=ax.transAxes,fontsize=7)
    ax.tick_params(direction='out',length=2.6,width=.5,pad=3)
    ax.spines[['right','top']].set_visible(False)

def global_key(fig):
    handles = [plt.Line2D([],[],color=ORANGE,marker='s',linewidth=1,label='Dedicated'),
               plt.Line2D([],[],color=BLUE,marker='o',linewidth=1,label='Shared')]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.52,.985),
               ncol=2,frameon=False,handlelength=1.8,columnspacing=2.1)

def layout_audit(fig):
    """Conservative rendered bbox and DATA-path audits, including side titles.

    Borders/ticks are not data. Text inside a workflow card is intentional.
    No legend, label or annotation is allowed to obscure a plotted data path.
    """
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    text=[]
    for obj in fig.findobj(match=Text):
        if obj.get_visible() and obj.get_text().strip():
            bb=obj.get_window_extent(renderer)
            if bb.width>0 and bb.height>0:
                text.append((obj,bb))
    problems=[]
    for i,(obj,bb) in enumerate(text):
        if not fig.bbox.contains(bb.x0,bb.y0) or not fig.bbox.contains(bb.x1,bb.y1):
            problems.append({'kind':'clipped_text','text':obj.get_text()})
        for obj2,bb2 in text[i+1:]:
            if bb.overlaps(bb2):
                overlap=Bbox.intersection(bb,bb2)
                if overlap and overlap.width>0.6 and overlap.height>0.6:
                    problems.append({'kind':'text_overlap','texts':[obj.get_text(),obj2.get_text()]})
        for ax in fig.axes:
            for ln in ax.lines:
                if not ln.get_visible(): continue
                path=ln.get_transform().transform_path(ln.get_path())
                if path.intersects_bbox(bb.padded(1.2),filled=False):
                    problems.append({'kind':'text_data_line_overlap','text':obj.get_text()})
            for patch in ax.patches:
                if patch.get_gid()=='workflow-card': continue
                path=patch.get_transform().transform_path(patch.get_path())
                # Compound arrow paths contain disjoint shaft/head polygons;
                # treating their union as one filled path creates false hits.
                if isinstance(patch,FancyArrowPatch):
                    from matplotlib.path import Path as MplPath
                    hit=any(MplPath(poly).intersects_bbox(bb.padded(.8),filled=False)
                            for poly in path.to_polygons(closed_only=False))
                else:
                    hit=path.intersects_bbox(bb.padded(.8),filled=True)
                if hit:
                    problems.append({'kind':'text_graphic_overlap','text':obj.get_text()})
    return {'visible_text_objects':len(text),'findings':problems,
            'zero_detected_overlap_or_clipping':not problems}

def save(fig,name,data):
    audit=layout_audit(fig)
    if audit['findings']:
        print(json.dumps({name:audit},ensure_ascii=False))
        raise RuntimeError('Figure layout audit failed')
    files=[]
    for ext in ['pdf','svg','png','tiff']:
        dest=OUT/f'{name}.{ext}'
        opts={'dpi':1200} if ext in ['png','tiff'] else {}
        if ext=='tiff': opts['pil_kwargs']={'compression':'tiff_lzw'}
        fig.savefig(dest,facecolor='white',**opts)
        files.append({'file':dest.name,'sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),
                      'raster_dpi':1200 if ext in ['png','tiff'] else None})
    fig.savefig(RESULT/f'{name}_preview.png',dpi=180,facecolor='white')
    MANIFEST.append({'name':name,'width_mm':fig.get_size_inches()[0]*25.4,
                     'height_mm':fig.get_size_inches()[1]*25.4,'data':data,
                     'audit':audit,'outputs':files})
    plt.close(fig)

def main():
    OUT.mkdir(parents=True,exist_ok=True); RESULT.mkdir(parents=True,exist_ok=True)
    plt.style.use([str(STYLE/n) for n in ['science.mplstyle','nature.mplstyle','no-latex.mplstyle']])
    plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
        'font.size':7,'axes.labelsize':7,'xtick.labelsize':7,'ytick.labelsize':7,
        'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.grid':False,
        'axes.spines.top':False,'axes.spines.right':False,'xtick.top':False,'ytick.right':False,
        'text.usetex':False,'savefig.bbox':None,'figure.dpi':120})
    # Recompute every lifetime point through the independent exact optimizer.
    from lifecycle_v38_path_verify import assignment_optimum, job_key
    jobs=[dict(id=i,owner=i%3,start=i,ready=i+1,lost=False) for i in range(9)]
    lifetimes=list(range(1,10)); records=[]
    for L in lifetimes:
        row={'life':L}
        for mode in ('dedicated','pooled'):
            row[mode]=assignment_optimum(job_key(jobs),5,L,mode)[0]
            m=min(5,1+(L-1)//(3 if mode=='dedicated' else 1))
            predicted=sum(math.ceil(3/m) for _ in range(3)) if mode=='dedicated' else math.ceil(9/m)
            assert row[mode]==predicted
        records.append(row)
    fig=plt.figure(figsize=(7.0,2.8))
    a=fig.add_axes([.075,.22,.43,.49]); b=fig.add_axes([.63,.22,.33,.49])
    for mode,c,mark in [('dedicated',ORANGE,'s'),('pooled',BLUE,'o')]:
        a.step(lifetimes,[r[mode] for r in records],where='post',color=c,linewidth=1)
        a.plot(lifetimes,[r[mode] for r in records],linestyle='none',marker=mark,color=c,markersize=3)
    a.set(xlabel='Ready-end asset life (planning intervals)',ylabel='New assets (count)',
          xlim=(.7,9.3),ylim=(0,10),xticks=[1,3,5,7,9],yticks=[0,3,6,9])
    panel(a,'a','Lifetime limits attainable service')
    # Long-unavailability counterexample independently checked in v38.
    values=[2,2,3]
    b.bar(np.arange(3),values,color=[GREY,GREY,BLUE],width=.55)
    b.set(xticks=np.arange(3),xticklabels=['Peak\noccupied','Global\ncycle rule','Exact\nnew entries'],
          ylabel='Assets (count)',ylim=(0,4),yticks=[0,1,2,3,4])
    panel(b,'b','Stock is not cumulative manufacture')
    global_key(fig)
    save(fig,'figure1_lifecycle',{'equal_spacing_exact':records,
         'conditions':{'jobs':9,'owners':3,'K':5,'delta':1,'ready_delay':1,'synthetic':True},
         'seven_job_long_unavailability':{'peak_occupied':2,'global_division_lower_bound':2,'exact_births':3}})

    data={r['scenario']:r for r in load_csv(INPUTS[0])}
    names=['flat_tau1','staggered_tau1','staggered_tau3']
    labels=['Flat demand\n1-week return','Staggered demand\n1-week return','Staggered demand\n3-week return']
    fig=plt.figure(figsize=(7.0,3.3))
    a=fig.add_axes([.22,.20,.28,.54]); b=fig.add_axes([.67,.20,.28,.54])
    y=np.arange(3)
    for ax,suffix,divisor in [(a,'initial_assets',1000),(b,'maximum_upper_impact_ratio',1)]:
        for j,(prefix,c) in enumerate([('dedicated',ORANGE),('pooled',BLUE)]):
            ax.barh(y+(j-.5)*.26,[float(data[n][prefix+'_'+suffix])/divisor for n in names],
                    height=.22,color=c)
        ax.set_yticks(y,labels if ax==a else ['']*3); ax.invert_yaxis()
    a.set(xlabel='Opening stock (thousand crates)',xlim=(0,225),xticks=[0,70,140,210])
    b.set(xlabel='Maximum upper source-impact ratio',xlim=(0,1.4),xticks=[0,.5,1,1.3])
    b.axvline(1,color=GREY,linestyle='--',linewidth=.7)
    panel(a,'a','Calendar reserve'); panel(b,'b','Environmental passage')
    global_key(fig)
    save(fig,'figure2_calendar',{'historical_scenario_rows':[data[n] for n in names],
          'scope':'Historical source-response illustration; declared demand, primary-packaging geometry proxy; not observed EU crate operation'})

    age=json.loads(INPUTS[1].read_text())['age_pair_analysis']['age_only_profiles'][0]
    assert age['fixed_prior']==0
    fig=plt.figure(figsize=(7.0,2.95))
    a=fig.add_axes([.08,.23,.40,.48]); b=fig.add_axes([.62,.23,.34,.48])
    # Periodic points use only the continuing-pool boundary, not startup stock.
    delayed=[]
    periodic={r['case']:r for r in load_csv(INPUTS[3])}
    for t in [1,3,5]:
        row=periodic[f'main_uniform_c0_tau{t}']
        d=float(row['dedicated_fleet']); p=float(row['pooled_fleet'])
        delayed.append({'return_weeks':t,'dedicated':d,'shared':p,'saving_percent':100*(d-p)/d})
    a.plot([r['return_weeks'] for r in delayed],[r['saving_percent'] for r in delayed],
           color=GREY,marker='o',linestyle='none')
    a.set(xlabel='Return delay (weeks)',ylabel='Standing-fleet reduction with sharing (%)',
          xlim=(.5,5.5),ylim=(0,40),xticks=[1,3,5],yticks=[0,10,20,30,40])
    panel(a,'a','Reserve saving in continuing pools')
    profile=age['integer_age_profile']
    for mode,c,mark in [('dedicated',ORANGE,'s'),('pooled',BLUE,'o')]:
        b.step([r['age'] for r in profile],[r[mode] for r in profile],where='post',color=c)
        b.plot([r['age'] for r in profile],[r[mode] for r in profile],color=c,marker=mark,linestyle='none')
    b.set(xlabel='Initial asset age (planning intervals)',ylabel='New assets (count)',
          xlim=(-.3,8.3),ylim=(0,7),xticks=[0,2,4,6,8],yticks=[0,2,4,6])
    panel(b,'b','Ageing changes the sharing benefit')
    global_key(fig)
    save(fig,'figure3_delay_age',{'periodic_delay':delayed,'synthetic_age_only':age,
         'warning':'Panels have different scale, purpose and data status; no cohort probability or combined EU estimate'})

    fig=plt.figure(figsize=(7.0,1.65)); ax=fig.add_axes([.01,.08,.98,.83]);ax.axis('off')
    ax.set(xlim=(0,5),ylim=(0,1))
    headings=['Asset evidence','Calendar test','Fleet requirement','Environmental flows','Operating decision']
    lines=['Age and prior uses\nReady and exit dates','Serve each job\nRespect expiry','New introductions\nConcurrent reserve','Manufacture and wash\nEmpty transport','Compare matched service\nRetain failure bounds']
    for i,(heading,line) in enumerate(zip(headings,lines)):
        card=Rectangle((i+.02,.30),.86,.58,facecolor='#F4F6F7',edgecolor='#C2C7CA',linewidth=.6)
        card.set_gid('workflow-card');ax.add_patch(card)
        ax.text(i+.45,.73,heading,ha='center',va='center',fontsize=6.8,fontweight='bold')
        ax.text(i+.45,.48,line,ha='center',va='center',fontsize=6.3,linespacing=1.6)
        if i<4:
            arr=FancyArrowPatch((i+.90,.59),(i+1+.005,.59),arrowstyle='-|>',mutation_scale=6,
                               linewidth=.6,color=GREY);ax.add_patch(arr)
    ax.text(2.45,.09,'Unidentified survival is reported as bounded scenarios, not a fitted probability.',
            ha='center',va='center',fontsize=6.6)
    save(fig,'figure4_workflow',{'scope':'Engineering implementation sequence; no measured data encoded'})

    # Sized separately for the journal's graphical-abstract dimensions.
    fig=plt.figure(figsize=(60/25.4,50/25.4)); ax=fig.add_axes([.02,.02,.96,.96]);ax.axis('off')
    ax.set(xlim=(0,1),ylim=(0,1))
    for y,label,sub in [(.84,'Asset lifetime','Age and remaining uses'),(.53,'Return calendar','Availability and replacement'),(.22,'Environmental benefit','Matched service and full flows')]:
        card=Rectangle((.05,y-.10),.90,.21,facecolor='#F4F6F7',edgecolor='#C2C7CA',linewidth=.6)
        card.set_gid('workflow-card');ax.add_patch(card)
        ax.text(.5,y+.035,label,ha='center',va='center',fontsize=7,fontweight='bold')
        ax.text(.5,y-.045,sub,ha='center',va='center',fontsize=6.2)
    for y in [.735,.425]:
        ax.add_patch(FancyArrowPatch((.5,y-.012),(.5,y-.086),arrowstyle='-|>',mutation_scale=7,
                                    color=BLUE,linewidth=.8))
    save(fig,'graphical_abstract',{'width_mm':60,'height_mm':50,'scope':'Conceptual sequence only'})
    for p in INPUTS:
        assert p.exists()
    provenance={'source_commit':'b9b16959570bd2fbc9ff5118bacc423c3bddd592',
        'source_repository':'https://github.com/garrettj403/SciencePlots','license':'MIT',
        'nature_guide':'https://research-figure-guide.nature.com/figures/preparing-figures-our-specifications/',
        'inputs_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in INPUTS},
        'native_raster_export_dpi':1200,'vectors_have_no_dpi':True,
        'text_is_editable_in_svg_and_embedded_truetype_in_pdf':True,'figures':MANIFEST}
    (RESULT/'figure_provenance_v39.json').write_text(json.dumps(provenance,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({'figures':len(MANIFEST),'layout_audits_pass':True,'exports':len(MANIFEST)*4}))

if __name__=='__main__':main()
