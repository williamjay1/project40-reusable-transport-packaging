"""Modern foreground bridge, never a complete LCIA. No empirical demand."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from run_lifecycle_v38 import Job, path_masks, subset_dp

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/lifecycle_v38/foreground_v2"


def main():
    OUT.mkdir(exist_ok=True,parents=True)
    data_path=OUT/"foreground_bridge.json"
    assert not data_path.exists()
    facts=json.loads((ROOT/"results/lifecycle_v38/modern_lci_numeric_facts.json").read_text(encoding="utf-8"))
    mass=next(x for x in facts["physical_inputs"] if x["id"]=="rpc_weight")
    electricity=[x for x in facts["physical_inputs"]
                 if x["stage"]=="birth" and x["system"]=="RPC" and x["unit"]=="kWh"]
    assert mass["value"]==1.2 and len(electricity)==1 and electricity[0]["value"]==0.31
    wash_e=next(x for x in facts["physical_inputs"] if x["id"]=="rpc_wash_electricity_context")
    wash_w=next(x for x in facts["physical_inputs"] if x["id"]=="rpc_wash_water_context")
    assert wash_e["value"]==.29 and wash_w["value"]==1.2
    jobs=[Job(i,0 if i<2 else 1,t,t+1) for i,t in enumerate((0,1,4,5))]
    rows=[]
    for L in (2,101):
        births={mode:subset_dp(path_masks(jobs,4,L,mode),[],4)
                for mode in ("pooled","dedicated")}
        delta=births["dedicated"]-births["pooled"]
        rows.append(dict(calendar_ready_life=L,new_births=births,
                         equal_services=4,fish_kg_per_service=10,
                         final_ready_washes_each_mode=4,
                         gross_new_crate_mass_kg={m:mass["value"]*v for m,v in births.items()},
                         manufacture_electricity_kWh={m:electricity[0]["value"]*v for m,v in births.items()},
                         conditional_manufacture_plus_wash_electricity_kWh={m:electricity[0]["value"]*v+4*wash_e["value"] for m,v in births.items()},
                         conditional_wash_water_kg_each_mode=4*wash_w["value"],
                         ice_kg_each_mode=8,
                         avoided_gross_new_crate_mass_kg=delta*mass["value"],
                         avoided_manufacture_electricity_kWh=delta*electricity[0]["value"],
                         new_calendar_and_births_are_synthetic=True,
                         not_complete_materials_energy_or_lcia=True,
                         wash_context_assignment_common_and_cancels=True))
    result=dict(all_passed=True,source_doi="10.1016/j.scitotenv.2024.175452",
                source_locator="Table1 gross RPC mass; Table2 RPC manufacture electricity; article p.4",
                source_foreground_years=[2021,2022],
                service_boundary="four unit RPC deliveries, each10kg fish plus2kg ice, final-ready endpoint",
                washing_assignment_condition=wash_e["condition"],
                unknown_terms=["steam", "material/masterbatch recipe", "transport and complete washing technology",
                               "terminal gross burdens and substitution", "food protection and refrigeration"],
                no_single_use_superiority_or_current_EU_effect=True,rows=rows)
    with data_path.open("x",encoding="utf-8") as stream:
        json.dump(result,stream,indent=2);stream.write("\n")

    lives=np.arange(1,10)
    shared=[int(np.ceil(9/min(5,L))) for L in lives]
    dedicated=[3*int(np.ceil(3/min(5,1+(L-1)//3))) for L in lives]
    for L,bp,bd in zip(lives,shared,dedicated):
        equal_jobs=[Job(i,i%3,i,i+1) for i in range(9)]
        assert subset_dp(path_masks(equal_jobs,5,int(L),"pooled"),[],9)==bp
        assert subset_dp(path_masks(equal_jobs,5,int(L),"dedicated"),[],9)==bd
    fig,axes=plt.subplots(1,2,figsize=(7.1,3.05),layout="constrained")
    axes[0].step(lives,dedicated,where="post",label="Dedicated",color="#8b4a34",lw=1.9)
    axes[0].step(lives,shared,where="post",label="Shared",color="#245f8a",lw=1.9)
    axes[0].set(xlabel="Ready-end lifespan (intervals)",ylabel="Minimum new asset entries",
                xticks=[1,3,5,7,9],yticks=[1,3,5,7,9],ylim=(0.5,9.5),xlim=(1,9))
    axes[0].legend(frameon=False,fontsize=8)
    axes[0].set_title("a  Equal-spacing boundary",loc="left",fontsize=10)
    for y,upper,color in ((1,2,"#245f8a"),(0,3,"#8b4a34")):
        axes[1].plot([0,upper],[y,y],color=color,lw=2.5)
        axes[1].scatter([0,upper],[y,y],color=color,zorder=3)
    axes[1].set(xlabel="New entries over conditional state support",
                yticks=[0,1],yticklabels=["Dedicated","Shared"],xticks=[0,1,2,3],
                xlim=(-.2,3.2),ylim=(-.45,1.45))
    axes[1].set_title("b  Sharp support ranges",loc="left",fontsize=10)
    for ax in axes:
        ax.spines[["top","right"]].set_visible(False)
        ax.tick_params(labelsize=8)
        ax.xaxis.label.set_size(8.5);ax.yaxis.label.set_size(8.5)
    fig.savefig(OUT/"lifecycle_boundary_v38.pdf")
    fig.savefig(OUT/"lifecycle_boundary_v38.png",dpi=900)
    plt.close(fig)
    print(json.dumps({"all_passed":True,"rows":rows},indent=2))


if __name__=="__main__": main()
