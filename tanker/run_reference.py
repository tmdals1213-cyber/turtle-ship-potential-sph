"""Run qualified reference inputs; preserve existing output and continue on failure."""
from pathlib import Path
import json,traceback
from prepare_reference import prepare_all
from run_coupled import run
from analyze_coupled import analyze,fluid_qa
from extract_surface_loads import extract
from make_summary import summarize
ROOT=Path(__file__).resolve().parent

def main():
    records=[]
    for folder in prepare_all():
        item={'case':folder.name}
        try:
            item['run']=run(folder.name)
            if item['run']['status']=='complete':
                item['motion']=analyze(folder.name);item['fluid_qa']=fluid_qa(folder.name)
                if item['fluid_qa']['passed']:item['loads']=extract(folder.name)
            item['qualified']=bool(item.get('fluid_qa',{}).get('passed') and item.get('loads',{}).get('qualified_for_load_visualization'))
        except Exception as error:
            item.update(qualified=False,error=str(error));traceback.print_exc()
        records.append(item)
        report={'input':{'braking_g':.4,'lateral_g':.69,'lateral_start_s':1.5,'lateral_release_s':4.4,'lateral_end_s':5},
                'sound_speed_m_s':220,'cfl_number':.2,'duration_s':8,'cases':records,
                'policy':'No required improvement outcome. Preserve failed runs, continue remaining cases.'}
        (ROOT/'sph/verified_comparison.json').write_text(json.dumps(report,indent=2))
    if all(c['qualified'] for c in records):summarize()
    else:raise RuntimeError('Some cases failed qualification; inspect reports before rendering.')

if __name__=='__main__':main()
