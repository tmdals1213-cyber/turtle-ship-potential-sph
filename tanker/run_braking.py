"""Run audited tanker cases sequentially; preserve outputs and continue after failures."""
from pathlib import Path
import csv,json,os,subprocess,time

ROOT=Path(__file__).resolve().parent
from tanker_paths import BIN

def rows(path):
    if not path.exists():return []
    try:
        with path.open() as f:return [r for r in csv.DictReader(f,delimiter=';') if r.get('Part','').isdigit()]
    except (OSError,ValueError):return []

def run(name):
    folder=ROOT/'sph'/name
    if (folder/'run.json').exists():return json.loads((folder/'run.json').read_text())
    if (folder/'output').exists():raise FileExistsError(f'Existing incomplete output: {name}')
    assert json.loads((folder/'initial_audit.json').read_text())['passed']
    setup=json.loads((folder/'setup.json').read_text())
    env=os.environ.copy();env['PATH']=str(BIN)+';'+env.get('PATH','');env['OMP_NUM_THREADS']='4'
    cmd=[str(BIN/'DualSPHysics5.4_win64.exe'),'-gpu:0',str(folder/'generated/Case'),str(folder/'output'),'-dirdataout','data','-sv:binx,info','-svres']
    start=time.perf_counter();reason=None
    with (folder/'solver.log').open('w') as log:
        proc=subprocess.Popen(cmd,cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT)
        while proc.poll() is None:
            rr=rows(folder/'output/RunPARTs.csv');last=rr[-1] if rr else {}
            try:
                loss=sum(float(r['NpOut']) for r in rr)
                if loss>setup['fluid_particles']*.001:reason='More than 0.1 percent fluid particles excluded'
                if len(rr)>=4 and all(float(r['DTsMin'])>10 for r in rr[-3:]):reason='Sustained minimum time-step clamping'
            except (KeyError,ValueError):loss=None
            state=dict(case=name,running=True,pid=proc.pid,wall_s=time.perf_counter()-start,
                       last_saved_time_s=float(last.get('TimeStep [s]',0)),particle_losses=loss,stop_reason=reason)
            (ROOT/'sph/running.json').write_text(json.dumps(state,indent=2))
            if reason:proc.terminate();break
            time.sleep(1)
        code=proc.wait()
    rr=rows(folder/'output/RunPARTs.csv')
    last=float(rr[-1]['TimeStep [s]']) if rr else 0
    native=(folder/'output/Run.out').read_text(errors='replace') if (folder/'output/Run.out').exists() else ''
    complete=code==0 and reason is None and last>=setup['duration_s']-1e-6 and 'Finished execution (code=0)' in native
    report=dict(case=name,status='complete' if complete else 'failed_or_incomplete',returncode=code,stop_reason=reason,
                last_saved_time_s=last,wall_s=time.perf_counter()-start,
                particle_losses=sum(float(r['NpOut']) for r in rr),dt_min_clamps=sum(float(r['DTsMin']) for r in rr),
                min_dt_s=min((float(r['DtMin [s]']) for r in rr if float(r['DtMin [s]'])>0),default=None))
    (folder/'run.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report),flush=True)
    return report

if __name__=='__main__':
    manifest=json.loads((ROOT/'sph/braking_manifest.json').read_text())
    reports=[]
    for entry in manifest['cases']:
        try:reports.append(run(entry['case']))
        except Exception as e:
            report={'case':entry['case'],'status':'runner_error','error':str(e)};reports.append(report);print(json.dumps(report),flush=True)
        (ROOT/'sph/batch_results.json').write_text(json.dumps(reports,indent=2))
    (ROOT/'sph/running.json').write_text(json.dumps({'running':False,'cases':reports},indent=2))
