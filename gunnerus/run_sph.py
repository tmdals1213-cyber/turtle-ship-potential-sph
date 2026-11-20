"""Run the prepared GenCase and optional GPU solver; preserve logs per case."""
from pathlib import Path
import argparse,subprocess,os,json,time,csv,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('case',nargs='?',default='smoke_hull');ap.add_argument('--generate-only',action='store_true');ap.add_argument('--solver-only',action='store_true');ap.add_argument('--guard',action='store_true',help='Stop after sustained dt-min clamping or more than 1 percent particle loss.')
ap.add_argument('--resume-part',type=int);ap.add_argument('--resume-from',default='output');ap.add_argument('--resume-data');ap.add_argument('--output-dir',default='output');a=ap.parse_args()
folder=ROOT/'sph'/a.case
output=(folder/a.output_dir).resolve();assert output.is_relative_to(folder.resolve())
bindir=Path(os.environ['DUALSPHYSICS_BIN']).expanduser().resolve()
env=os.environ.copy();env['PATH']=str(bindir)+';'+env.get('PATH','');env['OMP_NUM_THREADS']='4'
generated=folder/'generated';generated.mkdir(exist_ok=True)
commands=[[str(bindir/'GenCase_win64.exe'),'Case_Def',str(generated/'Case'),'-save:all']]
if not a.generate_only:commands.append([str(bindir/'DualSPHysics5.4_win64.exe'),'-gpu:0',str(generated/'Case'),str(output),'-dirdataout','data','-sv:binx,info,vtk','-svres','-svextraparts:1'])
if not a.generate_only and output.exists():
    raise FileExistsError(f'Preserve existing solver output: {output}')
history=[]
if a.resume_part is not None:
    assert a.solver_only and not a.generate_only and a.resume_part>0
    previous=(folder/a.resume_from).resolve();assert previous.is_relative_to(folder.resolve())
    restart_data=(folder/a.resume_data).resolve() if a.resume_data else previous/'data'
    assert restart_data.is_relative_to(folder.resolve())
    assert not output.exists(), 'Resume into a new directory; never overwrite a saved segment.'
    assert not json.loads((folder/'running.json').read_text()).get('running')
    for name in [f'Part_{a.resume_part:04d}.bi4',f'PartExtra_{a.resume_part:04d}.bi4','Part_Head.ibi4','PartFloatInfo.ibi4']:
        assert (restart_data/name).is_file(), f'Missing restart data: {name}'
    commands[-1]+=[f'-partbegin:{a.resume_part}:{a.resume_part}',str(restart_data)]
    manifest=folder/'output_segments.json'
    segments=json.loads(manifest.read_text()) if manifest.exists() else {'outputs':[a.resume_from]}
    assert a.output_dir not in segments['outputs'];segments['outputs'].append(a.output_dir)
    manifest.write_text(json.dumps(segments,indent=2))
    history=json.loads((folder/'run.json').read_text())
for i,cmd in enumerate(commands):
    if i==0 and a.solver_only:
        assert (generated/'Case.xml').exists(), 'Generate and audit the case first.'
        continue
    start=time.perf_counter();log=folder/('gencase.log' if i==0 else ('solver.log' if a.resume_part is None else f'solver_{a.output_dir}.log'))
    reason=None
    with log.open('w') as f:
        proc=subprocess.Popen(cmd,cwd=folder,env=env,stdout=f,stderr=subprocess.STDOUT)
        while proc.poll() is None:
            live=dict(case=a.case,pid=proc.pid,phase='solver' if i else 'gencase',running=True,elapsed_wall_s=time.perf_counter()-start,output_directory=a.output_dir,resume_part=a.resume_part)
            if i and a.guard:
                parts=output/'RunPARTs.csv'
                if parts.exists():
                    try:
                        rows=[r for r in csv.DictReader(parts.open(),delimiter=';') if r.get('Part','').isdigit()]
                        if rows:
                            last=rows[-1];number=lambda r,k:float(r[k].replace(',',''))
                            losses=sum(number(r,'NpOut') for r in rows)
                            live.update(part=int(last['Part']),physical_time_s=number(last,'TimeStep [s]'),particle_losses=losses)
                            if len(rows)>=4 and all(number(r,'DTsMin')>10 for r in rows[-3:]):reason='Sustained time-step minimum clamping over three outputs.'
                            if len(rows)>=8:
                                reference=max(number(r,'DtMax [s]') for r in rows[1:6])
                                if all(number(r,'DtMax [s]')<.25*reference for r in rows[-3:]):reason='Time step collapsed below one quarter of startup values over three outputs.'
                            if losses>0.01*number(rows[0],'NpfSim'):reason='More than one percent of initial fluid particles excluded cumulatively.'
                    except (ValueError,KeyError,TypeError,OSError):pass # A partially flushed row is retried.
            if reason:
                proc.terminate();live['stop_reason']=reason
            (folder/'running.json').write_text(json.dumps(live,indent=2))
            if reason:break
            time.sleep(2)
        code=proc.wait()
    if i and code==0 and reason is None:
        # Native MinFluidStop can terminate early with exit code zero.
        finished=[r for r in csv.DictReader((output/'RunPARTs.csv').open(),delimiter=';') if r.get('Part','').isdigit()]
        target=float(ET.parse(folder/'Case_Def.xml').getroot().find('.//parameter[@key="TimeMax"]').get('value'))
        interval=float(ET.parse(folder/'Case_Def.xml').getroot().find('.//parameter[@key="TimeOut"]').get('value'))
        if not finished or float(finished[-1]['TimeStep [s]'])<target-interval-1e-5:
            reason='Solver exited before the requested physical duration; inspect native stopping criteria.'
            code=2
    (folder/'running.json').write_text(json.dumps(dict(case=a.case,pid=proc.pid,running=False,returncode=code,stop_reason=reason,output_directory=a.output_dir,resume_part=a.resume_part),indent=2))
    history.append(dict(command=cmd,returncode=code,elapsed_s=time.perf_counter()-start,log=str(log),stop_reason=reason,output_directory=a.output_dir,resume_part=a.resume_part))
    (folder/'run.json').write_text(json.dumps(history,indent=2))
    print(json.dumps(history[-1]),flush=True)
    if code:raise SystemExit(code)
