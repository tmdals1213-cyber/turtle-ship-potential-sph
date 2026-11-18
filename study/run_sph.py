# SPDX-License-Identifier: GPL-3.0-or-later
from project_paths import ASSET_DIR, SOLVER_BIN, FONT_DIR
"""Run the prepared GenCase and optional GPU solver; preserve logs per case."""
from pathlib import Path
import argparse,subprocess,os,json,time,csv,xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('case',nargs='?',default='smoke_hull');ap.add_argument('--generate-only',action='store_true');ap.add_argument('--solver-only',action='store_true');ap.add_argument('--guard',action='store_true',help='Stop after sustained dt-min clamping or more than 1 percent particle loss.');a=ap.parse_args()
folder=ROOT/'sph'/a.case
bindir=SOLVER_BIN
env=os.environ.copy();env['PATH']=str(bindir)+';'+env.get('PATH','');env['OMP_NUM_THREADS']='4'
if not folder.is_dir(): raise SystemExit('Prepare a case first; see docs/REPRODUCE.md.')
if (folder/'output').exists() and not a.generate_only: raise SystemExit('Output already exists; preserve it and choose a new run name.')
for binary in ['GenCase_win64.exe','DualSPHysics5.4_win64.exe']:
    if not (bindir/binary).exists(): raise SystemExit('Missing solver tool: '+str(bindir/binary)+'; set DUALSPHYSICS_BIN.')
generated=folder/'generated';generated.mkdir(exist_ok=True)
commands=[[str(bindir/'GenCase_win64.exe'),'Case_Def',str(generated/'Case'),'-save:all']]
if not a.generate_only:commands.append([str(bindir/'DualSPHysics5.4_win64.exe'),'-gpu:0',str(generated/'Case'),str(folder/'output'),'-dirdataout','data','-sv:binx,info,vtk','-svres'])
history=[]
for i,cmd in enumerate(commands):
    if i==0 and a.solver_only:
        assert (generated/'Case.xml').exists(), 'Generate and audit the case first.'
        continue
    start=time.perf_counter();log=folder/('gencase.log' if i==0 else 'solver.log')
    reason=None
    with log.open('w') as f:
        proc=subprocess.Popen(cmd,cwd=folder,env=env,stdout=f,stderr=subprocess.STDOUT)
        while proc.poll() is None:
            live=dict(case=a.case,pid=proc.pid,phase='solver' if i else 'gencase',running=True,elapsed_wall_s=time.perf_counter()-start)
            if i and a.guard:
                parts=folder/'output/RunPARTs.csv'
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
        finished=[r for r in csv.DictReader((folder/'output/RunPARTs.csv').open(),delimiter=';') if r.get('Part','').isdigit()]
        target=float(ET.parse(folder/'Case_Def.xml').getroot().find('.//parameter[@key="TimeMax"]').get('value'))
        interval=float(ET.parse(folder/'Case_Def.xml').getroot().find('.//parameter[@key="TimeOut"]').get('value'))
        if not finished or float(finished[-1]['TimeStep [s]'])<target-interval-1e-5:
            reason='Solver exited before the requested physical duration; inspect native stopping criteria.'
            code=2
    (folder/'running.json').write_text(json.dumps(dict(case=a.case,pid=proc.pid,running=False,returncode=code,stop_reason=reason),indent=2))
    history.append(dict(command=cmd,returncode=code,elapsed_s=time.perf_counter()-start,log=str(log),stop_reason=reason))
    (folder/'run.json').write_text(json.dumps(history,indent=2))
    print(json.dumps(history[-1]),flush=True)
    if code:raise SystemExit(code)
