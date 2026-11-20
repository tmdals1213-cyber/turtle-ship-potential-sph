"""Postprocess completed wave pilots; optionally follow the existing solver jobs.

Does not launch simulations, modify a running case, or invent missing data.
"""
from pathlib import Path
import argparse, json, time, traceback
import numpy as np
from scipy.spatial.transform import Rotation, Slerp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from report_gunnerus_startup import report
from read_native_float import read_float
from sph_segments import native_history

ROOT=Path(__file__).resolve().parent
NAMES=['gunnerus_waves_050_10s','gunnerus_calm_050_10s']


def process(name):
    folder=ROOT/'sph'/name
    summary,motion,_=report(name,motion_stride=1)
    data,meta=native_history(folder)
    t=data['time_s']; pt=np.array([q['time_s'] for q in motion])
    orientation=Slerp(pt,Rotation.from_matrix([q['rotation_body_to_world'] for q in motion]))
    rot=orientation(np.clip(t,pt[0],pt[-1]))
    data['force_body_N']=rot.inv().apply(data['fluforcelin'])
    data['moment_body_Nm']=rot.inv().apply(data['fluforceang'])
    data['euler_xyz_deg']=rot.as_euler('xyz',degrees=True)
    data['net_vertical_N']=data['fluforcelin'][:,2]-meta['mass_kg']*9.81
    data['heave_m']=data['center'][:,2]-data['center'][0,2]
    data['speed_body_m_s']=rot.inv().apply(data['fvel'])[:,0]
    np.savez_compressed(folder/'monitors.npz',**data)
    meta.update(orientation_source='Rigid fit to actual hull particles at 0.1 s; SLERP to native load timestamps at the recorded cadence.',
        units='Model units: m, s, N, N m, degrees',positive_pitch='Right-handed +y; positive pitch lowers the bow for x-forward/z-up.',
        not_included=['Calibrated propeller operating point','Actuator-disk wake','Validated slam pressure','Structural bending moment'])
    (folder/'monitors.json').write_text(json.dumps(meta,indent=2))
    use=t>0
    summary.update(native_monitor=meta,
        final_native_speed_body_m_s=float(data['speed_body_m_s'][-1]),
        heave_range_m=[float(data['heave_m'].min()),float(data['heave_m'].max())],
        pitch_range_deg=[float(data['euler_xyz_deg'][:,1].min()),float(data['euler_xyz_deg'][:,1].max())],
        body_pitch_moment_range_Nm=[float(data['moment_body_Nm'][use,1].min()),float(data['moment_body_Nm'][use,1].max())])
    (folder/'wave_pilot_summary.json').write_text(json.dumps(summary,indent=2))
    return summary


def figures(summaries):
    plt.style.use('dark_background')
    fig,axes=plt.subplots(3,2,figsize=(14,10),facecolor='#080e18')
    for name,color,label in [(NAMES[0],'#64cbe7','Head waves'),(NAMES[1],'#e3b75e','Calm water')]:
        if name not in summaries: continue
        with np.load(ROOT/'sph'/name/'monitors.npz') as q:
            use=q['time_s']>0; t=q['time_s'][use]
            values=[q['speed_body_m_s'],q['force_body_N'][:,0],q['heave_m']*1000,
                    q['net_vertical_N'],q['euler_xyz_deg'][:,1],q['moment_body_Nm'][:,1]]
            for ax,v in zip(axes.flat,values): ax.plot(t,v[use],color=color,label=label,lw=1.2)
    for ax,title,unit in zip(axes.flat,['Forward speed','Surge load | body x','Heave',
        'Net vertical load | world z','Pitch | bow-down positive','Pitch moment | body y'],['m/s','N','mm','N','deg','N m']):
        ax.set(title=title,xlabel='Model time [s]',ylabel=unit); ax.grid(alpha=.15)
    axes[0,0].legend(); fig.suptitle('GUNNERUS | MATCHED INPUT THRUST',fontsize=22)
    fig.text(.5,.012,'1:10 model units | Free 6 DOF | Assumed 30 N force | Transient pilot, no propeller wake',ha='center',color='#b0c3d5')
    fig.tight_layout(rect=[0,.035,1,.95]); fig.savefig(ROOT/'media/gunnerus_head_wave_comparison.png',dpi=140); plt.close(fig)


def wave_figure():
    folder=ROOT/'sph/gunnerus_empty_050_10s'
    plt.style.use('dark_background'); fig,ax=plt.subplots(figsize=(12,5),facecolor='#080e18')
    for x in [10,6,2,0]:
        q=np.genfromtxt(folder/f'output/GaugesSWL_x{x}.csv',names=True,delimiter=';')
        ax.plot(q['time_s'],q['swlz_m'],label=f'x = {x} m',lw=1.4)
    ax.set(title='INCIDENT WAVE | EMPTY-TANK CHECK',xlabel='Model time [s]',ylabel='Measured surface elevation [m]')
    ax.grid(alpha=.18); ax.legend(ncol=4)
    fig.text(.5,.015,'1:10 model | Piston target H = 0.4 m, T = 2.2 s | Ramp and propagation are included',ha='center',color='#b0c3d5')
    fig.tight_layout(rect=[0,.04,1,1]); fig.savefig(ROOT/'media/gunnerus_incident_wave_check.png',dpi=140); plt.close(fig)


def main(watch):
    wave_figure(); summaries={}; failed={}
    while True:
        states={}
        for name in NAMES:
            folder=ROOT/'sph'/name
            finished=folder/'wave_pilot_summary.json'
            if finished.exists(): summaries[name]=json.loads(finished.read_text())
            try:
                state=json.loads((folder/'running.json').read_text()) if (folder/'running.json').exists() else {}
            except (json.JSONDecodeError,OSError):
                # The existing solver wrapper refreshes this small live file.
                state={'running':True,'note':'Live status being refreshed'}
            states[name]=state
            # GenCase completion is not solver completion.
            done=(folder/'output/Run.out').exists() and state.get('running') is False
            if done and name not in summaries and name not in failed:
                try:
                    if state.get('returncode')!=0 or state.get('stop_reason'):
                        raise RuntimeError(str(state))
                    summaries[name]=process(name); figures(summaries)
                    print('Postprocessed '+name,flush=True)
                except Exception:
                    failed[name]=traceback.format_exc(); print(failed[name],flush=True)
        complete=len(summaries)+len(failed)==len(NAMES)
        status=dict(status='complete' if complete and not failed else ('completed_with_errors' if complete else 'running'),
            cases=summaries,errors=failed,job_states=states,
            comparison='Matched external force and initial speed, not maintained equal speed or validated added resistance.',
            next_propulsion='Actuator disk agreed as next model; absent from these runs.')
        target=ROOT/'results/gunnerus_head_wave_status.json'; temp=target.with_suffix('.tmp')
        temp.write_text(json.dumps(status,indent=2)); temp.replace(target)
        if complete or not watch: break
        time.sleep(10)


if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--watch',action='store_true'); main(ap.parse_args().watch)
