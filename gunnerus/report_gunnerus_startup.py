"""Report actual startup motion and applied/recovered forces in model units."""
from pathlib import Path
import csv,json,xml.etree.ElementTree as ET
import numpy as np
import pyvista as pv
from scipy.spatial.transform import Rotation
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_propulsion import analyze
from sph_segments import part_rows,vtk_path,output_dirs

ROOT=Path(__file__).resolve().parent

def report(name,motion_stride=4):
    f=ROOT/'sph'/name;s=json.loads((f/'setup.json').read_text());p=s['physical_model']
    state=json.loads((f/'running.json').read_text())
    if state.get('running') or state.get('returncode')!=0 or state.get('stop_reason'):
        raise ValueError('Require a completed, unstopped case.')
    rows=part_rows(f)
    times={int(r['Part']):float(r['TimeStep [s]']) for r in rows}
    assert max(times.values())>=s['duration_s']-1e-6
    assert 'Finished execution (code=0)' in (output_dirs(f)[-1]/'Run.out').read_text(errors='replace')
    analyze(name);force=json.loads((f/'force_history.json').read_text())
    paths=[vtk_path(row) for row in rows];first=pv.read(paths[0])
    mask=first['Type']==2;order=np.argsort(first['Idp'][mask]);ids=first['Idp'][mask][order]
    initial=first.points[mask][order].astype(float);mean0=initial.mean(0);cg0=np.array(p['center_of_mass_m'])
    selected=paths[::motion_stride]
    if selected[-1]!=paths[-1]:selected.append(paths[-1])
    motion=[]
    for path in selected:
        d=pv.read(path);mask=d['Type']==2;order=np.argsort(d['Idp'][mask]);assert np.array_equal(d['Idp'][mask][order],ids)
        current=d.points[mask][order].astype(float);mean1=current.mean(0);u,_,vt=np.linalg.svd((initial-mean0).T@(current-mean1));rot=vt.T@u.T
        if np.linalg.det(rot)<0:vt[-1]*=-1;rot=vt.T@u.T
        cg=rot@(cg0-mean0)+mean1;fluid=d['Type']==3;rho=d['Rho'][fluid];speed=np.linalg.norm(d['Vel'][fluid],axis=1)
        assert np.isfinite(rho).all() and np.isfinite(speed).all()
        motion.append(dict(time_s=times[int(path.stem.split('_')[-1])],cg_m=cg.tolist(),
            rotation_body_to_world=rot.tolist(),
            heave_m=float(cg[2]-cg0[2]),surge_m=float(cg[0]-cg0[0]),
            pitch_deg=float(Rotation.from_matrix(rot).as_euler('xyz',degrees=True)[1]),
            rho_min=float(rho.min()),rho_max=float(rho.max()),max_fluid_speed_m_s=float(speed.max())))
    number=lambda r,k:float(r[k].replace(',',''))
    total_time=sum(float(q['elapsed_s']) for q in json.loads((f/'run.json').read_text()))
    lastspeed=(motion[-1]['cg_m'][0]-motion[-2]['cg_m'][0])/(motion[-1]['time_s']-motion[-2]['time_s'])
    summary=dict(case=name,status='startup_complete_not_equilibrium_or_resistance_validation',
        model_duration_s=max(times.values()),wall_s=total_time,total_particles=int(number(rows[0],'NpSim')) if 'NpSim' in rows[0] else json.loads((f/'initial_audit.json').read_text())['total_particles'],
        particle_losses=int(sum(number(r,'NpOut') for r in rows)),dt_min_clamp_events=int(sum(number(r,'DTsMin') for r in rows)),
        min_dt_s=min(number(r,'DtMin [s]') for r in rows if number(r,'DtMin [s]')>0),
        peak_solver_gpu_mib=max(number(r,'MemGPU [MiB]') for r in rows),
        final_heave_m=motion[-1]['heave_m'],final_surge_m=motion[-1]['surge_m'],final_pitch_deg=motion[-1]['pitch_deg'],
        final_interval_mean_surge_speed_m_s=lastspeed,last_fluid_force_N=force['last_fluid_force_N'],
        last_applied_force_N=force['last_external_force_N'],
        sampled_density_range_kg_m3=[min(q['rho_min'] for q in motion),max(q['rho_max'] for q in motion)],
        sampled_fields_finite=True,limitations=s['limitations'])
    (f/'startup_summary.json').write_text(json.dumps(summary,indent=2));(f/'startup_motion.json').write_text(json.dumps(motion,indent=2))
    return summary,motion,force

if __name__=='__main__':
    names=['gunnerus_float_050_2s','gunnerus_thrust_050_2s'];out=[report(n) for n in names]
    result=dict(status='two_startup_runs_complete',scale_ratio=10,cases=[q[0] for q in out],
        interpretation='Prescribed force is world-fixed at CG. No fluid actuator disk, resolved propeller, incident wave or resistance validation.')
    (ROOT/'results/gunnerus_startup_summary.json').write_text(json.dumps(result,indent=2))
    plt.style.use('dark_background');fig,ax=plt.subplots(2,2,figsize=(14,8),facecolor='#080e18')
    for (summary,motion,force),color,label in zip(out,['#59c6dc','#f1b45b'],['Float','Forward + assumed thrust']):
        t=[q['time_s'] for q in motion];cg=np.array([q['cg_m'] for q in motion]);ft=np.array(force['time_s']);ff=np.array(force['recovered_fluid_force_N']);ef=np.array(force['assumed_external_force_N'])
        ax[0,0].plot(t,cg[:,0]-cg[0,0],color=color,label=label)
        ax[0,1].plot(t,1000*(cg[:,2]-cg[0,2]),color=color)
        ax[1,0].plot(ft,ff[:,0],color=color,label=label+' | fluid')
        if summary['last_applied_force_N'][0]:ax[1,0].plot(ft,ef[:,0],color='white',ls='--',label='Applied force')
        ax[1,1].plot(t,[q['pitch_deg'] for q in motion],color=color)
    for a,title,unit in zip(ax.flat,['Forward displacement','Heave','Longitudinal load','Pitch'],['m','mm','N','deg']):
        a.set(title=title,xlabel='Model time [s]',ylabel=unit);a.grid(alpha=.18)
    ax[0,0].legend();ax[1,0].legend(fontsize=8)
    fig.suptitle('GUNNERUS | FIRST GPU STARTUP',fontsize=22)
    fig.text(.5,.012,'1:10 model units | Calm water | 6 DOF | Assumed external force, no propeller wake',ha='center',color='#b1cad8')
    fig.tight_layout(rect=[0,.035,1,.95]);fig.savefig(ROOT/'media/gunnerus_startup_comparison.png',dpi=130);plt.close(fig)
    print(json.dumps(result,indent=2))
