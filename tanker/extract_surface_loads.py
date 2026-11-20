"""Separate saved boundary-particle loads with official ComputeForces.

Source binary and license provenance: computeforces_source.json. No solver
rerun, pressure-traction approximation, or invented per-baffle forces is used.
Loads are recomputed at saved snapshots, then checked against native body loads.
"""
from pathlib import Path
import argparse,json,subprocess
import numpy as np
import pyvista as pv
from scipy.signal import savgol_filter
from analyze_coupled import native
from run_braking import rows

ROOT=Path(__file__).resolve().parent
from tanker_paths import BIN
EXE=BIN/'ComputeForces_win64.exe'
GROUP_NAMES=['Tank wall','Rear transverse baffle','Front transverse baffle','Longitudinal baffle']


def extract(name):
    f=ROOT/'sph'/name;out=f/'surface_loads';out.mkdir(exist_ok=True)
    assert json.loads((f/'fluid_boundary_qa.json').read_text())['passed']
    setup=json.loads((f/'setup.json').read_text());pose=np.load(f/'body_pose.npz')
    rr=rows(f/'output/RunPARTs.csv');last=int(rr[-1]['Part'])
    if not (out/f'BoundaryLoad_{last:04d}.vtk').exists():
        command=[str(EXE),'-dirdata',str(f/'output/data'),'-onlymk:12','-threads:4','-viscoauto',
                 '-savecsv',str(out/'total'),'-savevtk',str(out/'BoundaryLoad')]
        with (out/'extract.log').open('w') as log:
            p=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
        if p.returncode:raise RuntimeError('Official ComputeForces failed')
    m=pv.read(out/'BoundaryLoad_0000.vtk');order=np.argsort(m['Id']);ids=m['Id'][order];x0=m.points[order].astype(float)
    groups=np.zeros(len(ids),np.int8);half=setup.get('numerical_plate_thickness_m',.08)/2
    inside=(abs(x0[:,0])<3.-1e-5)&(np.linalg.norm(x0[:,1:],axis=1)<1.-1e-5)
    if setup['layout']=='D_cross_baffles':groups[inside&(abs(x0[:,1])<half+1e-5)]=3
    if setup['layout']!='A_unpartitioned':
        for j,cx in enumerate([-1.,1.]):groups[inside&(abs(x0[:,0]-cx)<half+1e-5)]=j+1
    forces=[];moments=[];bodyforces=[];bodymoments=[];percentiles=[]
    for i,row in enumerate(rr):
        m=pv.read(out/f'BoundaryLoad_{int(row["Part"]):04d}.vtk');order=np.argsort(m['Id'])
        assert np.array_equal(m['Id'][order],ids)
        ff=m['Force'][order].astype(float);xyz=m.points[order].astype(float)
        mm=np.cross(xyz-pose['center'][i],ff)
        fg=np.array([ff[groups==j].sum(0) for j in range(4)])
        mg=np.array([mm[groups==j].sum(0) for j in range(4)])
        forces.append(fg);moments.append(mg)
        bodyforces.append(fg@pose['rotation'][i]);bodymoments.append(mg@pose['rotation'][i])
        percentiles.append(np.percentile(np.linalg.norm(ff,axis=1),[95,99,99.5,100]))
    force=np.asarray(forces);moment=np.asarray(moments);fb=np.asarray(bodyforces);mb=np.asarray(bodymoments)
    data,meta=native(f/'output/data/PartFloatInfo.ibi4');assert len(data['time_s'])==len(rr)
    smooth=lambda a:savgol_filter(a,9,3,axis=0)
    use=(pose['time_s']>.5)&(pose['time_s']<7.75)
    def compare(a,b):
        aa=smooth(a);bb=smooth(b);scale=np.sqrt(np.mean(bb[use]**2))
        return dict(filtered_nrmse=float(np.sqrt(np.mean((aa[use]-bb[use])**2))/max(1.,scale)),
                    max_filtered_absolute_error=float(abs(aa[use]-bb[use]).max()))
    fc=compare(force.sum(1),data['fluforcelin']);mc=compare(moment.sum(1),data['fluforceang'])
    np.savez_compressed(out/'loads.npz',time_s=pose['time_s'],idp=ids,initial_points=x0,group=groups,
                        force_world_N=force,moment_world_Nm=moment,force_body_N=fb,moment_body_Nm=mb,
                        force_body_filtered_N=smooth(fb),moment_body_filtered_Nm=smooth(mb),
                        per_particle_force_percentiles_N=np.asarray(percentiles))
    header=['time_s']
    columns=[pose['time_s']]
    for j in range(4):
        for label,a in [('F',fb),('M',mb)]:
            for k,axis in enumerate('xyz'):header.append(f'group{j}_{label}{axis}_body');columns.append(a[:,j,k])
    np.savetxt(out/'group_load_history.csv',np.column_stack(columns),delimiter=',',header=','.join(header),comments='')
    report=dict(case=name,groups=GROUP_NAMES,group_counts=np.bincount(groups,minlength=4).tolist(),
                source='Official ComputeForces v5.4.266.02, per-boundary-particle Force vector from saved BI4 states.',
                units='Force N; moment N*m about instantaneous dry-vehicle CG. Body axes rotate with truck.',
                partition='Joint cells at transverse/longitudinal intersections assigned to transverse plates; groups disjoint and sum exactly to whole tank.',
                filtering='9 saved 0.05s states, cubic Savitzky-Golay; 0.4s window for displayed integrated loads. Local impact peaks are not validation targets.',
                local_scalar='Force magnitude in N per boundary particle at DP=40mm. This is a local load visualization, not pressure or grid-independent stress.',
                native_force_comparison=fc,native_moment_comparison=mc,
                qualified_for_load_visualization=bool(fc['filtered_nrmse']<.1 and mc['filtered_nrmse']<.1))
    (out/'analysis.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('cases',nargs='+');a=p.parse_args()
    for name in a.cases:extract(name)
