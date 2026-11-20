"""Measure body pose from native boundary particles and loads from native output.

Rigid-fit pose avoids integration drift from angular velocity. Wheel clearance
uses the same polygonal wheel geometry as contact. It is a geometric indicator,
not a reconstructed per-tire load. Results describe this rigid truck surrogate.
"""
from pathlib import Path
import argparse, csv, json, subprocess, sys
import numpy as np
import pyvista as pv
from scipy.spatial.transform import Rotation
from run_braking import rows

ROOT = Path(__file__).resolve().parent
from tanker_paths import BIN
from read_native_float import Reader


def native(path):
    r = Reader(path)
    header = r.item()
    assert header['values']['FtCount'] == 1
    series = {}
    while r.pos < len(r.data):
        rec = r.item(); a = rec['arrays']; v = rec['values']
        series.setdefault('time_s', []).append(np.asarray(a.get('TimeStep', [v.get('TimeStep')]), float))
        for key in ['center','fvel','fomega','facelin','faceang','extforcelin','extforceang',
                    'fluforcelin','fluforceang','preacelin','preaceang']:
            series.setdefault(key, []).append(a[key])
    data = {k:np.concatenate(v).astype(float) for k,v in series.items()}
    assert all(np.isfinite(v).all() for v in data.values())
    assert np.all(np.diff(data['time_s']) > 0)
    mass = float(header['arrays']['Mass'][0]); use = data['time_s'] > 0
    residual = mass*data['preacelin'][use]-data['fluforcelin'][use]-data['extforcelin'][use]-[0,0,-9.81*mass]
    error = float(abs(residual).max())
    # Native force arrays are float32; scale the check to the largest force.
    tolerance = max(.02, 2e-6*max(mass*9.81, abs(data['fluforcelin']).max(), abs(data['extforcelin']).max()))
    assert error < tolerance, (error, tolerance)
    return data, dict(mass_kg=mass, max_balance_residual_N=error, balance_tolerance_N=tolerance)


def wheel_points():
    # Same 32-sided cylinders as the native contact mesh.
    angle = np.arange(32)*2*np.pi/32
    circle = np.column_stack((.5*np.cos(angle), np.zeros(32), .5*np.sin(angle)))
    return [np.concatenate([circle+[x,y+d,-1.6] for d in [-.125,.125]])
            for x in [-2.1,-.8,3.8] for y in [-1.,1.]]


def convert(folder, output, label='Body', fluid=False, stride=1):
    dest = folder/('vtks_smoke' if output.name.endswith('smoke') else 'vtks')
    dest.mkdir(exist_ok=True)
    rr = rows(output/'RunPARTs.csv')
    # PartVTK supports selecting saved files; all output is small at this DP.
    if not (dest/f'{label}_{int(rr[-1]["Part"]):04d}.vtk').exists():
        cmd=[str(BIN/'PartVTK_win64.exe'), '-dirdata',str(output/'data'),'-savevtk',str(dest/label),
             '-onlytype:-all,+'+('fluid' if fluid else 'floating'),'-vars:idp,vel,rhop','-threads:4']
        with (folder/f'convert_{label}_{output.name}.log').open('w') as log:
            p=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
        if p.returncode:raise RuntimeError('PartVTK failed')
    return dest, rr


def analyze(name, smoke=False):
    folder=ROOT/'sph'/name; suffix='_smoke' if smoke else ''
    output=folder/('output'+suffix)
    run=json.loads((folder/('run'+suffix+'.json')).read_text())
    assert run['status']=='complete'
    setup=json.loads((folder/'setup.json').read_text())
    data,meta=native(output/'data/PartFloatInfo2.ibi4')
    np.savez_compressed(folder/('native_float'+suffix+'.npz'),**data)
    dest,rr=convert(folder,output)
    first=pv.read(dest/'Body_0000.vtk'); order=np.argsort(first['Idp'])
    # Even subsampling of IDs keeps the fit fast and spreads landmarks over the tank.
    ids=np.asarray(first['Idp'])[order][::19]; initial=first.points[order][::19].astype(float)
    ini_mean=initial.mean(0); ini=initial-ini_mean
    centers=[];matrices=[];euler=[];errors=[];clearance=[];times=[];tilt=[]
    wheels=wheel_points();cg0=np.asarray(setup['dry_cg_m']);road=setup['road_z_m']
    for row in rr:
        mesh=pv.read(dest/f'Body_{int(row["Part"]):04d}.vtk');order=np.argsort(mesh['Idp'])
        assert np.array_equal(mesh['Idp'][order][::19],ids)
        xyz=mesh.points[order][::19].astype(float);mean=xyz.mean(0)
        u,_,vh=np.linalg.svd(ini.T@(xyz-mean));rot=vh.T@u.T
        assert np.linalg.det(rot)>.999
        translation=mean-rot@ini_mean
        error=np.linalg.norm(initial@rot.T+translation-xyz,axis=1).max()
        centers.append(rot@cg0+translation);matrices.append(rot);errors.append(error)
        euler.append(Rotation.from_matrix(rot).as_euler('xyz',degrees=True))
        clearance.append([float((w@rot.T+translation)[:,2].min()-road) for w in wheels])
        tilt.append(float(np.degrees(np.arccos(np.clip(rot[2,2],-1,1)))))
        times.append(float(row['TimeStep [s]']))
    tt=np.asarray(times);cc=np.asarray(centers);ee=np.asarray(euler);ww=np.asarray(clearance);tilt=np.asarray(tilt)
    assert max(errors)<.005, 'Moving particle cloud is no longer a rigid transform'
    # Positive clearance at all three wheels on one side, > 20 mm, is wheel lift.
    side_lift=(ww[:,[0,2,4]].min(1)>.02)|(ww[:,[1,3,5]].min(1)>.02)
    overturned=tilt>=75
    held=np.convolve(overturned.astype(int),np.ones(5,dtype=int),'valid')>=5 if len(tt)>=5 else np.array([False])
    result=dict(case=name, status='Measured native two-way SPH/Chrono result',
                physical_time_s=float(tt[-1]),max_abs_roll_deg=float(abs(ee[:,0]).max()),
                max_abs_pitch_deg=float(abs(ee[:,1]).max()),max_tilt_deg=float(tilt.max()),
                final_tilt_deg=float(tilt[-1]),overturned=bool(held.any()),
                overturn_definition='Body-up tilt >= 75 deg at 5 consecutive saved states (~0.2 s).',
                first_75deg_time_s=float(tt[np.flatnonzero(overturned)[0]]) if overturned.any() else None,
                first_side_lift_time_s=float(tt[np.flatnonzero(side_lift)[0]]) if side_lift.any() else None,
                max_center_horizontal_shift_m=float(np.linalg.norm(cc[:,:2]-cc[0,:2],axis=1).max()),
                min_wheel_clearance_m=float(ww.min()),max_rigid_fit_error_m=float(max(errors)),
                pre_brake_max_tilt_deg=float(tilt[tt<1.5].max()),
                native_force_check=meta,run=run,
                interpretation='Rigid fixed-wheel truck on an accelerating support platform, with frictional contact. No suspension, tire slip curve, steering linkage or road-vehicle calibration.')
    np.savez_compressed(folder/('body_pose'+suffix+'.npz'),time_s=tt,center=cc,rotation=np.asarray(matrices),
                        euler_deg=ee,tilt_deg=tilt,wheel_clearance_m=ww,side_lift=side_lift)
    np.savetxt(folder/('body_history'+suffix+'.csv'),np.column_stack((tt,cc,ee,tilt,ww)),delimiter=',',
               header='time_s,cg_x_m,cg_y_m,cg_z_m,roll_deg,pitch_deg,yaw_deg,tilt_deg,w1_clear_m,w2_clear_m,w3_clear_m,w4_clear_m,w5_clear_m,w6_clear_m',comments='')
    (folder/('coupled_analysis'+suffix+'.json')).write_text(json.dumps(result,indent=2))
    print(json.dumps(result),flush=True)
    return result


def fluid_qa(name):
    folder=ROOT/'sph'/name; setup=json.loads((folder/'setup.json').read_text())
    pose=np.load(folder/'body_pose.npz');dest,rr=convert(folder,folder/'output',label='Fluid',fluid=True)
    origin=None;ids=None;max_x=0.;max_radius=0.;crossed=set();cross_y=set();rho_min=1e9;rho_max=0.;speed=0.
    cg0=np.asarray(setup['dry_cg_m']);half=setup.get('numerical_plate_thickness_m',.08)/2
    for i,row in enumerate(rr):
        mesh=pv.read(dest/f'Fluid_{int(row["Part"]):04d}.vtk');order=np.argsort(mesh['Idp']);current=mesh['Idp'][order]
        if ids is None:ids=current.copy()
        assert np.array_equal(ids,current), 'Fluid count/IDs changed'
        points=mesh.points[order].astype(float)
        local=(points-pose['center'][i])@pose['rotation'][i]+cg0
        if origin is None:origin=local.copy()
        max_x=max(max_x,float(abs(local[:,0]).max()));max_radius=max(max_radius,float(np.linalg.norm(local[:,1:],axis=1).max()))
        rho=np.asarray(mesh['Rhop']);rho_min=min(rho_min,float(rho.min()));rho_max=max(rho_max,float(rho.max()))
        speed=max(speed,float(np.linalg.norm(mesh['Vel'],axis=1).max()))
        if setup['layout']=='C_separate_compartments':
            changed=np.searchsorted([-1,1],origin[:,0])!=np.searchsorted([-1,1],local[:,0])
            beyond=np.minimum(abs(local[:,0]-1),abs(local[:,0]+1))>half+.02
            crossed.update(current[changed&beyond].tolist())
        if setup['layout']=='D_cross_baffles':
            changed=np.signbit(origin[:,1])!=np.signbit(local[:,1])
            cross_y.update(current[changed&(abs(local[:,1])>half+.02)].tolist())
    report=dict(case=name,samples=len(rr),particle_count=len(ids),max_abs_local_x_m=max_x,max_local_radius_m=max_radius,
                outer_escape=max_x>3.16 or max_radius>1.16,
                unique_particles_crossing_sealed_transverse_plates=len(crossed),
                unique_particles_crossing_longitudinal_plate=len(cross_y),
                density_range_kg_m3=[rho_min,rho_max],max_frame_speed_m_s=speed,
                passed=not crossed and not cross_y and max_x<=3.16 and max_radius<=1.16,
                note='Particle-ID isolation checked at every saved state in the actual instantaneous body frame. This is boundary integrity, not spatial convergence.')
    (folder/'fluid_boundary_qa.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case');p.add_argument('--smoke',action='store_true');p.add_argument('--fluid-qa',action='store_true')
    a=p.parse_args();analyze(a.case,a.smoke)
    if a.fluid_qa:fluid_qa(a.case)
