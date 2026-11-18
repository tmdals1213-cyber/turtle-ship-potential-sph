# SPDX-License-Identifier: GPL-3.0-or-later
"""Read actual SPH VTK/CSV outputs; separate startup health from wave validation."""
from pathlib import Path
import csv,json,sys
import numpy as np
import pyvista as pv
from scipy.spatial.transform import Rotation
ROOT=Path(__file__).resolve().parent
def summarize(name):
    folder=ROOT/'sph'/name;out=folder/'output'
    with (out/'Run.csv').open() as f:
        rows=list(csv.DictReader(f,delimiter=';'));run=rows[-1]
    num=lambda key:float(run[key].replace(',',''))
    paths=sorted((out/'vtks').glob('PartVtk_*.vtk'))
    end=pv.read(paths[-1]);start=pv.read(paths[0]);fluid=end['Type']==3
    part=int(paths[-1].stem.split('_')[-1]);snapshot_time=num('PhysicalTime')
    with (out/'RunPARTs.csv').open() as pf:
        partrows=list(csv.DictReader(pf,delimiter=';'))
    match=[r for r in partrows if r.get('Part','').isdigit() and int(r['Part'])==part]
    if match:snapshot_time=float(match[-1]['TimeStep [s]'].replace(',',''))
    speed=np.linalg.norm(end['Vel'][fluid],axis=1);rho=end['Rho'][fluid]
    res=dict(case=name,status='startup_completed_not_wave_validated',physical_time_s=num('PhysicalTime'),solver_time_s=num('TSimul'),
        total_time_s=num('TTotal'),steps=int(num('Steps')),initial_particles=int(num('Np')),max_particles=int(num('MaxParticles')),
        excluded_particles=int(num('PartsOut')),density_exclusions=int(num('PartsOutRho')),velocity_exclusions=int(num('PartsOutVel')),
        gpu_reported_particle_and_aux_memory_MiB=num('MemGpu')/1024**2,hardware=run['Hardware'],
        fluid_speed_max_m_s=float(speed.max()),fluid_density_range_kg_m3=[float(rho.min()),float(rho.max())],
        wall_seconds_per_physical_second=num('TSimul')/num('PhysicalTime'),
        mean_dt_s=num('PhysicalTime')/num('Steps'),final_vtk=str(paths[-1].relative_to(ROOT)))
    res['last_saved_snapshot_time_s']=snapshot_time
    invalid_hollow={'smoke_hull','pilot_hull','rest_dbc_30','rest_dbc_15'}
    if name in invalid_hollow:
        res['status']='invalid_hydrostatics_internal_fluid_in_original_setup'
        res['invalid_reason']='Original STL shell enclosed fluid particles. Do not use buoyancy or motion as physical validation.'
    elif name.startswith('rest_solid') or name.startswith('rest_defaultmass') or name.startswith('rest_mdbc_advanced'):
        res['status']='corrected_solid_hull_hydrostatic_diagnostic'
    elif name.startswith('piston_2d'):
        res['status']='completed_wave_maker_diagnostic_not_ship_response'
    # Completion and numerical health are separate. These conservative screens
    # flag diagnostics for review; they are not an accuracy validation criterion.
    health=[]
    if rho.min()<900 or rho.max()>1100:health.append('Fluid density deviates more than 10% from reference.')
    if name.startswith('piston_2d') and speed.max()>3:health.append('Wave-tank fluid speed exceeds 3 m/s; inspect localized boundary instability.')
    if health:
        res['status']='completed_but_numerically_unhealthy_do_not_use_as_validation'
        res['health_flags']=health
    if np.any(start['Type']==2):
        setup=json.loads((folder/'setup.json').read_text(encoding='utf-8'))
        cg=np.array(setup.get('center_of_mass_m',json.loads((ROOT/'case.json').read_text(encoding='utf-8'))['hull']['center_of_mass_m']))
        def floats(mesh):
            k=mesh['Type']==2;ids=mesh['Idp'][k];return mesh.points[k][np.argsort(ids)]
        a,b=floats(start),floats(end);aa,bb=a.mean(0),b.mean(0)
        u,_,vt=np.linalg.svd((a-aa).T@(b-bb));rot=vt.T@u.T
        if np.linalg.det(rot)<0:vt[-1]*=-1;rot=vt.T@u.T
        cg1=rot@(cg-aa)+bb
        res['hull_final_cg_m']=cg1.tolist();res['hull_heave_m']=float(cg1[2]-cg[2]);res['hull_pitch_deg']=float(np.degrees(np.arctan2(rot[0,2],rot[0,0])))
        res['rigid_fit_rms_m']=float(np.sqrt(np.mean(((a-aa)@rot.T+bb-b)**2)))
        if name.startswith('tow_'):
            time_s=snapshot_time;expected=.125*time_s**2 if time_s<=2 else .5+.5*(time_s-2)
            res['hull_surge_displacement_m']=float(cg1[0]-cg[0]);res['prescribed_surge_integral_m']=expected
            res['surge_integral_error_m']=float(cg1[0]-cg[0]-expected)
            res['status']='completed_3d_startup_not_seakeeping_validation' if not health else res['status']
        elif name.startswith('free6_') or 'propulsion' in setup:
            res['hull_surge_displacement_m']=float(cg1[0]-cg[0]);res['hull_euler_xyz_deg']=Rotation.from_matrix(rot).as_euler('xyz',degrees=True).tolist()
            res['prescribed_body_velocity']=False;res['free_dofs']=['Surge','Sway','Heave','Roll','Pitch','Yaw']
            res['status']='completed_free_6dof_demonstration_not_seakeeping_validation' if not health else res['status']
            if 'propulsion' in setup and not health:
                res['status']='completed_force_driven_startup_not_resistance_validation' if 'pilot' in name else 'completed_free_6dof_propulsion_diagnostic_not_resistance_validation'
    (folder/'summary.json').write_text(json.dumps(res,indent=2));return res
if __name__=='__main__':
    results=[summarize(n) for n in (sys.argv[1:] or ['smoke_hull','pilot_hull','pilot_flume'])]
    aggregate=ROOT/'results/sph_startup_summary.json'
    previous=json.loads(aggregate.read_text(encoding='utf-8')) if aggregate.exists() else []
    merged={r['case']:r for r in previous};merged.update({r['case']:r for r in results})
    aggregate.write_text(json.dumps(list(merged.values()),indent=2));print(json.dumps(results,indent=2))
