"""Check actual generated occupancy, lattice phase, body inputs and wall joins."""
import json, sys, xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
import pyvista as pv

ROOT=Path(__file__).resolve().parent
for name in sys.argv[1:]:
    folder=ROOT/'sph'/name; setup=json.loads((folder/'setup.json').read_text())
    d=pv.read(folder/'generated/Case_All.vtk'); dp=setup['dp_m']
    fluid=d.points[d['Type']==3]; body=d.points[d['Type']==2]
    result=dict(case=name,total_particles=d.n_points,fluid_particles=len(fluid),hull_particles=len(body),
        moving_particles=int(sum(d['Type']==1)),fixed_particles=int(sum(d['Type']==0)),
        top_water_particle_z_m=float(fluid[:,2].max()))
    passed=np.isfinite(d.points).all() and abs(result['top_water_particle_z_m']+dp/2)<dp*.01
    # The piston must not own stationary floor particles at the shared corner.
    moving=d.points[d['Type']==1]; fixed=d.points[d['Type']==0]
    depth=setup['water_depth_m']; x1=setup['domain_bounds_m'][1]
    result['moving_particles_below_floor_top']=int(sum(moving[:,2]<-depth))
    bottom=fixed[(fixed[:,0]>x1-.2)&(fixed[:,0]<x1+.5)&(fixed[:,2]<-depth)]
    result['fixed_floor_particles_under_piston']=len(bottom)
    passed=passed and result['moving_particles_below_floor_top']==0 and len(bottom)>0
    if setup['physical_model'] is not None:
        p=setup['physical_model']; bounds=np.array(p['bounds_m'])
        near=fluid[np.all((fluid>=bounds[0]-dp)&(fluid<=bounds[1]+dp),axis=1)]
        surface=pv.read(p['hull_file']).clean().compute_normals(auto_orient_normals=True,consistent_normals=True)
        dist=pv.PolyData(near).compute_implicit_distance(surface)['implicit_distance']
        deep=int(sum(dist < -dp)); result['deep_interior_fluid_particles']=deep
        xml=ET.parse(folder/'generated/Case.xml').getroot()
        native=xml.find('.//particles/floating'); decl=xml.find('.//floatings/floating')
        mass=float(native.find('massbody').get('value'))
        rel=abs(mass-p['mass_kg'])/p['mass_kg']
        dofs=[int(decl.find(tag).get(axis)) for tag in ['translationDOF','rotationDOF'] for axis in 'xyz']
        result.update(native_mass_kg=mass,mass_serialization_relative_error=rel,all_six_dofs_free=all(dofs))
        passed=passed and deep==0 and len(body)>0 and rel<1e-6 and all(dofs)
        passed=passed and np.allclose([float(decl.find('inertia').get(ax)) for ax in 'xyz'],p['inertia_kg_m2'])
    result['passed']=bool(passed)
    (folder/'initial_audit.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result),flush=True)
    if not passed: raise RuntimeError('Generated input audit failed; inspect initial_audit.json')
