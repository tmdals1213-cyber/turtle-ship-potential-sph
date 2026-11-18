# SPDX-License-Identifier: GPL-3.0-or-later
from project_paths import ASSET_DIR, SOLVER_BIN, FONT_DIR
"""Audit the actual generated six-DOF case against its transformed hull."""
from pathlib import Path
import json,sys,xml.etree.ElementTree as ET
import numpy as np
import pyvista as pv
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parent
name=sys.argv[1] if len(sys.argv)>1 else 'free6_history_40'
folder=ROOT/'sph'/name;setup=json.loads((folder/'setup.json').read_text(encoding='utf-8'));p=json.loads((ROOT/'results/turtle_history_model.json').read_text(encoding='utf-8'))
geo=ASSET_DIR;fit=np.load(geo/'turtle_extended_half_breadth.npz')
breadth=RegularGridInterpolator((fit['x'],fit['z']),fit['half_breadth'],bounds_error=False,fill_value=0)
m=pv.read(folder/'generated/Case_All.vtk');offset=np.array(setup['center_of_mass_m'])-p['center_of_mass_m'];dp=setup['dp_m']
points=m.points[m['Type']==3]-offset
cg0=np.array(p['cg_before_equilibrium_m']);q=np.radians(p['initial_trim_deg']);c,s=np.cos(q),np.sin(q);R=np.array([[c,0,s],[0,1,0],[-s,0,c]])
original=((points-cg0-[0,0,p['equilibrium_shift_z_m']])@R+cg0)/p['original_to_model_scale']+[0,0,p['original_keel_z']]
b=breadth(original[:,[0,2]]);inside=(b>0)&(abs(original[:,1])<b)
# Distances to densely sampled surface vertices provide a conservative exclusion
# band for tiny differences between the half-breadth raster and marching mesh.
hm=pv.read(geo/'turtle_history_closed.stl');dist=cKDTree(hm.points).query(points[inside],workers=2)[0]
interior=int(sum(dist>dp));nf=int(sum(m['Type']==2));vol=nf*dp**3
root=ET.parse(folder/'generated/Case.xml').getroot();f=root.find('.//particles/floating')
assert f is not None
# GenCase omits the optional element when all three components are zero.
ini=f.find('linearvelini');velocity={ax:float(ini.get(ax)) if ini is not None else 0.0 for ax in 'xyz'}
assert np.allclose([velocity[ax] for ax in 'xyz'],setup['initial_velocity_m_s'])
definition=root.find('.//casedef/floatings/floating')
if definition is None:definition=root.find('.//floatings/floating')
q=dict(case=name,total_particles=m.n_points,interior_fluid_particles_excluding_one_dp_surface_band=interior,floating_particles=nf,
    particle_volume_m3=vol,closed_volume_m3=setup['closed_hull_volume_m3'],particle_volume_error_percent=100*(vol/setup['closed_hull_volume_m3']-1),
    massbody_kg=float(f.find('massbody').get('value')),masspart_kg=float(f.find('masspart').get('value')),initial_velocity_m_s=velocity,
    six_dofs_enabled=all(definition.find(tag).get(ax)=='1' for tag in ['translationDOF','rotationDOF'] for ax in 'xyz'),
    no_prescribed_velocity=definition.find('linearvel') is None,
    meaning='Occupancy, rigid-body input and discrete-volume audit only, not hydrodynamic validation.')
(folder/'initial_audit.json').write_text(json.dumps(q,indent=2));print(json.dumps(q,indent=2))
assert interior==0 and q['six_dofs_enabled'] and q['no_prescribed_velocity']
