# SPDX-License-Identifier: GPL-3.0-or-later
from project_paths import ASSET_DIR, SOLVER_BIN, FONT_DIR
"""A declared, history-informed loading scenario, not an exact 1592 reconstruction.

Licensed geometry and every derivative remain in the local work directory.
Mass is specified first; sinkage and trim then follow from buoyancy equilibrium.
"""
from pathlib import Path
import json
import numpy as np
import trimesh
from scipy.optimize import root, brentq

ROOT = Path(__file__).resolve().parent
GEO = ASSET_DIR
RHO, GRAV, SCALE = 1025., 9.81, 10.
source = trimesh.load(GEO / 'turtle_extended_closed.stl', process=True)
# 35.27 - 9.21 is a derived non-tail reference, NOT a reported Lpp/Lwl.
reference_length = 35.27 - 9.21
s = reference_length / SCALE / np.ptp(source.bounds[:, 0])
keel = source.bounds[0, 2]
mesh = source.copy()
mesh.vertices = (source.vertices - [0, 0, keel]) * s

loads = [
    dict(item='structure', tonnes=140.4, kg_above_keel_m=2., lcg_from_midship_m=1.,
         basis='Mass: 2026 official turtle-ship reconstruction announcement. CG: Panokseon analogue assumption, not measured turtle CG.'),
    dict(item='crew_and_belongings', tonnes=9.8, kg_above_keel_m=46.69/9.8, lcg_from_midship_m=-1.435/9.8,
         basis='140 people at 70 kg; loading example transcribed from 2021 Panokseon report. Analogue, not turtle complement fact.'),
    dict(item='weapons_and_ammunition', tonnes=1.2, kg_above_keel_m=3., lcg_from_midship_m=0.,
         basis='Mass follows Panokseon example; placement assumed. Exact equipment included in published self-weight is unresolved.'),
    dict(item='stores', tonnes=3., kg_above_keel_m=1., lcg_from_midship_m=0., basis='Demonstration assumption.'),
]
tonnes = sum(q['tonnes'] for q in loads)
kg_full = sum(q['tonnes']*q['kg_above_keel_m'] for q in loads)/tonnes
lcg_full = sum(q['tonnes']*q['lcg_from_midship_m'] for q in loads)/tonnes
mass = tonnes*1000/SCALE**3

def rotate(pitch):
    c, s = np.cos(pitch), np.sin(pitch)
    return np.array([[c,0,s],[0,1,0],[-s,0,c]])

def wet_at(z):
    m = mesh.copy(); m.apply_translation([0,0,z])
    return m.slice_plane([0,0,0],[0,0,-1],cap=False)

z0 = brentq(lambda z: wet_at(z).volume-mass/RHO, -.55, -.02)
untrimmed = wet_at(z0)
midship = untrimmed.bounds[:,0].mean()
cg_keel = np.array([midship+lcg_full/SCALE, 0, kg_full/SCALE])

def posed(z, pitch):
    m = mesh.copy(); m.vertices = (mesh.vertices-cg_keel) @ rotate(pitch).T + cg_keel + [0,0,z]
    return m

def residual(q):
    m=posed(*q); wet=m.slice_plane([0,0,0],[0,0,-1],cap=False)
    return [wet.volume/(mass/RHO)-1, (wet.center_mass[0]-cg_keel[0])/reference_length*SCALE]

sol = root(residual,[z0,0],tol=1e-9)
assert sol.success and np.linalg.norm(residual(sol.x))<1e-7, sol.message
shift, pitch=sol.x; closed=posed(shift,pitch); wet=closed.slice_plane([0,0,0],[0,0,-1],cap=False)
cg=cg_keel+[0,0,shift]; cb=wet.center_mass; volume=float(wet.volume)

section=closed.section(plane_origin=[0,0,0],plane_normal=[0,0,1])
polys=section.discrete
assert len(polys)==1, 'Expected a single simple waterplane polygon.'
poly=polys[0][:,:2]-cg[:2]
if np.allclose(poly[0],poly[-1]):poly=poly[:-1]
x,y=poly.T; xx,yy=np.roll(poly,-1,axis=0).T; cross=x*yy-xx*y
if cross.sum()<0:
    poly=poly[::-1];x,y=poly.T;xx,yy=np.roll(poly,-1,axis=0).T;cross=x*yy-xx*y
area=cross.sum()/2; sx=((x+xx)*cross).sum()/6; sy=((y+yy)*cross).sum()/6
ix=((y*y+y*yy+yy*yy)*cross).sum()/12
iy=((x*x+x*xx+xx*xx)*cross).sum()/12
ixy=((2*x*y+x*yy+xx*y+2*xx*yy)*cross).sum()/24
C=np.zeros((6,6)); sub=np.array([[area,sy,-sx],[sy,ix+volume*(cb[2]-cg[2]),-ixy],[-sx,-ixy,iy+volume*(cb[2]-cg[2])]])*RHO*GRAV
C[np.ix_([2,3,4],[2,3,4])]=sub
assert np.linalg.eigvalsh(sub).min()>0
# Inertia is an explicitly assumed diagonal approximation in the initial axes.
beam=np.ptp(closed.bounds[:,1]); radii=np.array([.35*beam,.25*reference_length/SCALE,.25*reference_length/SCALE])
inertia=mass*radii**2
gm_roll=(ix-sy*sy/area)/volume+cb[2]-cg[2]
gm_pitch=(iy-sx*sx/area)/volume+cb[2]-cg[2]
closed.export(GEO/'turtle_history_closed.stl');wet.export(GEO/'turtle_history_wetted.stl')
# STL float32 rounding can collapse tiny marching-cubes triangles. Remove only
# those zero-area/duplicate faces and verify the exported, reloaded shell.
disk=trimesh.load(GEO/'turtle_history_closed.stl',process=True)
disk.update_faces(disk.nondegenerate_faces());disk.update_faces(disk.unique_faces());disk.remove_unreferenced_vertices()
assert disk.is_watertight and abs(disk.volume/closed.volume-1)<1e-6
disk.export(GEO/'turtle_history_closed.stl')
# Decimate a closed surface first; preserving thousands of waterline boundary
# edges while decimating an open mesh unnecessarily inflates the BEM system.
coarse=closed.simplify_quadric_decimation(face_count=2000);coarse.fix_normals()
small=coarse.slice_plane([0,0,0],[0,0,-1],cap=False);small.export(GEO/'turtle_history_bem.stl')
assert abs(small.volume/volume-1)<.01
visual=np.load(GEO/'turtle_visual.npz')
vv=(visual['vertices']-[0,0,keel])*s
vv=(vv-cg_keel)@rotate(pitch).T+cg_keel+[0,0,shift]
np.savez_compressed(GEO/'turtle_history_visual.npz',vertices=vv,faces=visual['faces'],colors=visual['colors'],waterplane=poly+cg[:2])

assumptions=[
    '2026 reconstruction studies 1795 documents; it does not establish the exact 1592 vessel.',
    'Published 35.27 m overall length includes a 9.21 m tail. The difference 26.06 m is only a derived scaling reference, not Lpp or Lwl.',
    'Artistic visual mesh and fitted symmetric envelope are not the official reconstruction CAD.',
    'Mass budget, CG and radii of gyration are declared plausible loading assumptions; exact historical mass distribution is unknown.',
    'Model is watertight; openings, flooding, oars, sails, rudder forces and wind are excluded.',
    '1:10 Froude scaling represents gravity-wave motions, not simultaneous full-scale Reynolds similarity.',
]
meta=dict(status='equilibrium_computed',scale_ratio=SCALE,rho_kg_m3=RHO,mass_kg=mass,full_scale_mass_t=tonnes,
    loads=loads,reference_hull_length_m=reference_length/SCALE,reference_full_scale_length_m=reference_length,
    published_LOA_m=35.27,published_tail_m=9.21,published_self_weight_t=140.4,
    estimated_full_scale_KG_m=kg_full,estimated_full_scale_LCG_m=lcg_full,
    original_to_model_scale=s,original_keel_z=keel,cg_before_equilibrium_m=cg_keel.tolist(),equilibrium_shift_z_m=float(shift),
    initial_trim_deg=float(np.degrees(pitch)),center_of_mass_m=cg.tolist(),center_of_buoyancy_m=cb.tolist(),
    volume_m3=volume,inertia_kg_m2=inertia.tolist(),inertia_policy='Assumed radii: roll 0.35 B, pitch/yaw 0.25 reference L; diagonal in initial axes.',
    hydrostatic_stiffness_6dof=C.tolist(),waterplane_area_m2=float(area),GM_roll_m=float(gm_roll),GM_pitch_m=float(gm_pitch),
    draft_m=float(-wet.bounds[0,2]),beam_m=float(beam),closed_volume_m3=float(closed.volume),
    equilibrium_residual=residual(sol.x),closed_watertight=bool(closed.is_watertight),bem_panel_count=len(small.faces),bem_volume_relative_error=float(small.volume/volume-1),assumptions=assumptions,
    KG_sensitivity_full_scale_m=[2.,2.5,3.],GM_roll_fixed_geometry_full_scale_m=[float(gm_roll*SCALE+kg_full-z) for z in [2.,2.5,3.]],
    sources=[
        dict(role='Official reconstruction announcement reproduced by KCDF',url='https://www.kculture.or.kr/brd/board/399/L/menu/498?bbIdx=15970&brdType=R',access='Public announcement read; report download server returned an HTML error.'),
        dict(role='Official report catalogue',url='https://www.seamuse.go.kr/resources/academicreport/info/4452'),
        dict(role='2021 Panokseon analogue; public text transcription, original PDF table not visually inspected',url='https://www.scribd.com/document/676848883/판옥선-학술-복원-보고서'),
    ])
(ROOT/'results/turtle_history_model.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
c=dict(study='History-informed turtle-ship wave demonstration',hull=dict(type='fitted_turtle_envelope',length_m=reference_length/SCALE,center_of_mass_m=cg.tolist()),
    environment=dict(rho_kg_m3=RHO,g_m_s2=GRAV,water_depth_m=1.),
    motion=dict(initial_velocity_m_s=[.5,0,0],forward_speed_m_s=.5,wave_direction_rad=float(np.pi),free_dofs=['Surge','Sway','Heave','Roll','Pitch','Yaw'],prescribed_or_constrained_dofs=[],propulsion='None; SPH has initial speed only. Potential is a linear response at reference mean speed, not free-coasting transient.'),
    waves=dict(baseline_wavelength_over_length=3/(reference_length/SCALE),baseline_height_m=.24),
    full_scale_equivalent=dict(initial_speed_m_s=.5*np.sqrt(SCALE),wave_height_m=2.4,wavelength_m=30,water_depth_m=10),
    accuracy_claim='History-informed educational scenario, not a historically validated ship or performance prediction.')
(ROOT/'case_turtle_history.json').write_text(json.dumps(c,indent=2))
print(json.dumps({k:meta[k] for k in ['mass_kg','full_scale_mass_t','draft_m','beam_m','initial_trim_deg','estimated_full_scale_KG_m','GM_roll_m','bem_panel_count','equilibrium_residual']},indent=2))
