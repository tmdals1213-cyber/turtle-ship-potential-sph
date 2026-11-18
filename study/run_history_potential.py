# SPDX-License-Identifier: GPL-3.0-or-later
from project_paths import ASSET_DIR, SOLVER_BIN, FONT_DIR
"""Independent six-DOF linear BEM at a reference speed, NOT coasting SPH."""
from pathlib import Path
import os,sys,json,time,copy
os.environ.setdefault('OMP_NUM_THREADS','4');os.environ.setdefault('OPENBLAS_NUM_THREADS','4')
import numpy as np
import capytaine as cpt
import trimesh
from matplotlib.path import Path as PolygonPath
from capytaine.bem.airy_waves import froude_krylov_force,airy_waves_potential,airy_waves_velocity
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from common import ROOT,wave
c=json.loads((ROOT/'case_turtle_history.json').read_text(encoding='utf-8'));p=json.loads((ROOT/'results/turtle_history_model.json').read_text(encoding='utf-8'));w=wave(c=c)
geo=ASSET_DIR;out=ROOT/'results/potential_turtle_history';out.mkdir(exist_ok=True)
tri=trimesh.load(geo/'turtle_history_bem.stl',process=True)
mesh=cpt.Mesh(vertices=tri.vertices,faces=np.c_[tri.faces,tri.faces[:,-1]],name='history_informed_turtle')
cg=p['center_of_mass_m'];names=['Surge','Sway','Heave','Roll','Pitch','Yaw']
body=cpt.FloatingBody(mesh=mesh,dofs=cpt.rigid_body_dofs(rotation_center=cg),center_of_mass=cg,mass=p['mass_kg']).with_only_dofs(names)
solver=cpt.BEMSolver();U=c['motion']['forward_speed_m_s'];g=9.81;rho=c['environment']['rho_kg_m3']
M=np.diag([p['mass_kg']]*3+p['inertia_kg_m2']);C=np.array(p['hydrostatic_stiffness_6dof'])
t0=time.perf_counter();kw=dict(body=body,omega=w['omega_rad_s'],wave_direction=np.pi,forward_speed=U,water_depth=1,rho=rho,g=g)
rr=[]
for n in names:
    rr.append(solver.solve(cpt.RadiationProblem(radiating_dof=n,**kw)));print('Radiation',n,round(time.perf_counter()-t0,1),flush=True)
dd=solver.solve(cpt.DiffractionProblem(**kw))
A=np.array([[r.added_mass[n] for r in rr] for n in names]);B=np.array([[r.radiation_damping[n] for r in rr] for n in names])
fk=froude_krylov_force(dd.problem);F=np.array([dd.forces[n]+fk[n] for n in names]);we=dd.encounter_omega
Z=C-we**2*(M+A)-1j*we*B;xi=np.linalg.solve(Z,F)
res=float(np.linalg.norm(Z@xi-F)/np.linalg.norm(F));assert res<1e-10
combined=copy.copy(dd);combined.sources=dd.sources+sum(xi[j]*rr[j].sources for j in range(6))
np.savez_compressed(geo/'turtle_history_solution_panels.npz',vertices=mesh.vertices,faces=mesh.faces,centers=mesh.faces_centers,normals=mesh.faces_normals,sources=combined.sources)
np.savez_compressed(out/'baseline.npz',A=A,B=B,F=F,xi=xi,omega=w['omega_rad_s'],encounter_omega=we,M=M,C=C)
amps=abs(xi)*w['amplitude_m']
row=dict(amplitude_m=dict(zip(names[:3],amps[:3].tolist())),amplitude_deg=dict(zip(names[3:],np.degrees(amps[3:]).tolist())),equation_residual=res,solve_wall_s=time.perf_counter()-t0)
meta=dict(status='response_complete_field_running',panel_count=mesh.nb_faces,dofs=names,reference_mean_speed_m_s=U,rows=[row],
    assumptions=p['assumptions']+['Linear harmonic response about a fixed mean speed. This does not predict free-coasting deceleration, viscous resistance or startup; SPH is a separate initial-value problem.','Second-order incident shape plus linear scattered waves is a visualization approximation, not a full nonlinear solution.','The approximately 12.6 degree pitch amplitude at this large wave height exceeds a comfortable small-motion regime; use qualitatively, not as a quantitatively reliable large-wave prediction.'])
(out/'baseline.json').write_text(json.dumps(meta,indent=2));print(json.dumps(row),flush=True)
poly=np.load(geo/'turtle_history_visual.npz')['waterplane'];boundary=PolygonPath(poly)
x=np.linspace(-2.8,2.8,131);y=np.linspace(-1.25,1.25,71);xx,yy=np.meshgrid(x,y,indexing='ij');pts=np.c_[xx.ravel(),yy.ravel(),np.zeros(xx.size)]
valid=~boundary.contains_points(pts[:,:2],radius=.05);eta=np.full(len(pts),np.nan+0j);indices=np.flatnonzero(valid)
for st in range(0,len(indices),384):
    ii=indices[st:st+384];eta[ii]=solver.compute_pressure(pts[ii],combined)/(rho*g)
    if st%1920==0:print('surface',st,'/',len(indices),flush=True)
inc=np.exp(-1j*w['k_rad_m']*pts[:,0]);phi=airy_waves_potential(pts[:50],dd.problem);vel=airy_waves_velocity(pts[:50],dd.problem)
identity=(1j*we*phi+U*vel[:,0])/g;error=float(np.max(abs(identity-inc[:50])));assert error<1e-10
np.savez_compressed(out/'surface.npz',x=x,y=y,eta_scattered=eta.reshape(xx.shape),eta_incident=inc.reshape(xx.shape),valid=valid.reshape(xx.shape))
meta.update(status='complete',total_wall_s=time.perf_counter()-t0,incident_bernoulli_identity_error=error,max_unit_scattered_surface=float(np.nanmax(abs(eta))),mesh_volume_m3=float(body.volume))
(out/'baseline.json').write_text(json.dumps(meta,indent=2));print(json.dumps(meta,indent=2))
