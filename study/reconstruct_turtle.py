# SPDX-License-Identifier: GPL-3.0-or-later
from project_paths import ASSET_DIR, SOLVER_BIN, FONT_DIR
"""Fit a symmetric wetted outer envelope to the local licensed visual hull.
All derived meshes remain in work/, not in the redistributable project.
"""
from pathlib import Path
import json,sys,argparse
import numpy as np
import trimesh
from scipy.ndimage import binary_closing,distance_transform_edt,gaussian_filter
from scipy.interpolate import RegularGridInterpolator
from skimage.measure import marching_cubes
ROOT=Path(__file__).resolve().parent;out=ASSET_DIR
ap=argparse.ArgumentParser();ap.add_argument('--extended',action='store_true');args=ap.parse_args()
a=np.load(out/'turtle_source_mesh.npz');v=(a['vertices']-[-1.70,.0198,-4.90])*3/14.6
names=a['names'];selected=['Cube.000','Cube.003','Cube.004','Cube.006','Cube.008','Cube.009']+(['Cube.001','Cube.005'] if args.extended else []);ids=[i for i,n in enumerate(names) if n in selected]
tri=v[a['faces'][np.isin(a['object_ids'],ids)]]
dx=.008;dz=.004;x=np.arange(-2.40,1.55+dx/2,dx);z=np.arange(-.26,(.75 if args.extended else .21)+dz/2,dz);b=np.full((len(x),len(z)),-np.inf)
for t in tri:
    if t[:,1].max()<0:continue
    px,pz,py=t[:,0],t[:,2],np.abs(t[:,1]);den=(pz[1]-pz[2])*(px[0]-px[2])+(px[2]-px[1])*(pz[0]-pz[2])
    if abs(den)<1e-14:continue
    i0=max(0,int(np.ceil((px.min()-x[0])/dx)));i1=min(len(x)-1,int(np.floor((px.max()-x[0])/dx)))
    j0=max(0,int(np.ceil((pz.min()-z[0])/dz)));j1=min(len(z)-1,int(np.floor((pz.max()-z[0])/dz)))
    if i1<i0 or j1<j0:continue
    xx=x[i0:i1+1,None];zz=z[None,j0:j1+1]
    aa=((pz[1]-pz[2])*(xx-px[2])+(px[2]-px[1])*(zz-pz[2]))/den
    bb=((pz[2]-pz[0])*(xx-px[2])+(px[0]-px[2])*(zz-pz[2]))/den;cc=1-aa-bb
    val=aa*py[0]+bb*py[1]+cc*py[2];valid=(aa>=-1e-6)&(bb>=-1e-6)&(cc>=-1e-6)
    part=b[i0:i1+1,j0:j1+1];part[:]=np.maximum(part,np.where(valid,val,-np.inf))
valid=np.isfinite(b)&(b>1e-6);filled=binary_closing(valid,structure=np.ones((5,5)))|valid
holes=filled&~valid;_,near=distance_transform_edt(~valid,return_indices=True);b[holes]=b[tuple(near[:,holes])];b[~filled]=0
# Millimetre smoothing removes artistic plank seams, not the flat-bottom form.
b=gaussian_filter(b,.65);b[b<.002]=0
b[:, -2:]=0  # Close the computational deck above the waterline.
np.savez_compressed(out/('turtle_extended_half_breadth.npz' if args.extended else 'turtle_half_breadth.npz'),x=x,z=z,half_breadth=b)
dy=.008;y=np.arange(-.6,.61,dy);field=np.where(b>0,b,-.008)[:,None,:]-np.abs(y)[None,:,None]
vv,ff,_,_=marching_cubes(field.astype('f4'),0,spacing=(dx,dy,dz));vv+=[x[0],y[0],z[0]]
m=trimesh.Trimesh(vv,ff,process=True);m.fix_normals(multibody=True);m=max(m.split(only_watertight=False),key=lambda q:abs(q.volume));m.fix_normals()
print('mesh audit',m.is_watertight,m.volume,m.bounds,flush=True)
assert m.is_watertight and m.volume>0
if args.extended:
    m.export(out/'turtle_extended_closed.stl');print('Extended topside envelope saved; baseline model preserved.');raise SystemExit()
m.export(out/'turtle_fitted_closed.stl');wet=m.slice_plane([0,0,0],[0,0,-1],cap=False);wet.export(out/'turtle_fitted_wetted.stl')
vol=wet.volume;cb=wet.center_mass;cg=np.array([cb[0],0,.08]);rho=1000;g=9.81
interp=RegularGridInterpolator((x,z),b,bounds_error=False,fill_value=0);wp=2*interp(np.c_[x,np.zeros_like(x)])
A=np.trapezoid(wp,x);S=np.trapezoid(wp*(x-cg[0]),x);Iw=np.trapezoid(wp*(x-cg[0])**2,x)
mass=vol*rho;C=rho*g*np.array([[A,-S],[-S,Iw+vol*(cb[2]-cg[2])]])
assert np.linalg.eigvalsh(C).min()>0
hyd=dict(volume_m3=float(vol),mass_kg=float(mass),center_of_buoyancy_m=cb.tolist(),center_of_mass_m=cg.tolist(),inertia_kg_m2=(mass*np.array([.35, .75,.75])**2).tolist(),waterplane_area_m2=float(A),hydrostatic_stiffness=C.tolist(),wetted_bounds_m=wet.bounds.tolist(),closed_watertight=m.is_watertight,assumptions=['Symmetric sampled outer envelope of six main hull objects; rudder, oars and planking seams excluded.','Selected waterline, 3 m nominal model length, CGz=0.08 m and radii of gyration are illustrative; not a historical reconstruction.'])
(ROOT/'results/turtle_fitted_hydrostatics.json').write_text(json.dumps(hyd,indent=2))
for count in [1800,2800]:
    sm=wet.simplify_quadric_decimation(face_count=count);sm.fix_normals();sm.export(out/f'turtle_bem_{count}.stl')
c=json.loads((ROOT/'case_turtle_proxy.json').read_text(encoding='utf-8'));c['waves']['baseline_height_m']=.24;c['hull']['center_of_mass_m']=cg.tolist();c['description']='Fitted turtle-ship main-hull envelope; illustrative mass distribution and selected waterline.'
(ROOT/'case_turtle_fitted.json').write_text(json.dumps(c,indent=2));print(json.dumps(hyd,indent=2))
