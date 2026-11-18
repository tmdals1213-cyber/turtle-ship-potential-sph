# SPDX-License-Identifier: GPL-3.0-or-later
from project_paths import ASSET_DIR, SOLVER_BIN, FONT_DIR
"""Measured rigid poses and a visualization surface reconstructed from SPH particles.

No potential-flow samples are used in this SPH processing path. Surface smoothing
is for display only; velocities, densities and poses use the actual solver output.
"""
from pathlib import Path
import argparse,csv,json,time
import numpy as np
import pyvista as pv
from scipy.ndimage import gaussian_filter
from scipy.spatial.transform import Rotation
ROOT=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('--case',default='free6_history_40');ap.add_argument('--latest-only',action='store_true');a=ap.parse_args()
folder=ROOT/'sph'/a.case;setup=json.loads((folder/'setup.json').read_text(encoding='utf-8'));model=json.loads((ROOT/'results/turtle_history_model.json').read_text(encoding='utf-8'))
cache=ASSET_DIR/f'{a.case}_render';cache.mkdir(exist_ok=True)
rows=[r for r in csv.DictReader((folder/'output/RunPARTs.csv').open(),delimiter=';') if r.get('Part','').isdigit()]
times={int(r['Part']):float(r['TimeStep [s]']) for r in rows}
start=pv.read(folder/'output/vtks/PartVtk_0000.vtk');sel=start['Type']==2
order=np.argsort(start['Idp'][sel]);initial=start.points[sel][order].astype(float);ids=start['Idp'][sel][order];mean0=initial.mean(0);cg0=np.array(setup['center_of_mass_m'])
dx=setup['dp_m'];axes=[np.arange(-2.8,2.8+dx/2,dx),np.arange(-1.2,1.2+dx/2,dx),np.arange(-.4,.64+dx/2,dx)]
extended=[np.r_[v[0]-dx*np.arange(3,0,-1),v,v[-1]+dx*np.arange(1,4)] for v in axes]
edges=[np.r_[v-dx/2,v[-1]+dx/2] for v in extended]
paths=sorted((folder/'output/vtks').glob('PartVtk_*.vtk'))
if a.latest_only:paths=paths[:1]+paths[-1:]
for path in paths:
    part=int(path.stem.split('_')[-1]);target=cache/f'frame_{part:04d}.npz'
    if target.exists():
        old=np.load(target)
        if 'visualization_version' in old and int(old['visualization_version'])==2:continue
    if part not in times or time.time()-path.stat().st_mtime<2:continue
    m=pv.read(path);sel=m['Type']==2;order=np.argsort(m['Idp'][sel]);assert np.array_equal(m['Idp'][sel][order],ids)
    current=m.points[sel][order].astype(float);mean1=current.mean(0)
    u,sv,vt=np.linalg.svd((initial-mean0).T@(current-mean1));R=vt.T@u.T
    if np.linalg.det(R)<0:vt[-1]*=-1;R=vt.T@u.T
    trans=mean1-R@mean0;cg=R@cg0+trans;fit=(initial@R.T+trans-current)
    fluid=m['Type']==3;pts=m.points[fluid].astype(float);vel=m['Vel'][fluid];rho=m['Rho'][fluid]
    follow=np.array([cg[0]-model['center_of_mass_m'][0],0,0]);local=pts-follow
    hist=np.histogramdd(local,bins=edges)[0];occupancy=gaussian_filter(hist.astype('f4'),sigma=1.1,mode='constant')[3:-3,3:-3,3:-3]
    # A padded occupancy isosurface exposes the water volume at the cropped view edges.
    np.savez_compressed(target,occupancy=occupancy.astype('f2'),x=axes[0],y=axes[1],z=axes[2],R=R,translation=trans,cg=cg,follow=follow,
        time=times[part],visualization_sigma_cells=1.1,visualization_version=2,euler_deg=Rotation.from_matrix(R).as_euler('xyz',degrees=True),fit_rms=np.sqrt(np.mean(fit**2)),
        fluid_speed_max=np.linalg.norm(vel,axis=1).max(),density_min=rho.min(),density_max=rho.max(),part=part)
    if part%10==0:print('Processed',part,times[part],flush=True)
allfiles=sorted(cache.glob('frame_*.npz'));data=[np.load(p) for p in allfiles];tt=np.array([float(q['time']) for q in data]);cg=np.array([q['cg'] for q in data])
velocity=np.gradient(cg,tt,axis=0) if len(tt)>1 else np.array([[.5,0,0]])
if tt[0]==0:velocity[0]=setup['initial_velocity_m_s']
result=dict(case=a.case,status='sampled_actual_sph_motion',times_s=tt.tolist(),cg_m=cg.tolist(),velocity_m_s=velocity.tolist(),
    euler_xyz_deg=[q['euler_deg'].tolist() for q in data],max_rigid_fit_rms_m=max(float(q['fit_rms']) for q in data),
    density_min_kg_m3=min(float(q['density_min']) for q in data),density_max_kg_m3=max(float(q['density_max']) for q in data),
    max_fluid_speed_m_s=max(float(q['fluid_speed_max']) for q in data),
    surface_method='Gaussian occupancy isosurface from actual SPH particles; visual smoothing only. Model-scale dimensions.',
    limitations=setup['limitations'])
(folder/'motion_6dof.json').write_text(json.dumps(result,indent=2))
print(json.dumps({k:v for k,v in result.items() if k not in ['times_s','cg_m','velocity_m_s','euler_xyz_deg']},indent=2));print('Last time',tt[-1],'CG',cg[-1],'velocity',velocity[-1],flush=True)
