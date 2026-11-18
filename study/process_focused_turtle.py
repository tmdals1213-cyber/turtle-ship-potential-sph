# SPDX-License-Identifier: GPL-3.0-or-later
from project_paths import ASSET_DIR, SOLVER_BIN, FONT_DIR
"""Cache actual SPH volumes and rigid poses for the focused-wave film."""
from pathlib import Path
import argparse,csv,json,time
import numpy as np
import pyvista as pv
from scipy.ndimage import gaussian_filter
from scipy.spatial.transform import Rotation
from scipy.interpolate import RegularGridInterpolator

ROOT=Path(__file__).resolve().parent;GEO=ASSET_DIR
ap=argparse.ArgumentParser();ap.add_argument('case');ap.add_argument('--start',type=float,default=0)
ap.add_argument('--stride',type=int,default=1);a=ap.parse_args()
folder=ROOT/'sph'/a.case;cache=GEO/f'{a.case}_focus_render';cache.mkdir(exist_ok=True)
setup=json.loads((folder/'setup.json').read_text(encoding='utf-8'))
roofpath=GEO/'turtle_history_roof_map.npz'
if not roofpath.exists():
    hull=pv.read(GEO/'turtle_history_closed.stl');rx=np.arange(-1.60,1.061,.02);ry=np.arange(-.42,.421,.02)
    rz=np.full((len(rx),len(ry)),np.nan)
    for ix,x in enumerate(rx):
        for iy,y in enumerate(ry):
            hit,_=hull.ray_trace((x,y,1.2),(x,y,-.4),first_point=True)
            if len(hit):rz[ix,iy]=hit[2]
    np.savez_compressed(roofpath,x=rx,y=ry,z=rz)
with np.load(roofpath) as q:roofheight=RegularGridInterpolator((q['x'],q['y']),q['z'],bounds_error=False,fill_value=np.nan)
rows=[r for r in csv.DictReader((folder/'output/RunPARTs.csv').open(),delimiter=';') if r.get('Part','').isdigit()]
times={int(r['Part']):float(r['TimeStep [s]']) for r in rows}
initial=pv.read(folder/'output/vtks/PartVtk_0000.vtk');sel=initial['Type']==2
order=np.argsort(initial['Idp'][sel]);body0=initial.points[sel][order].astype(float);ids=initial['Idp'][sel][order]
mean0=body0.mean(0);cg0=np.array(setup['center_of_mass_m'])
dx=setup['dp_m'];axes=[np.arange(-5,5+dx/2,dx),np.arange(-2.4,2.4+dx/2,dx),np.arange(-.6,1.6+dx/2,dx)]
ext=[np.r_[v[0]-dx*np.arange(3,0,-1),v,v[-1]+dx*np.arange(1,4)] for v in axes]
edges=[np.r_[v-dx/2,v[-1]+dx/2] for v in ext]
for path in sorted((folder/'output/vtks').glob('PartVtk_*.vtk')):
    part=int(path.stem.split('_')[-1]);target=cache/f'frame_{part:04d}.npz'
    if target.exists():
        with np.load(target) as old:
            if 'water_above_upper_hull_particles' in old:continue
    if part not in times or times[part]<a.start or part%a.stride or time.time()-path.stat().st_mtime<2:continue
    m=pv.read(path);sel=m['Type']==2;order=np.argsort(m['Idp'][sel]);assert np.array_equal(m['Idp'][sel][order],ids)
    body=m.points[sel][order].astype(float);mean1=body.mean(0)
    u,sv,vt=np.linalg.svd((body0-mean0).T@(body-mean1));R=vt.T@u.T
    if np.linalg.det(R)<0:vt[-1]*=-1;R=vt.T@u.T
    trans=mean1-R@mean0;cg=R@cg0+trans;fit=body0@R.T+trans-body
    fluid=m['Type']==3;p=m.points[fluid];vel=m['Vel'][fluid];rho=m['Rho'][fluid]
    counts=np.histogramdd(p,bins=edges)[0]
    occ=gaussian_filter(counts.astype('f4'),sigma=.8,mode='constant')[3:-3,3:-3,3:-3]
    # Water above the closed roof's central footprint, in the moving body frame.
    local=(p-trans)@R
    over=(local[:,0]>-.8)&(local[:,0]<.4)&(abs(local[:,1])<.2)&(local[:,2]>.49)&(local[:,2]<.9)
    candidates=(local[:,0]>-1.6)&(local[:,0]<1.06)&(abs(local[:,1])<.42)&(local[:,2]>.2)&(local[:,2]<1.1)
    lp=local[candidates];roof_z=roofheight(lp[:,:2]);clearance=lp[:,2]-roof_z
    upper=(roof_z>.2)&(clearance>0)&(clearance<3*dx)
    above_upper=p[candidates][upper]
    near=(abs(p[:,0])<5)&(abs(p[:,1])<2.4)&(p[:,2]>.05)
    render=(abs(p[:,0])<5.8)&(abs(p[:,1])<3.2)&(p[:,2]>-1.0)&(p[:,2]<2.0)
    temporary=target.with_name('.'+target.name)
    np.savez_compressed(temporary,occupancy=occ.astype('f2'),x=axes[0],y=axes[1],z=axes[2],R=R,translation=trans,cg=cg,
        time=times[part],part=part,fit_rms=np.sqrt(np.mean(fit**2)),euler_deg=Rotation.from_matrix(R).as_euler('xyz',degrees=True),
        density_min=rho.min(),density_max=rho.max(),speed_max=np.linalg.norm(vel,axis=1).max(),
        water_above_roof_particles=int(over.sum()),surface_points=p[near].astype('f4'),surface_vel=vel[near].astype('f4'),
        render_points=p[render].astype('f4'),render_ids=m['Idp'][fluid][render],visualization_version=3,
        water_above_upper_hull_particles=int(upper.sum()),upper_hull_water_points=above_upper.astype('f4'))
    temporary.replace(target)
    if part%10==0:print('Cached',part,round(times[part],3),flush=True)
files=sorted(cache.glob('frame_*.npz'));records=[]
for path in files:
    with np.load(path) as q:
        records.append(dict(time_s=float(q['time']),cg_m=q['cg'].tolist(),euler_deg=q['euler_deg'].tolist(),
            fit_rms_m=float(q['fit_rms']),density_range_kg_m3=[float(q['density_min']),float(q['density_max'])],
            max_speed_m_s=float(q['speed_max']),water_above_roof_particles=int(q['water_above_roof_particles']),
            water_above_upper_hull_particles=int(q['water_above_upper_hull_particles'])))
report=dict(case=a.case,frames=len(records),records=records,
    surface='3-D smoothed particle occupancy isosurface; overturning surfaces supported; no analytic waves.',
    roof_indicator='Particle count in body coordinates x=(-0.8,0.4), |y|<0.2, z=(0.49,0.9) m. Indicates water above the central roof footprint, not validated wetting or impact pressure.',
    upper_hull_indicator='Water particles 0..3dp above the highest closed-hull surface at their body-frame (x,y), where the roof height exceeds 0.2 m. Roof height is ray-intersected at 20 mm spacing; a geometrical indicator, not a wetting-pressure model.',
    particle_losses=sum(int(r['NpOut'].replace(',','')) for r in rows),
    status='completed' if times[max(times)]>=setup['duration_s']-.081 else 'partial')
(folder/'focus_motion.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k!='records'},indent=2))
