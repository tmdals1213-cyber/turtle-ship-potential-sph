"""Synchronized portrait replay of the completed calm/head-wave GPU pilots.

Water particles are matched by their native IDs between saved outputs. Linear
particle interpolation and rigid-pose SLERP are for replay only, not new CFD.
"""
from pathlib import Path
import os
from functools import lru_cache
import argparse, csv, json, subprocess, time
import numpy as np
import pyvista as pv
import vtk
from PIL import Image, ImageDraw, ImageFont
from scipy.ndimage import gaussian_filter,map_coordinates
from scipy.spatial.transform import Rotation, Slerp
from sph_segments import part_rows,vtk_path

ROOT=Path(__file__).resolve().parent
NAMES=['gunnerus_waves_050_10s','gunnerus_calm_050_10s']
COLORS=['#64cbe7','#e3b75e']; BG='#080d16'
W,H,FPS=1080,1920,30
CACHE=ROOT/'render_cache/gunnerus_comparison'
OUT=ROOT/'media'


def prepare(name,clip_x=5.2,clip_y=3.01,z_low=-.95,z_high=1.5):
    folder=ROOT/'sph'/name; cache=CACHE/name; cache.mkdir(parents=True,exist_ok=True)
    setup=json.loads((folder/'setup.json').read_text()); p=setup['physical_model']
    state=json.loads((folder/'running.json').read_text())
    assert not state['running'] and state['returncode']==0 and not state['stop_reason']
    rows=part_rows(folder);paths={int(r['Part']):vtk_path(r) for r in rows}
    times={int(r['Part']):float(r['TimeStep [s]']) for r in rows}
    initial=pv.read(paths[0]); mask=initial['Type']==2
    order=np.argsort(initial['Idp'][mask]); ids=initial['Idp'][mask][order]
    b0=initial.points[mask][order].astype(float); m0=b0.mean(0)
    records=[]
    for part,t in times.items():
        if t<2.85: continue
        target=cache/f'frame_{part:04d}.npz'
        if not target.exists():
            d=pv.read(paths[part]); mask=d['Type']==2
            order=np.argsort(d['Idp'][mask]); assert np.array_equal(d['Idp'][mask][order],ids)
            b1=d.points[mask][order].astype(float); m1=b1.mean(0)
            u,_,vt=np.linalg.svd((b0-m0).T@(b1-m1)); rot=vt.T@u.T; assert np.linalg.det(rot)>0
            trans=m1-rot@m0; cg=rot@p['center_of_mass_m']+trans
            fit=float(np.sqrt(np.mean((b0@rot.T+trans-b1)**2))); assert fit<setup['dp_m']*.02
            mask=d['Type']==3; xyz=d.points[mask]; fids=d['Idp'][mask]
            take=(abs(xyz[:,0]-cg[0])<clip_x)&(abs(xyz[:,1])<clip_y)&(xyz[:,2]>z_low)&(xyz[:,2]<z_high)
            order=np.argsort(fids[take])
            np.savez_compressed(target,points=xyz[take][order],ids=fids[take][order],R=rot,translation=trans,cg=cg,time=t,fit_rms=fit)
        records.append(dict(part=part,time_s=t,file=target.name))
        if part%20==0: print(f'Cached {name} {t:.1f} s',flush=True)
    (cache/'index.json').write_text(json.dumps(records,indent=2))


@lru_cache(maxsize=32)
def font(size,bold=False):
    return ImageFont.truetype((os.environ.get('SHIP_FONT_DIR', 'C:/Windows/Fonts').rstrip('/\\') + '/')+('segoeuib.ttf' if bold else 'segoeui.ttf'),size)


def text(draw,xy,value,size=26,color='#e6eef3',bold=False):
    draw.text(xy,value,font=font(size,bold),fill=color)


class Scene:
    def __init__(self,name,window_size=(1080,420),parallel_scale=1.55,water_x=4.6,water_y=2.95,water_bottom=-.75,water_top=1.25,show_dilute_water=False,surface_level=.5):
        self.name=name; self.cache=CACHE/name
        self.state=lru_cache(maxsize=2)(self.state)
        self.pair=lru_cache(maxsize=1)(self.pair)
        self.index=json.loads((self.cache/'index.json').read_text())
        self.times=np.array([r['time_s'] for r in self.index])
        self.model=json.loads((ROOT/'sph'/name/'setup.json').read_text())['physical_model']
        self.parallel_scale=parallel_scale
        self.pl=pv.Plotter(off_screen=True,window_size=window_size); self.pl.set_background(BG)
        importer=vtk.vtkGLTFImporter(); importer.SetFileName(str(ROOT/'geometry/candidates/gunnerus_vesseljs/Gunnerus.glb'))
        importer.SetRenderWindow(self.pl.render_window); importer.Update(); assert importer.GetRenderer()==self.pl.renderer
        self.importer=importer; actors=importer.GetRenderer().GetActors(); actors.InitTraversal(); self.actors=[]
        visual=np.eye(4); visual[:3,:3]=np.array([[0,0,-1],[-1,0,0],[0,1,0]])/10
        visual[2,3]=self.model['visual_to_model']['vertical_shift_before_scale_m']/10
        visual[:3,3]+=np.array(self.model['visual_to_model'].get('translation_model_m',[0,0,0]))
        for _ in range(actors.GetNumberOfItems()):
            actor=actors.GetNextActor()
            old=np.array([[actor.GetMatrix().GetElement(i,j) for j in range(4)] for i in range(4)])
            self.actors.append((actor,visual@old))
            actor.SetUserTransform(None); actor.SetPosition(0,0,0); actor.SetOrientation(0,0,0); actor.SetScale(1,1,1)
            actor.GetProperty().SetMetallic(.05); actor.GetProperty().SetRoughness(.65)
        self.pl.remove_all_lights(); self.lights=[]
        for pos,intensity in [((3,-5,8),1.05),((-4,4,5),.65),((4,3,2),.4)]:
            light=pv.Light(position=pos,focal_point=(0,0,0),intensity=intensity)
            self.pl.add_light(light); self.lights.append((light,np.array(pos)))
        self.pl.camera.parallel_projection=True; self.pl.camera.parallel_scale=self.parallel_scale
        self.pl.enable_anti_aliasing('ssaa'); self.water=None
        self.dx=.05
        self.show_dilute_water=show_dilute_water;self.surface_level=surface_level;self.parcels=None
        if show_dilute_water:
            self.parcels=self.pl.add_points(pv.PolyData(np.array([[0.,0.,-99.]])),color='#8ab6c6',point_size=3,
                render_points_as_spheres=True,lighting=True,specular=.4,reset_camera=False)
            self.parcels.SetVisibility(False)
        self.axes=[np.arange(-water_x,water_x+.001,self.dx),np.arange(-water_y,water_y+.001,self.dx),np.arange(water_bottom,water_top+.001,self.dx)]
        ext=[np.r_[a[0]-self.dx*np.arange(3,0,-1),a,a[-1]+self.dx*np.arange(1,4)] for a in self.axes]
        self.edges=[np.r_[a-self.dx/2,a[-1]+self.dx/2] for a in ext]

    def state(self,index):
        with np.load(self.cache/self.index[index]['file']) as q: return {k:q[k] for k in q.files}

    def pair(self,index):
        a,b=self.state(index),self.state(index+1)
        _,i,j=np.intersect1d(a['ids'],b['ids'],assume_unique=True,return_indices=True)
        return a['points'][i],b['points'][j]

    def render(self,t):
        j=int(np.clip(np.searchsorted(self.times,t)-1,0,len(self.times)-2))
        f=float(np.clip((t-self.times[j])/(self.times[j+1]-self.times[j]),0,1))
        a,b=self.state(j),self.state(j+1)
        rot=Slerp([0,1],Rotation.from_matrix([a['R'],b['R']]))(f).as_matrix()
        trans=(1-f)*a['translation']+f*b['translation']; cg=rot@self.model['center_of_mass_m']+trans
        pose=np.eye(4); pose[:3,:3]=rot; pose[:3,3]=trans
        self.pose=pose; self.cg=cg
        for actor,original in self.actors:
            matrix=pose@original; vm=vtk.vtkMatrix4x4()
            for i in range(4):
                for k in range(4): vm.SetElement(i,k,float(matrix[i,k]))
            actor.SetUserMatrix(vm)
        p0,p1=self.pair(j); points=(1-f)*p0+f*p1; points[:,0]-=cg[0]
        occ=gaussian_filter(np.histogramdd(points,bins=self.edges)[0].astype('f4'),sigma=.8)[3:-3,3:-3,3:-3]
        if self.show_dilute_water:
            near=(abs(points[:,0])<3.3)&(abs(points[:,1]-cg[1])<1.2)&(points[:,2]>-.15)&(points[:,2]<2.5)
            possible=points[near]
            support=map_coordinates(occ,((possible-np.array([v[0] for v in self.axes]))/self.dx).T,order=1,mode='constant')
            thin=support<self.surface_level
            world=possible[thin]+[cg[0],0,0]
            self.parcels.SetVisibility(bool(len(world)))
            if len(world):self.parcels.mapper.SetInputData(pv.PolyData(world))
            vz=((p1[near,2]-p0[near,2])/(self.times[j+1]-self.times[j]))
            self.last_water_stats=dict(time_s=t,near_surface_parcels=len(possible),below_old_surface_threshold=int(sum(support<.5)),
                visible_dilute_parcels=len(world),downward_dilute_parcels=int(sum(thin&(vz<-.1))),
                velocity_note='Downward classification uses displacement of matched saved particle IDs, not an invented trajectory.')
        origin=np.array([a[0] for a in self.axes]); origin[0]+=cg[0]
        grid=pv.ImageData(dimensions=occ.shape,spacing=(self.dx,)*3,origin=origin)
        grid['water']=occ.ravel(order='F'); surface=grid.contour([self.surface_level],scalars='water')
        surface.compute_normals(inplace=True,split_vertices=False)
        if self.water is None:
            self.water=self.pl.add_mesh(surface,color='#0b3855',smooth_shading=True,diffuse=.75,specular=.45,specular_power=40)
        else: self.water.mapper.SetInputData(surface)
        for light,offset in self.lights:
            light.position=offset+[cg[0],0,0]; light.focal_point=cg
        self.pl.camera_position=[(cg[0]+6.8,-8.6,4.72),(cg[0],0,.75),(0,0,1)]
        self.pl.camera.parallel_projection=True; self.pl.camera.parallel_scale=self.parallel_scale
        self.pl.render()
        return Image.fromarray(self.pl.screenshot(return_img=True))

    def capture(self,camera,scale,size):
        """Another camera on the same updated state; no repeated water reconstruction."""
        self.pl.window_size=size; self.pl.camera_position=camera
        self.pl.camera.parallel_projection=True; self.pl.camera.parallel_scale=scale
        self.pl.render()
        return Image.fromarray(self.pl.screenshot(return_img=True))


class Film:
    def __init__(self):
        self.scenes=[Scene(n) for n in NAMES]; self.data=[]
        for name in NAMES:
            with np.load(ROOT/'sph'/name/'monitors.npz') as q: self.data.append({k:q[k] for k in q.files})
        self.cards=[('SURGE LOAD','N','force_body_N',0,1),('PITCH MOMENT','N m','moment_body_Nm',1,1),
                    ('HEAVE','mm','heave_m',None,1000),('PITCH','deg','euler_xyz_deg',1,1)]

    def chart(self,draw,xy,card,t):
        x,y=xy; title,unit,key,col,factor=card; curves=[]; vals=[]
        for d in self.data:
            v=(d[key] if col is None else d[key][:,col])*factor
            mask=(d['time_s']>=3)&(d['time_s']<=10.001)
            curves.append((d['time_s'][mask],v[mask])); vals.append(float(np.interp(t,d['time_s'],v)))
        allv=np.concatenate([q[1] for q in curves]); low,high=float(allv.min()),float(allv.max())
        pad=max((high-low)*.1,.1); low-=pad; high+=pad
        text(draw,(x,y),title,23,'#c2d5df',True); text(draw,(x+365,y),unit,21,'#8197a8')
        text(draw,(x,y+32),f'{vals[0]:+.1f}',28,COLORS[0],True)
        text(draw,(x+225,y+32),f'{vals[1]:+.1f}',28,COLORS[1],True)
        left,right,top,bottom=x,x+450,y+78,y+153
        py=lambda v:bottom-(v-low)/(high-low)*(bottom-top)
        draw.line((left,py(0),right,py(0)),fill='#233849',width=1)
        draw.line((left,bottom,right,bottom),fill='#233849',width=1)
        for (tt,v),color in zip(curves,COLORS):
            take=tt<=t
            points=list(zip((left+(tt[take]-3)/7*(right-left)).tolist(),py(v[take]).tolist()))
            if len(points)>1: draw.line(points,fill=color,width=3)
        cursor=left+(t-3)/7*(right-left); draw.line((cursor,top,cursor,bottom),fill='#758997',width=1)

    def frame(self,t):
        image=Image.new('RGB',(W,H),BG)
        for scene,y in zip(self.scenes,[198,692]): image.paste(scene.render(t),(0,y))
        draw=ImageDraw.Draw(image)
        text(draw,(60,42),'GUNNERUS',48,'#f0f5f7',True)
        text(draw,(795,57),'1:10 MODEL',25,'#9fb7c8')
        text(draw,(60,109),'SAME THRUST  /  DIFFERENT SEAS',26,'#9fb7c8')
        text(draw,(800,109),f't  {t:05.2f} s',25,'#d0dfe8')
        for i,(label,y) in enumerate([('HEAD WAVES',165),('CALM WATER',659)]):
            text(draw,(60,y),label,28,COLORS[i],True)
            d=self.data[i]; speed=np.interp(t,d['time_s'],d['speed_body_m_s'])
            text(draw,(735,y),f'U  {speed:+.2f} m/s',25,COLORS[i])
        text(draw,(60,615),'SHIP  >                                      <  WAVES',21,'#8197a8')
        for xy,card in zip([(60,1150),(570,1150),(60,1330),(570,1330)],self.cards): self.chart(draw,xy,card,t)
        text(draw,(60,1550),f'MODEL TIME  {t:05.2f} s',28,'#d0dfe8',True)
        text(draw,(60,1600),'1:10 model  /  30 N assumed thrust  /  Free 6 DOF',23,'#849cad')
        text(draw,(60,1641),'Pitch: bow-down positive  /  Loads: body axes',21,'#849cad')
        text(draw,(60,1693),'Actual SPH replay  /  No propeller wake',23,'#849cad')
        text(draw,(60,1780),'Visual model: IHB NTNU / vessel.js | MIT',21,'#687e8f')
        return image


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--prepare',action='store_true'); ap.add_argument('--stills-only',action='store_true')
    args=ap.parse_args()
    if args.prepare:
        for name in NAMES: prepare(name)
        return
    film=Film(); OUT.mkdir(exist_ok=True)
    for t in [3.,6.8,8.6,9.5]:
        target=OUT/f'gunnerus_comparison_t{t:.1f}.png'; film.frame(t).save(target); print(target,flush=True)
    if args.stills_only: return
    path=OUT/'gunnerus_calm_vs_headwaves_portrait.mp4'; seconds=14.; n=int(seconds*FPS)
    command=['ffmpeg','-v','error','-y','-f','rawvideo','-pix_fmt','rgb24','-s',f'{W}x{H}',
        '-r',str(FPS),'-i','-','-an','-c:v','libx264','-preset','fast','-crf','19','-pix_fmt','yuv420p','-movflags','+faststart',str(path)]
    started=time.perf_counter(); encoder=subprocess.Popen(command,stdin=subprocess.PIPE)
    try:
        for i,t in enumerate(np.linspace(3.,10.,n)):
            encoder.stdin.write(np.asarray(film.frame(float(t))).tobytes())
            if i%30==0: print(f'Rendered {i}/{n}, model t={t:.2f} s',flush=True)
    finally:
        encoder.stdin.close()
    assert encoder.wait()==0
    meta=dict(video=str(path),video_seconds=seconds,fps=FPS,frames=n,dimensions=[W,H],model_interval_s=[3.,10.],
        cases=NAMES,render_wall_s=time.perf_counter()-started,
        time_semantics='3..10 model seconds shown over 14 video seconds; no claim of real-time CFD.',
        water='Actual water particle positions matched by ID and linearly interpolated between 0.1 s solver outputs; occupancy isosurface 0.5 with sigma 0.8 dp.',
        pose='SLERP of rigid fits to actual hull particles; rendered superstructure remains display-only.',
        monitors='Native 0.01 s solver force/moment/CG samples; body projections use saved rigid orientations. No load filtering.',
        presentation='Common fixed chart scales for both cases over the stated 3..10 s interval. Complete 0..10 s records remain in the static comparison.',
        limitations=['Assumed external force at CG, not actuator disk or resolved propeller.','No validated slam pressure or propulsive efficiency.'])
    path.with_suffix('.json').write_text(json.dumps(meta,indent=2)); print(json.dumps(meta,indent=2),flush=True)


if __name__=='__main__': main()
