# SPDX-License-Identifier: GPL-3.0-or-later
from project_paths import ASSET_DIR, SOLVER_BIN, FONT_DIR
"""Render a measured SPH impact interval; never invent a breaking surface."""
from pathlib import Path
import argparse,hashlib,json,subprocess
from functools import lru_cache
import numpy as np
import pyvista as pv
from PIL import Image,ImageDraw,ImageFont
from scipy.spatial.transform import Rotation,Slerp
from scipy.ndimage import gaussian_filter,map_coordinates

ROOT=Path(__file__).resolve().parent;GEO=ASSET_DIR;OUT=ROOT/'media'
ap=argparse.ArgumentParser();ap.add_argument('case');ap.add_argument('--start',type=float,default=4.)
ap.add_argument('--end',type=float,default=6.);ap.add_argument('--seconds',type=float,default=3.)
ap.add_argument('--center-x',type=float,default=-.9)
ap.add_argument('--name',default='turtle_focused_intro');ap.add_argument('--stills-only',action='store_true')
a=ap.parse_args();W,H,FPS=1080,1920,30;BG='#080b12'
cache=GEO/f'{a.case}_focus_render'
files=sorted(cache.glob('frame_*.npz'))
times=[]
for path in files:
    with np.load(path) as q:times.append(float(q['time']))
times=np.array(times);assert len(times)>1
@lru_cache(maxsize=2)
def state(index):
    with np.load(files[index]) as q:return {key:q[key] for key in ['R','translation','x','y','z','render_points','render_ids']}
assert times[0]<=a.start and times[-1]>=a.end, 'Requested interval has not been simulated/cached.'
skin=np.load(GEO/'turtle_history_visual.npz');materials=np.load(GEO/'turtle_history_materials.npz')
mesh=pv.PolyData(skin['vertices'],np.c_[np.full(len(skin['faces']),3),skin['faces']].ravel())
mesh['rgb']=materials['vertex_colors']
pl=pv.Plotter(off_screen=True,window_size=(1080,1450));pl.set_background(BG)
ha=pl.add_mesh(mesh,scalars='rgb',rgb=True,smooth_shading=True,specular=.2,ambient=.28)
mesh=ha.mapper.dataset;vertices=mesh.points.copy();normals=mesh['Normals'].copy()
pl.remove_all_lights()
for pos,power,color in [((0,-3,6),1.2,'#f2dec0'),((-3,3,4),1.4,'#b0d1e1'),((5,-6,3),.4,'#c7d8dc')]:
    pl.add_light(pv.Light(position=pos,focal_point=(0,0,.2),intensity=power,color=color))
pl.enable_anti_aliasing('ssaa');pl.enable_ssao(radius=.12,bias=.003,kernel_size=32,blur=True)
pl.camera_position=[(a.center_x+4.8,-8,3.25),(a.center_x,0,-.1),(0,0,1)]
pl.camera.parallel_projection=True;pl.camera.parallel_scale=2.65
pl.add_mesh(pv.lines_from_points(np.array([[-5,-2.4,-.6],[5,-2.4,-.6],[5,2.4,-.6],[-5,2.4,-.6],[-5,-2.4,-.6]])),color='#314656',line_width=1)
pl.add_mesh(pv.Plane(center=(0,0,-.6),direction=(0,0,1),i_size=10,j_size=4.8),color='#041626')
axes=[state(0)[key] for key in ['x','y','z']];dx=float(axes[0][1]-axes[0][0]);wa=None
spray_actor=pl.add_points(pv.PolyData(np.array([[0.,0.,-20.]])),color='#86aebd',point_size=5,
    render_points_as_spheres=True,lighting=True,specular=.5)
spray_actor.SetVisibility(False)
ext=[np.r_[v[0]-dx*np.arange(3,0,-1),v,v[-1]+dx*np.arange(1,4)] for v in axes]
edges=[np.r_[v-dx/2,v[-1]+dx/2] for v in ext]
@lru_cache(maxsize=1)
def paired(index):
    s0,s1=state(index),state(index+1)
    ids,i0,i1=np.intersect1d(s0['render_ids'],s1['render_ids'],assume_unique=True,return_indices=True)
    return s0['render_points'][i0],s1['render_points'][i1]

def text(draw,pos,s,size=34,color='#8dc6d7',bold=False):
    font=ImageFont.truetype(str(FONT_DIR)+'/'+('segoeuib.ttf' if bold else 'segoeui.ttf'),size)
    draw.text(pos,s,font=font,fill=color)

def frame(t):
    global wa
    j=int(np.clip(np.searchsorted(times,t)-1,0,len(times)-2));f=float(np.clip((t-times[j])/(times[j+1]-times[j]),0,1))
    s0,s1=state(j),state(j+1)
    R=Slerp([0,1],Rotation.from_matrix(np.array([s0['R'],s1['R']])))(f).as_matrix()
    trans=(1-f)*s0['translation']+f*s1['translation']
    mesh.points[:]=vertices@R.T+trans;mesh['Normals']=normals@R.T
    # Follow the SAME fluid IDs between actual outputs. Blending two Eulerian
    # surfaces would smear a rapidly travelling or overturning crest.
    p0,p1=paired(j);points=(1-f)*p0+f*p1
    hist=np.histogramdd(points,bins=edges)[0]
    occ=gaussian_filter(hist.astype('f4'),sigma=1.0,mode='constant')[3:-3,3:-3,3:-3]
    # Keep actual dilute water parcels that a continuous isosurface omits.
    # No random foam or decorative spray is generated.
    possible=points[(points[:,2]>.12)&(abs(points[:,0])<4.8)&(abs(points[:,1])<2.2)]
    coordinates=((possible-np.array([v[0] for v in axes]))/dx).T
    support=map_coordinates(occ,coordinates,order=1,mode='constant')
    spray=possible[support<.28]
    spray_actor.SetVisibility(len(spray)>0)
    if len(spray):spray_actor.mapper.SetInputData(pv.PolyData(spray))
    padded=np.pad(occ,1);grid=pv.ImageData(dimensions=padded.shape,spacing=(dx,dx,dx),origin=[v[0]-dx for v in axes])
    grid['water']=padded.ravel(order='F')
    surface=grid.contour([.38],scalars='water').smooth(n_iter=10,relaxation_factor=.12)
    surface.compute_normals(inplace=True,split_vertices=False)
    if wa is None:wa=pl.add_mesh(surface,color='#06263c',smooth_shading=True,specular=.9,specular_power=35,ambient=.14,diffuse=.72)
    else:wa.mapper.SetInputData(surface)
    pl.render();shot=Image.fromarray(pl.screenshot(return_img=True))
    image=Image.new('RGB',(W,H),BG);image.paste(shot,(0,220));draw=ImageDraw.Draw(image)
    text(draw,(65,65),'WHEN THE WAVE HITS.',57,'#f1f3f2',True)
    text(draw,(68,150),'Turtle ship · Extreme-wave test',32)
    text(draw,(68,1730),'SPH  /  Free 6-DOF motion',29,'#91a4b1')
    text(draw,(68,1795),f't = {t:.2f} s  ·  Model scale',26,'#6d8491')
    return image

for label,t in [('before',a.start),('impact',(a.start+a.end)/2),('after',a.end)]:
    frame(t).save(OUT/f'{a.name}_{label}.png')
if a.stills_only:
    pl.close();raise SystemExit
path=OUT/f'{a.name}.mp4'
cmd=['ffmpeg','-v','error','-y','-f','rawvideo','-pix_fmt','rgb24','-s',f'{W}x{H}','-r',str(FPS),'-i','-',
    '-an','-c:v','libx264','-preset','fast','-crf','19','-pix_fmt','yuv420p','-movflags','+faststart',str(path)]
proc=subprocess.Popen(cmd,stdin=subprocess.PIPE)
count=round(a.seconds*FPS)
for i in range(count):
    t=a.start+(a.end-a.start)*i/max(count-1,1)
    proc.stdin.write(np.asarray(frame(t)).tobytes())
    if i%30==0:print('Rendered',i,'/',count,flush=True)
proc.stdin.close();code=proc.wait();pl.close();assert code==0
report=dict(video=str(path),case=a.case,physical_interval_s=[a.start,a.end],video_seconds=a.seconds,fps=FPS,frames=count,
    width=W,height=H,water='Same-ID SPH particle-position interpolation followed by a smoothed 3-D occupancy surface; no analytic replacement.',
    body='Measured rigid-particle poses applied to the licensed local turtle visual skin.',
    spray='Only dilute measured SPH water parcels omitted by the continuous surface; no random foam, no resolved droplet-size claim.',
    condition='Initially at rest, fully free six-DOF; separate extreme-wave condition from the existing coasting film.',
    fixed_camera_center_x_m=a.center_x,
    sha256=hashlib.sha256(path.read_bytes()).hexdigest())
path.with_suffix('.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
