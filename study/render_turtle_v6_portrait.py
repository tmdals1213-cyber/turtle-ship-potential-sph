# SPDX-License-Identifier: GPL-3.0-or-later
from project_paths import ASSET_DIR, SOLVER_BIN, FONT_DIR
"""Portrait film: lower SPH camera plus a synchronized fixed side view.

Simulation data are unchanged. Source credit and detailed limits accompany the
video in a text description rather than covering every rendered scene.
"""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
from scipy.spatial.transform import Rotation, Slerp

parser=argparse.ArgumentParser()
parser.add_argument('--stills-only',action='store_true')
args=parser.parse_args()
saved=sys.argv;sys.argv=[saved[0],'--library-only']
import render_turtle_v3 as s
sys.argv=saved
ROOT=s.ROOT;out=s.out;geo=s.geo
W,H,FPS,DURATION=1080,1920,30,32
BG='#080b12';WHITE='#f1f3f2';BLUE='#8dc6d7';GOLD='#c9a274';GRAY='#8fa1ae'

materials=np.load(geo/'turtle_history_materials.npz')
s.hull.point_data['material_rgb']=materials['vertex_colors']
s.ha.mapper.SetScalarModeToUsePointFieldData();s.ha.mapper.SelectColorArray('material_rgb');s.ha.mapper.SetColorModeToDirectScalars()
s.ha.GetProperty().SetAmbient(.27);s.ha.GetProperty().SetDiffuse(.86);s.ha.GetProperty().SetSpecular(.16);s.ha.GetProperty().SetSpecularPower(22)
s.pl.remove_all_lights()
for position,intensity,color in [((0,-3,6),1.1,'#f2dec0'),((-3,3,2.5),1.2,'#b0d1e1'),((4,-5,2),.35,'#c7d8dc')]:
    s.pl.add_light(pv.Light(position=position,focal_point=(0,0,.2),intensity=intensity,color=color))
s.wateractor.GetProperty().SetColor(pv.Color('#0a3454').float_rgb)
s.wateractor.GetProperty().SetSpecular(.9);s.wateractor.GetProperty().SetSpecularPower(35)
s.fa.GetProperty().SetColor(pv.Color('#062c38').float_rgb);s.sa.GetProperty().SetColor(pv.Color('#062c38').float_rgb)
s.pl.enable_ssao(radius=.12,bias=.003,kernel_size=32,blur=True)
old_frame=[a for a in s.pl.actors.values() if a not in [s.ha,s.wateractor,s.fa,s.sa,s.ea]]

fixed=[np.load(p) for p in sorted((geo/'free6_history_40_portrait_fixed').glob('frame_*.npz'))]
ft=np.array([float(q['time']) for q in fixed]);origin=fixed[0]['origin']
assert len(fixed)==len(s.states)
axes=[fixed[0][key] for key in ['x','y','z']];dx=float(axes[0][1]-axes[0][0])
initial_cg=fixed[0]['cg'];xstart=float(initial_cg[0]-origin[0])
fixed_frame=[s.pl.add_mesh(pv.Box(bounds=(-3.2,3.2,-1.2,1.2,-.4,.72)).outline(),color='#526675',line_width=1.5),
    s.pl.add_mesh(pv.Plane(center=(0,0,-.4),direction=(0,0,1),i_size=6.4,j_size=2.4),color='#041a30'),
    s.pl.add_mesh(pv.Line((xstart,-1.25,-.4),(xstart,-1.25,.12)),color=GOLD,line_width=3)]
for x in np.arange(-5,1)-origin[0]:
    fixed_frame.append(s.pl.add_mesh(pv.Line((x,-1.21,-.4),(x,-1.21,-.31)),color='#6c8795',line_width=2))
side_edge=pv.lines_from_points(np.c_[axes[0],np.full(len(axes[0]),-1.215),np.zeros(len(axes[0]))])
side_edge_actor=s.pl.add_mesh(side_edge,color=BLUE,line_width=2)
side_edge_actor.SetVisibility(False)

def text(draw,xy,message,size=30,color=GRAY,bold=False):
    font='segoeuib.ttf' if bold else 'segoeui.ttf'
    draw.text(xy,message,font=ImageFont.truetype(str(FONT_DIR)+'/'+font,size),fill=color)

def arrow(draw,p,q,color=BLUE,width=4):
    draw.line([p,q],fill=color,width=width)
    v=np.asarray(q,float)-p;v/=np.linalg.norm(v);n=np.array([-v[1],v[0]]);end=np.asarray(q,float)
    draw.polygon([tuple(end),tuple(end-15*v+7*n),tuple(end-15*v-7*n)],fill=color)

def page(title,subtitle=''):
    image=Image.new('RGB',(W,H),BG);draw=ImageDraw.Draw(image)
    text(draw,(65,65),title,62,WHITE,True)
    if subtitle:text(draw,(68,150),subtitle,32,BLUE)
    return image,draw

def mode(which):
    side_edge_actor.SetVisibility(False)
    for actor in old_frame:actor.SetVisibility(which=='potential')
    for actor in fixed_frame:actor.SetVisibility(which=='fixed')
    if which!='potential':
        for actor in [s.fa,s.sa,s.ea]:actor.SetVisibility(False)

def fixed_sph(t):
    mode('fixed')
    s.wateractor.GetProperty().SetOpacity(1)
    s.pl.window_size=(1080,810)
    s.pl.camera_position=[(3.5,-6,3.0),(0,0,.15),(0,0,1)]
    s.pl.camera.parallel_projection=True;s.pl.camera.parallel_scale=2.78
    j=int(np.clip(np.searchsorted(ft,t)-1,0,len(ft)-2));f=float(np.clip((t-ft[j])/(ft[j+1]-ft[j]),0,1))
    a,b=fixed[j],fixed[j+1]
    R=Slerp([0,1],Rotation.from_matrix(np.array([a['R'],b['R']])))(f).as_matrix()
    translation=(1-f)*a['translation']+f*b['translation']
    s.pose(R,R@s.offset+translation-origin)
    occ=(1-f)*a['occupancy'].astype('f4')+f*b['occupancy'].astype('f4')
    padded=np.pad(occ,1)
    grid=pv.ImageData(dimensions=padded.shape,spacing=(dx,dx,dx),origin=[v[0]-dx for v in axes])
    grid['water']=padded.ravel(order='F')
    surface=grid.contour([.45],scalars='water').smooth(n_iter=8,relaxation_factor=.12)
    surface.compute_normals(inplace=True,split_vertices=False);s.setwater(surface)
    center=(1-f)*a['cg']+f*b['cg']
    s.pl.render()
    s.pl.renderer.SetWorldPoint(xstart,-1.25,.12,1);s.pl.renderer.WorldToDisplay()
    px,py,_=s.pl.renderer.GetDisplayPoint()
    main=Image.fromarray(s.pl.screenshot(return_img=True))
    # Both views use this very same rigid pose and water surface. Only the
    # display camera changes. Geometry-scope colour coding belongs to the
    # methods page, not to this motion view.
    for actor in fixed_frame:actor.SetVisibility(False)
    s.pl.window_size=(940,410)
    s.pl.camera_position=[(-.1,-8,.68),(-.1,0,.68),(0,0,1)]
    s.pl.camera.parallel_projection=True;s.pl.camera.parallel_scale=1.22
    s.wateractor.GetProperty().SetOpacity(.78)
    # Near-side free-surface profile recovered from the same measured occupancy.
    zz=[]
    for column in occ[:,4,:]:
        indices=np.where(column>.45)[0]
        kk=int(indices[-1]) if len(indices) else 0
        z=float(axes[2][kk])
        if kk+1<len(column) and column[kk]!=column[kk+1]:
            z+=dx*(column[kk]-.45)/(column[kk]-column[kk+1])
        zz.append(z)
    side_edge.points[:,2]=zz;side_edge_actor.SetVisibility(True)
    s.pl.render();side=Image.fromarray(s.pl.screenshot(return_img=True))
    s.wateractor.GetProperty().SetOpacity(1);side_edge_actor.SetVisibility(False)
    return main,side,center,(px,810-py)

def hero(t):
    mode('hero');s.pl.window_size=(1080,1250)
    angle=np.radians(-61+4*t/3)
    s.pl.camera_position=[(6*np.cos(angle),6*np.sin(angle),2.7),(0,0,.38),(0,0,1)]
    s.pl.camera.parallel_scale=2.2
    shot,*_=s.sph(10.5+.6*t)
    image=Image.new('RGB',(W,H),BG);image.paste(shot,(0,250));draw=ImageDraw.Draw(image)
    text(draw,(65,70),'CAN WE PREDICT',60,WHITE,True)
    text(draw,(65,145),'THIS SHIP IN WAVES?',58,WHITE,True)
    text(draw,(70,1480),'Waves. Weight. Buoyancy.',36,BLUE)
    return image

# Explanation graphics retain the entire simplified upper boundary.
closed=pv.read(geo/'turtle_history_closed.stl')
display=s.hull.copy(deep=True);display.points=s.hb.copy();display.clear_data()
source=np.load(geo/'turtle_source_mesh.npz')
envelope_names=['Cube.000','Cube.003','Cube.004','Cube.006','Cube.008','Cube.009','Cube.001','Cube.005']
envelope_ids=[i for i,name in enumerate(source['names']) if name in envelope_names]
detail_faces=s.skin['faces'][~np.isin(source['object_ids'],envelope_ids)]
details=pv.PolyData(s.skin['vertices'],np.c_[np.full(len(detail_faces),3),detail_faces].ravel())
pp=pv.Plotter(off_screen=True,window_size=(960,490));pp.set_background(BG)
pp.add_mesh(display,color='#81909a',opacity=.12,lighting=False)
pp.add_mesh(s.pm,color=BLUE,show_edges=True,edge_color='#315465',line_width=.65,lighting=False)
pp.add_mesh(pv.Line((-1.75,0,0),(1.25,0,0)),color='#dde6e8',line_width=2)
pp.camera_position=[(2.5,-4,1.7),(-.22,0,.65),(0,0,1)]
pp.camera.parallel_projection=True;pp.camera.parallel_scale=1.09
panel=Image.fromarray(pp.screenshot(return_img=True));pp.close()

vm=pv.read(ROOT/'sph/free6_history_40/generated/Case_All.vtk');points=vm.points-s.offset;kind=vm['Type']
body=points[kind==2]
water=points[(kind==3)&(points[:,0]>-1.9)&(points[:,0]<1.3)&(abs(points[:,1])<.58)&(points[:,2]>-.25)]
pp=pv.Plotter(off_screen=True,window_size=(960,490));pp.set_background(BG)
pp.add_mesh(details,color='#81909a',opacity=.20,lighting=False)
pp.add_points(water,color='#377892',point_size=3.5,opacity=.65,render_points_as_spheres=True)
pp.add_points(body,color=GOLD,point_size=4.7,render_points_as_spheres=True,lighting=False)
pp.camera_position=[(2.5,-4,1.7),(-.22,0,.65),(0,0,1)]
pp.camera.parallel_projection=True;pp.camera.parallel_scale=1.09
particles=Image.fromarray(pp.screenshot(return_img=True));pp.close()

theory,draw=page('ONE SHIP.','Two independent methods')
text(draw,(65,300),'POTENTIAL FLOW',40,WHITE,True)
text(draw,(67,355),'Blue: computed wetted panels',29,BLUE)
theory.paste(panel,(60,405));draw=ImageDraw.Draw(theory)
text(draw,(65,970),'SPH',44,WHITE,True)
text(draw,(67,1030),'Gold: computed hull boundary',29,GOLD)
theory.paste(particles,(60,1080));draw=ImageDraw.Draw(theory)
text(draw,(67,1635),'Grey detail: display only',31,GRAY)

def potential_frame(t):
    mode('potential');s.pl.window_size=(1080,1080)
    s.pl.camera_position=[(3.5,-6,3.4),(0,0,.2),(0,0,1)]
    s.pl.camera.parallel_scale=2.98
    shot,motion,z=s.potential(t)
    image,draw=page('POTENTIAL FLOW','Linear wave response')
    image.paste(shot,(0,280));draw=ImageDraw.Draw(image)
    text(draw,(68,230),'SHIP-FIXED VIEW',25,GRAY)
    text(draw,(68,1390),'WAVE PROFILE',25,BLUE,True)
    xp=70+(s.xs-s.xs[0])/(s.xs[-1]-s.xs[0])*940;yp=1510-z*310
    draw.polygon([(70,1605)]+list(zip(xp,yp))+[(1010,1605)],fill='#0a2c3d')
    draw.line(list(zip(xp,yp)),fill=BLUE,width=4)
    draw.line([(70,1510),(1010,1510)],fill='#34545f',width=1)
    text(draw,(70,1745),'154 s · CPU',31,GRAY)
    text(draw,(70,1800),'Ryzen 9 9950X',24,GRAY)
    return image

def sph_frame(t):
    shot,side,center,marker=fixed_sph(t)
    image,draw=page('SPH','Coasting into oncoming waves')
    image.paste(shot,(0,270));draw=ImageDraw.Draw(image)
    text(draw,(70,230),'SHIP',25,GOLD,True);arrow(draw,(150,247),(300,247),GOLD)
    text(draw,(810,230),'WAVES',25,BLUE,True);arrow(draw,(795,247),(645,247),BLUE)
    mx,my=marker;text(draw,(mx-34,my+270-31),'START',22,GOLD,True)
    text(draw,(70,1050),'FIXED CAMERA · NO PROPULSION',27,BLUE)
    text(draw,(70,1155),'SIDE VIEW',29,WHITE,True)
    image.paste(side,(70,1210));draw=ImageDraw.Draw(image)
    text(draw,(70,1745),'11.5 min · GPU',31,GRAY)
    text(draw,(70,1800),'RTX 5080',24,GRAY)
    return image

def frame(t):
    if t<3:return hero(t)
    if t<9:return theory.copy()
    if t<17:return potential_frame(t-9)
    elapsed=t-17
    physical=elapsed/5*8 if elapsed<5 else 8+(elapsed-5)/10*(ft[-1]-8)
    return sph_frame(float(min(physical,ft[-1])))

for moment,name in [(1.5,'opening'),(6,'methods'),(11.5,'potential'),(17,'sph_start'),(21.5,'sph_forward'),(30.5,'sph_waves')]:
    frame(moment).save(out/f'turtle_v6_portrait_{name}.png')
if args.stills_only:
    s.pl.close();raise SystemExit()

video=out/'turtle_ship_story_v6_portrait.mp4'
cmd=[shutil.which('ffmpeg'),'-y','-f','rawvideo','-pix_fmt','rgb24','-s',f'{W}x{H}','-r',str(FPS),'-i','-','-an','-c:v','libx264','-preset','fast','-crf','18','-pix_fmt','yuv420p','-movflags','+faststart',str(video)]
with (out/'turtle_v6_portrait_render.log').open('w') as log:
    process=subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=log,stderr=log)
    for i in range(DURATION*FPS):
        process.stdin.write(np.asarray(frame(i/FPS)).tobytes())
        if i%120==0:print(f'Portrait frame {i} / {DURATION*FPS}',flush=True)
    process.stdin.close();assert process.wait()==0
s.pl.close()
meta=dict(duration_s=DURATION,size=[W,H],fps=FPS,source_simulation='free6_history_40',
    sections=[['question',0,3],['methods together',3,9],['potential result',9,17],['SPH result',17,32]],
    changes=['Lower 22 degree SPH camera, framed slightly farther back','Synchronized fixed side view with the display model and actual water','Computed versus visual geometry distinguished on the mesh/particle methods page','Main-view start marker retained; position bar removed for the side view','No artificial translation or changes to the simulation'],
    sph_time_mapping='First 5 video seconds: 0..8 model seconds; next 10 seconds: 8..15.900095 model seconds',
    sph_frame='Two fixed laboratory cameras; identical interpolated rigid pose and particle surface. No propulsion.',
    sph_main_camera=dict(position=[3.5,-6,3.0],target=[0,0,.15],parallel_scale=2.78,viewport=[1080,810]),
    sph_side_camera=dict(position=[-.1,-8,.68],target=[-.1,0,.68],parallel_scale=1.22,viewport=[940,410],water_opacity=.78),
    geometry_scope_display='Methods page only: blue wetted BEM panels; gold full SPH hull boundary; grey display details for context. Side view is a motion view, not a geometry-scope plot.',
    side_waterline='Near-side profile recovered from the same SPH occupancy, drawn in blue; not an analytical wave.',
    potential_frame='Ship-fixed view of a linear harmonic result at fixed reference mean speed',
    upper_boundary=dict(full_floating_particles=len(body),above_initial_waterline=int(sum(body[:,2]>0)),mesh='Closed simplified hull, not every decorative appendage'),
    physics='Unchanged calculations. Water-only free surface; no air flow, flooding or appendage hydrodynamics.',
    source_credit_policy='Creator page permits commercial rendered use and forbids mesh redistribution/sale; no explicit mandatory in-video attribution wording found on 2026-11-18.',
    sha256=hashlib.sha256(video.read_bytes()).hexdigest())
(out/(video.stem+'.json')).write_text(json.dumps(meta,indent=2))
print(video)
