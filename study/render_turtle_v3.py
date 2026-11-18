# SPDX-License-Identifier: GPL-3.0-or-later
from project_paths import ASSET_DIR, SOLVER_BIN, FONT_DIR
"""Five concise scenes: physical model, potential setup/result, SPH setup/result."""
from pathlib import Path
import argparse,json,sys,os,subprocess,shutil,hashlib
import numpy as np
import pyvista as pv
from PIL import Image,ImageDraw,ImageFont
from scipy.interpolate import RegularGridInterpolator,NearestNDInterpolator
from scipy.spatial.transform import Rotation,Slerp
from matplotlib.path import Path as PolygonPath
os.environ['MPLBACKEND']='Agg'
ROOT=Path(__file__).resolve().parent;geo=ASSET_DIR;out=ROOT/'media'
ap=argparse.ArgumentParser();ap.add_argument('--stills-only',action='store_true');ap.add_argument('--library-only',action='store_true');args=ap.parse_args()
model=json.loads((ROOT/'results/turtle_history_model.json').read_text(encoding='utf-8'));cg=np.array(model['center_of_mass_m'])
base=ROOT/'results/potential_turtle_history';co=np.load(base/'baseline.npz');xi=co['xi']*.12;we=float(co['encounter_omega']);k=2*np.pi/3
skin=np.load(geo/'turtle_history_visual.npz');field=np.load(base/'surface.npz');sx,sy=np.meshgrid(field['x'],field['y'],indexing='ij');valid=field['valid'];vals=field['eta_scattered'].copy()
envelope=np.load(geo/'turtle_extended_half_breadth.npz');breadth=RegularGridInterpolator((envelope['x'],envelope['z']),envelope['half_breadth'],bounds_error=False,fill_value=0)
trim=Rotation.from_euler('y',model['initial_trim_deg'],degrees=True).as_matrix();cgkeel=np.array(model['cg_before_equilibrium_m'])
vals[~valid]=NearestNDInterpolator(np.c_[sx[valid],sy[valid]],vals[valid])(sx[~valid],sy[~valid])
scatter=RegularGridInterpolator((field['x'],field['y']),vals,bounds_error=True)
a2=k*.12**2*(3-np.tanh(k)**2)/(4*np.tanh(k)**3)
def elevation(x,y,t):
    x,y=np.broadcast_arrays(x,y);phase=-k*x-we*t
    return .12*np.cos(phase)+a2*np.cos(2*phase)+.12*np.real(scatter(np.c_[x.ravel(),y.ravel()]).reshape(x.shape)*np.exp(-1j*we*t))

BG='#080b12';WHITE='#f3f5f7';BLUE='#8dd5ed';GRAY='#90a4b9'
def text(d,xy,s,n=22,col=GRAY,bold=False):
    fn='malgun.ttf' if any('\uac00'<=x<='\ud7a3' for x in s) else ('seguisym.ttf' if '∇' in s else ('segoeuib.ttf' if bold else 'segoeui.ttf'))
    d.text(xy,s,font=ImageFont.truetype(str(FONT_DIR)+'/'+fn,n),fill=col)
def arrow(d,p,q,color=BLUE,width=3):
    d.line([p,q],fill=color,width=width);v=np.array(q,float)-p;v/=np.linalg.norm(v);n=np.array([-v[1],v[0]]);q=np.array(q)
    d.polygon([tuple(q),tuple(q-12*v+5*n),tuple(q-12*v-5*n)],fill=color)
def page(title,subtitle,number):
    im=Image.new('RGB',(1440,810),BG);d=ImageDraw.Draw(im)
    text(d,(40,26),title,35,WHITE,True);text(d,(40,78),subtitle,23,BLUE)
    for j in range(5):d.rounded_rectangle((1260+j*28,42,1277+j*28,48),radius=3,fill=BLUE if j==number else '#26374a')
    return im,d
def footer(d,left,right=''):
    text(d,(40,758),left,18,GRAY);text(d,(40,786),'Turtle model: 예은이네 / 와이프잘모시기  ·  Capytaine / DualSPHysics',13,'#647b93')
    if right:text(d,(1060,786),right,13,'#647b93')
def style(pl):
    pl.set_background(BG);pl.remove_all_lights()
    pl.add_light(pv.Light(position=(-3,3,2.5),focal_point=(0,0,0),intensity=1.1,color='#d4e8ff'))
    pl.add_light(pv.Light(position=(2,-5,5),focal_point=(0,0,0),intensity=.6,color='#fff0d7'))
    pl.camera_position=[(3.5,-6,3.05),(0,0,-.02),(0,0,1)];pl.camera.parallel_projection=True;pl.camera.parallel_scale=1.91
    pl.enable_anti_aliasing('ssaa')
pl=pv.Plotter(off_screen=True,window_size=(1080,610));style(pl)
hull=pv.PolyData(skin['vertices'],np.c_[np.full(len(skin['faces']),3),skin['faces']].ravel());hull['colors']=skin['colors']
ha=pl.add_mesh(hull,scalars='colors',rgb=True,smooth_shading=True,specular=.2,ambient=.3);hull=ha.mapper.dataset;hb=hull.points.copy();hn=hull['Normals'].copy()
pl.add_mesh(pv.Box(bounds=(-2.8,2.8,-1.2,1.2,-.4,.65)).outline(),color='#5f7286',line_width=1.3)
pl.add_mesh(pv.Plane(center=(0,0,-.4),direction=(0,0,1),i_size=5.6,j_size=2.4),color='#041a30')
wateractor=None

xx,yy=np.meshgrid(np.linspace(-2.8,2.8,181),np.linspace(-1.2,1.2,81),indexing='ij');gr=pv.StructuredGrid(xx,yy,np.zeros_like(xx));centers=gr.cell_centers().points
keep=~PolygonPath(skin['waterplane']).contains_points(centers[:,:2],radius=.035)
potwater=gr.extract_cells(keep).extract_surface();pw=potwater.points.copy()
xs=np.linspace(-2.8,2.8,181);ys=np.linspace(-1.2,1.2,81)
def curtain(x,y,z):
    n=len(x);v=np.vstack([np.c_[x,y,z],np.c_[x,y,np.full(n,-.4)]]);faces=np.array([[4,i,i+1,n+i+1,n+i] for i in range(n-1)])
    return pv.PolyData(v,faces.ravel())
front=curtain(xs,np.full(len(xs),-1.2),np.zeros(len(xs)));side=curtain(np.full(len(ys),2.8),ys,np.zeros(len(ys)))
fa=pl.add_mesh(front,color='#06345a',opacity=.96);sa=pl.add_mesh(side,color='#06345a',opacity=.96)
edge=pv.lines_from_points(np.c_[xs,np.full(len(xs),-1.205),np.zeros(len(xs))]);ea=pl.add_mesh(edge,color=BLUE,line_width=2)
def setwater(mesh):
    global wateractor
    if wateractor is None:wateractor=pl.add_mesh(mesh,color='#075181',smooth_shading=True,specular=.8,specular_power=45,ambient=.1,diffuse=.85)
    else:wateractor.mapper.SetInputData(mesh)
def pose(R,translation):
    hull.points[:]=hb@R.T+translation;hull['Normals']=hn@R.T
def potential(t):
    global potwater
    motion=np.real(xi*np.exp(-1j*we*t));R=Rotation.from_euler('xyz',motion[3:]).as_matrix();pose(R,cg+motion[:3]-R@cg)
    # The linear field is evaluated about the mean hull, but its display mask
    # follows the actual visual pose so a fixed waterline hole is not left behind.
    points=gr.points;points[:,2]=elevation(points[:,0],points[:,1],t);gr.points=points
    centers=gr.cell_centers().points;local=(centers-(cg+motion[:3]-R@cg))@R
    original=((local-cgkeel-[0,0,model['equilibrium_shift_z_m']])@trim+cgkeel)/model['original_to_model_scale']+[0,0,model['original_keel_z']]
    b=breadth(original[:,[0,2]]);keep=(b==0)|(abs(original[:,1])>b+.012/model['original_to_model_scale'])
    potwater=gr.extract_cells(keep).extract_surface()
    potwater.compute_normals(inplace=True,split_vertices=False);setwater(potwater)
    z=elevation(xs,np.full(len(xs),-1.2),t);front.points[:len(xs),2]=z;edge.points[:,2]=z;side.points[:len(ys),2]=elevation(np.full(len(ys),2.8),ys,t)
    for ac in [fa,sa,ea]:ac.SetVisibility(True)
    pl.render();return Image.fromarray(pl.screenshot(return_img=True)),motion,z

# Actual BEM panels, shown separately from the decorative visual skin.
pan=np.load(geo/'turtle_history_solution_panels.npz');tri=pan['faces'][:,:3];pm=pv.PolyData(pan['vertices'],np.c_[np.full(len(tri),3),tri].ravel())
pp=pv.Plotter(off_screen=True,window_size=(830,470));style(pp);pp.add_mesh(pm,color='#8da5b9',show_edges=True,edge_color='#28435a',line_width=.6,lighting=False)
ii=np.arange(0,len(pan['centers']),max(1,len(pan['centers'])//45));sign=pan['sources'].real[ii]>0
for selected,col in [(ii[sign],BLUE),(ii[~sign],'#ffb783')]:
    pp.add_points(pan['centers'][selected],point_size=8,color=col,render_points_as_spheres=True,lighting=False)
pp.camera_position=[(2.8,-4,-1.7),(0,0,-.08),(0,0,1)];pp.camera.parallel_scale=.90
panelimg=Image.fromarray(pp.screenshot(return_img=True));pp.close()

cache=geo/'free6_history_40_render';files=sorted(cache.glob('frame_*.npz'));states=[np.load(p) for p in files];st=np.array([float(q['time']) for q in states])
assert len(st)>1
motiondata=json.loads((ROOT/'sph/free6_history_40/motion_6dof.json').read_text(encoding='utf-8'));mt=np.array(motiondata['times_s']);mv=np.array(motiondata['velocity_m_s'])
offset=np.array([-3,0,0]);gridaxes=[states[0][axis] for axis in ['x','y','z']];dx=float(gridaxes[0][1]-gridaxes[0][0])
def sph(t):
    j=int(np.clip(np.searchsorted(st,t)-1,0,len(st)-2));f=float(np.clip((t-st[j])/(st[j+1]-st[j]),0,1));a,b=states[j],states[j+1]
    R=Slerp([0,1],Rotation.from_matrix(np.array([a['R'],b['R']])))(f).as_matrix();trans=(1-f)*a['translation']+f*b['translation'];follow=(1-f)*a['follow']+f*b['follow']
    pose(R,R@offset+trans-follow)
    occ=(1-f)*a['occupancy'].astype('f4')+f*b['occupancy'].astype('f4');padded=np.pad(occ,1)
    grid=pv.ImageData(dimensions=padded.shape,spacing=(dx,dx,dx),origin=[v[0]-dx for v in gridaxes]);grid['water']=padded.ravel(order='F')
    surf=grid.contour([.45],scalars='water').smooth(n_iter=8,relaxation_factor=.12);surf.compute_normals(inplace=True,split_vertices=False);setwater(surf)
    for ac in [fa,sa,ea]:ac.SetVisibility(False)
    zz=[]
    for col in occ[:,4,:]:
        idx=np.where(col>.45)[0];kk=int(idx[-1]) if len(idx) else 0
        z=gridaxes[2][kk]
        if kk+1<len(col) and col[kk]!=col[kk+1]:z+=dx*(col[kk]-.45)/(col[kk]-col[kk+1])
        zz.append(z)
    cgworld=(1-f)*a['cg']+f*b['cg'];angles=Rotation.from_matrix(R).as_euler('xyz',degrees=True)
    speed=float(np.interp(t,mt,mv[:,0]));pl.render()
    return Image.fromarray(pl.screenshot(return_img=True)),cgworld,angles,speed,np.array(zz)

def profile(d,x,z,y=240):
    ox,W=1110,285;oy=y+75;xp=ox+(x-x[0])/(x[-1]-x[0])*W;yp=oy-z*380
    text(d,(ox,y-20),'SIDE WAVE PROFILE',18,BLUE,True)
    d.polygon([(ox,y+155)]+list(zip(xp,yp))+[(ox+W,y+155)],fill='#09395f');d.line(list(zip(xp,yp)),fill=BLUE,width=3)
    d.line([(ox,oy),(ox+W,oy)],fill='#486177');text(d,(ox,y+166),'Vertical scale enlarged',14,GRAY)

still,m,z=potential(0)
pose(np.eye(3),np.zeros(3));gr.points[:,2]=0;centers=gr.cell_centers().points;keep=~PolygonPath(skin['waterplane']).contains_points(centers[:,:2],radius=.025)
potwater=gr.extract_cells(keep).extract_surface();front.points[:len(xs),2]=0;side.points[:len(ys),2]=0;edge.points[:,2]=0
for attrs in [potwater.point_data,potwater.cell_data]:
    if 'Normals' in attrs:del attrs['Normals']
potwater.compute_normals(inplace=True,split_vertices=False);setwater(potwater);pl.render();still=Image.fromarray(pl.screenshot(return_img=True))
pl.renderer.SetWorldPoint(*cg,1);pl.renderer.WorldToDisplay();gx,gy,_=pl.renderer.GetDisplayPoint();gy=130+610-gy
pl.renderer.SetWorldPoint(*model['center_of_buoyancy_m'],1);pl.renderer.WorldToDisplay();bx,by,_=pl.renderer.GetDisplayPoint();by=130+610-by
physical,d=page('A TURTLE SHIP IN WAVES','History-informed loading  ·  Initial speed, then free motion',0);physical.paste(still,(0,130));d=ImageDraw.Draw(physical)
text(d,(1095,200),'154.4 t',44,WHITE,True);text(d,(1098,255),'Estimated loaded mass',21)
text(d,(1095,320),'H 2.4 m',31,WHITE,True);text(d,(1095,373),'Initial 1.58 m/s',27,WHITE)
text(d,(1095,435),'6 DOF FREE',28,BLUE,True);text(d,(1095,490),'Gravity + buoyancy',22)
arrow(d,(gx,gy-150),(gx,gy),'#ffb783',5);text(d,(gx+22,gy-135),'Weight',23,'#ffb783')
arrow(d,(bx,by+150),(bx,by),BLUE,5);text(d,(bx-145,by+100),'Buoyancy',23,BLUE)
text(d,(64,171),'AIR',21,GRAY,True);text(d,(92,654),'WATER',22,BLUE,True)
footer(d,'Full-scale equivalents  |  1:10 simulation  |  1795 records / 2026 reconstruction reference')
text(d,(70,704),'CROPPED WATER VOLUME',15,GRAY)

theory,d=page('POTENTIAL FLOW','A separate linear model at a reference speed',1);theory.paste(panelimg,(20,225));d=ImageDraw.Draw(theory)
text(d,(70,163),'FIT THE HULL → PANEL THE SURFACE',24,BLUE,True)
text(d,(900,207),'u = ∇φ     ∇²φ = 0',31,WHITE,True);text(d,(900,282),'No flow through the hull',24)
text(d,(900,357),'Solve source strengths',25,WHITE,True);text(d,(900,415),'Pressure → 6-DOF response',25,BLUE)
text(d,(900,490),'Incoming + scattered waves',23);text(d,(900,536),'Motion radiates waves',23)
for i,label in enumerate(['Waves','Surface pressure','Rigid-body response']):
    x=60+i*455;d.rounded_rectangle((x,668,x+385,725),radius=8,outline='#3c627a',width=2);text(d,(x+22,682),label,23,WHITE,True)
    if i<2:arrow(d,(x+395,696),(x+445,696))
footer(d,'Linear, inviscid response  |  Sources are mathematical, not water injection')

sphsetup,d=page('SPH: FREE WATER + FREE BODY','Initial velocity only  ·  No imposed ship speed',3)
surface=[(75+x,275+7*np.cos(x/32)) for x in range(1290)]
d.polygon([(75,415)]+surface+[(1364,415)],fill='#08273f');d.line(surface,fill=BLUE,width=3)
d.rectangle((65,201,1370,418),outline='#6c8299',width=2)
d.rectangle((1360,205,1367,415),fill='#c89562');d.rectangle((68,280,220,415),fill='#0b2030')
arrow(d,(1270,233),(1160,233),BLUE,3)
for xx0 in range(87,1350,20):
    for yy0 in range(294,400,20):
        if not (500<xx0<845 and yy0<350):d.ellipse((xx0-2,yy0-2,xx0+2,yy0+2),fill='#2f95cc')
d.polygon([(500,275),(845,275),(807,346),(545,346)],fill='#9a7044',outline='#ddc194');arrow(d,(670,245),(770,245),'#ffb783')
text(d,(542,161),'U(0) = 0.5 m/s',25,WHITE,True);text(d,(80,440),'Absorbing zone',20);text(d,(1160,440),'Wave piston',20)
text(d,(80,515),'Neighbor particles → pressure forces',26,WHITE,True);text(d,(80,569),'Pressure + gravity → acceleration',24,BLUE)
text(d,(865,515),'SURGE  SWAY  HEAVE',24,WHITE,True);text(d,(865,563),'ROLL    PITCH   YAW',24,WHITE,True)
text(d,(80,674),'Water forces + gravity move the rigid body',25)
for q in np.linspace(0,2*np.pi,9)[:-1]:
    x,y=1047+40*np.cos(q),349+40*np.sin(q);d.line([(1047,349),(x,y)],fill='#477c99',width=1)
d.ellipse((1003,305,1091,393),outline=BLUE,width=2);d.ellipse((1042,344,1052,354),fill=WHITE)
footer(d,'RTX 5080  |  0.94 M particles  |  40 mm spacing  |  Air is not simulated')

def frame(t):
    if t<3:return physical.copy()
    if t<6:return theory.copy()
    if t<14:
        tt=t-6;im,d=page('POTENTIAL: LINEAR WAVE RESPONSE','6-DOF harmonic response  ·  Reference speed 0.5 m/s',2)
        shot,m,z=potential(tt);im.paste(shot,(0,130));d=ImageDraw.Draw(im);profile(d,xs,z)
        text(d,(1100,465),f'Heave   {m[2]*100:+.1f} cm',25,WHITE);text(d,(1100,512),f'Pitch     {np.degrees(m[4]):+.1f}°',25,WHITE)
        text(d,(1100,579),'H = 24 cm',23,BLUE);text(d,(1100,617),'Wavelength 3 m',22)
        footer(d,'1:10 model  |  Large linear response: qualitative only  |  No coasting prediction')
    elif t<17:return sphsetup.copy()
    else:
        # Preserve the actual startup, wave arrival and late response in one scene.
        tt=min(st[-1],(t-17)/12*st[-1]);im,d=page('SPH: INITIAL SPEED, THEN FREE MOTION','Independent GPU simulation  ·  All six DOFs released',4)
        shot,center,angle,speed,z=sph(tt);im.paste(shot,(0,130));d=ImageDraw.Draw(im);profile(d,gridaxes[0],z)
        text(d,(1100,460),f'Surge U  {speed:+.2f} m/s',24,WHITE);text(d,(1100,507),f'Heave    {(center[2]-cg[2])*100:+.1f} cm',25,WHITE)
        text(d,(1100,554),f'Pitch      {angle[1]:+.1f}°',25,WHITE);text(d,(1100,601),f'Roll        {angle[0]:+.1f}°',25,WHITE)
        text(d,(1100,670),f't = {tt:.1f} s',23,BLUE)
        text(d,(255,710),'18 m TANK',15,GRAY);d.line([(385,723),(1025,723)],fill='#54758d',width=2)
        pos=385+(center[0]+9)/18*640;half=2.606/18*640/2
        d.rounded_rectangle((pos-half,716,pos+half,730),radius=4,fill='#b89461');d.line([(1025,712),(1025,733)],fill=BLUE,width=3)
        footer(d,'1:10 model  |  Actual SPH particles / body motion  |  Surface smoothed for display','RTX 5080')
    text(d,(64,168),'AIR',21,GRAY,True);text(d,(92,652),'WATER',22,BLUE,True);text(d,(64,704),'LOCAL VIEW',15,GRAY)
    return im

if not args.library_only:
    physical.save(out/'turtle_history_physical_model.png');theory.save(out/'turtle_history_potential_setup.png');sphsetup.save(out/'turtle_history_sph_setup.png')
    for t,name in [(8.2,'turtle_history_potential_preview.png'),(27.5,'turtle_history_sph_preview.png')]:frame(t).save(out/name)
    if args.stills_only:pl.close();raise SystemExit()
    video=out/'turtle_ship_story_v3_6dof_16x9.mp4'
    cmd=[shutil.which('ffmpeg'),'-y','-f','rawvideo','-pix_fmt','rgb24','-s','1440x810','-r','30','-i','-','-an','-c:v','libx264','-preset','fast','-crf','19','-pix_fmt','yuv420p','-movflags','+faststart',str(video)]
    with (out/'turtle_v3_render.log').open('w') as log:
        proc=subprocess.Popen(cmd,stdin=subprocess.PIPE,stdout=log,stderr=log)
        for i in range(870):
            proc.stdin.write(np.asarray(frame(i/30)).tobytes())
            if i%120==0:print('frame',i,'/ 870',flush=True)
        proc.stdin.close();assert proc.wait()==0
    pl.close()
    meta=dict(duration_s=29,size=[1440,810],sections=[['physical model',0,3],['potential principle/setup',3,6],['potential result',6,14],['SPH principle/setup',14,17],['SPH result',17,29]],
        potential='6-DOF harmonic response at a fixed reference speed; not a coasting transient.',sph='Initial velocity only, all six DOFs free; independent solver.',sph_last_saved_time_s=float(st[-1]),
        presentation='Interpolated saved poses and smoothed particle occupancy surfaces. No potential waves are supplied to SPH.',sha256=hashlib.sha256(video.read_bytes()).hexdigest())
    (out/(video.stem+'.json')).write_text(json.dumps(meta,indent=2));print(video)
