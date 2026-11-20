"""One-vessel film: setup / SPH schematic / close-up and diagnostic views."""
from pathlib import Path
import argparse,json,subprocess,time
import numpy as np
import pyvista as pv
from PIL import Image,ImageDraw
from scipy.spatial.transform import Rotation
from matplotlib import colormaps
from render_gunnerus_comparison import Scene,prepare,text,BG,W,H,FPS,ROOT,OUT
from prepare_gunnerus_surface_pressure import prepare_pressure

NAME='gunnerus_focus_A085_x2_050_16s'
FOLDER=ROOT/'sph'/NAME
PREFIX='gunnerus_focus'


def arrow(draw,start,end,color,width=7):
    a=np.array(start,float); b=np.array(end,float); d=b-a; d/=np.linalg.norm(d)
    n=np.array([-d[1],d[0]]); draw.line([tuple(a),tuple(b)],fill=color,width=width)
    draw.polygon([tuple(b),tuple(b-d*22+n*11),tuple(b-d*22-n*11)],fill=color)


def intro_frames():
    setup=json.loads((FOLDER/'setup.json').read_text()); p=setup['physical_model']
    initial=pv.read(FOLDER/'generated/Case_All.vtk'); xyz=initial.points; station=p['initial_station_x_m']
    # The section uses the actual initial lattice, not decorative water points.
    near=(abs(xyz[:,1])<.026)&(abs(xyz[:,0]-station)<2.3)&(xyz[:,2]>-.65)&(xyz[:,2]<1.1)
    method=Image.new('RGB',(W,H),BG); d=ImageDraw.Draw(method)
    text(d,(60,60),'SPH',72,bold=True); text(d,(60,165),'PARTICLES TO SHIP MOTION',32,'#a8c4d4')
    def point(v):return (int(540+(v[0]-station)*210),int(780-v[2]*210))
    for typ,color,radius in [(3,'#318ec2',3),(2,'#e7b769',3)]:
        for v in xyz[near&(initial['Type']==typ)]:
            x,y=point(v); d.ellipse((x-radius,y-radius,x+radius,y+radius),fill=color)
    text(d,(60,365),'WATER PARTICLES',28,'#64cbe7',True)
    text(d,(60,1010),'SOLID BOUNDARY',28,'#e7b769',True)
    text(d,(60,1060),'Hull + deckhouse',26,'#9ab0c0')
    # Circular inset illustrates the neighbor interaction used by SPH.
    cx,cy=210,1260; d.ellipse((cx-95,cy-95,cx+95,cy+95),outline='#38586d',width=3)
    for x,y in [(0,0),(-50,-30),(45,-45),(-35,60),(55,38),(5,72)]:
        if x or y:d.line((cx,cy,cx+x,cy+y),fill='#355c77',width=2)
        d.ellipse((cx+x-7,cy+y-7,cx+x+7,cy+y+7),fill='#64cbe7')
    text(d,(360,1210),'LOCAL FORCES',34,bold=True); text(d,(360,1270),'Pressure interactions',27,'#9ab0c0')
    text(d,(60,1465),'RIGID BODY  /  FREE 6 DOF',34,bold=True)
    text(d,(60,1530),'Fluid loads + weight',27,'#9ab0c0')
    text(d,(60,1770),'DualSPHysics  /  GPU SPH  /  mDBC hull',23,'#738d9e')
    method.save(OUT/f'{PREFIX}_method.png')

    pl=pv.Plotter(off_screen=True,window_size=(1080,700));pl.set_background(BG)
    pl.add_mesh(pv.Box(bounds=setup['domain_bounds_m']).outline(),color='#69869a',line_width=2)
    pl.add_mesh(pv.Plane(center=(0,0,0),direction=(0,0,1),i_size=36,j_size=10),color='#0b3855',opacity=.65)
    pl.add_mesh(pv.read(p['hull_file']),color='#e3b75e',smooth_shading=True)
    pl.add_mesh(pv.Box(bounds=(17.98,18.03,-5,5,-setup['water_depth_m'],3)),color='#64cbe7',opacity=.45)
    pl.add_mesh(pv.Box(bounds=(-18,-11,-5,5,-setup['water_depth_m'],.08)),color='#486071',opacity=.3)
    pl.add_mesh(pv.Arrow(start=(13,0,1),direction=(-1,0,0),scale=6),color='#64cbe7')
    pl.camera_position=[(25,-37,27),(0,0,-.2),(0,0,1)];pl.camera.parallel_projection=True;pl.camera.parallel_scale=13.5
    pl.enable_anti_aliasing('ssaa'); shot=Image.fromarray(pl.screenshot(return_img=True));pl.close()
    layout=Image.new('RGB',(W,H),BG);layout.paste(shot,(0,330));d=ImageDraw.Draw(layout)
    text(d,(60,60),'GUNNERUS',60,bold=True);text(d,(60,158),'FOCUSED HEAD WAVES',34,'#64cbe7',True)
    text(d,(60,246),'LARGER WAVE BASIN',29,'#a8c4d4')
    text(d,(60,1080),f'360 x 100 m  /  {setup["water_depth_m"]*10:g} m depth',39,bold=True)
    text(d,(60,1145),'Full-scale equivalent  /  1:10 simulation',26,'#9ab0c0')
    text(d,(60,1260),'MEASURED INCIDENT WAVE',28,'#64cbe7',True)
    hlow,hhigh=setup.get('incident_wave_height_equivalent_m',[8,11])
    text(d,(60,1320),f'H ~ {hlow:g}-{hhigh:g} m equivalent',39,bold=True)
    text(d,(60,1390),'Empty-tank crest-to-trough measurement',24,'#9ab0c0')
    text(d,(60,1525),'Initial speed 8 kn  /  Assumed thrust',27,'#a8c4d4')
    text(d,(60,1580),'Free 6 DOF  /  Waves drive ship motion',26,'#a8c4d4')
    particles=json.loads((FOLDER/'initial_audit.json').read_text())['total_particles']
    text(d,(60,1760),f'{particles/1e6:.2f} M particles  /  RTX 5080',24,'#738d9e')
    layout.save(OUT/f'{PREFIX}_setup.png')
    return method,layout


def add_bow_monitor():
    p=json.loads((FOLDER/'setup.json').read_text())['physical_model']
    with np.load(FOLDER/'monitors.npz') as q:data={k:q[k] for k in q.files}
    hull=pv.read(p['hull_file']); bx=p['bounds_m'][1][0]-.18
    point,_=hull.ray_trace((bx,0,3),(bx,0,-2),first_point=True);assert len(point)==3
    r=Rotation.from_euler('xyz',data['euler_xyz_deg'],degrees=True).apply(np.asarray(point)-p['center_of_mass_m'])
    acc=data['facelin']+np.cross(data['faceang'],r)+np.cross(data['fomega'],np.cross(data['fomega'],r))
    data['bow_accel_z_g']=acc[:,2]/9.81
    np.savez_compressed(FOLDER/'monitors.npz',**data)
    (FOLDER/'bow_monitor.json').write_text(json.dumps(dict(point_initial_m=point.tolist(),
        definition='World-z kinematic acceleration of a fixed foredeck point, divided by g; not accelerometer specific force.',
        formula='a_point = a_CG + alpha cross r + omega cross (omega cross r)',
        source='Native rigid-body facelin, faceang and fomega; saved particle-pose rotation.'),indent=2))


class Result:
    def __init__(self,show_dilute_water=False,surface_level=.5):
        self.scene=Scene(NAME,window_size=(1080,720),parallel_scale=1.9,water_x=6.8,water_y=4.9,water_bottom=-1.1,water_top=2.8,
            show_dilute_water=show_dilute_water,surface_level=surface_level)
        with np.load(FOLDER/'monitors.npz') as q:self.data={k:q[k] for k in q.files}
        self.arrow=self.scene.pl.add_mesh(pv.Arrow(start=(0,0,-99),direction=(0,0,1)),color='#64cbe7',reset_camera=False)
        pressure_path=prepare_pressure(NAME)
        with np.load(pressure_path) as q:
            self.pressure_times=q['time_s']; self.pressures=q['pressure_kPa']
        self.skin=pv.read(pressure_path.parent/'contact_mesh.vtp')
        self.skin_reference=self.skin.points.copy()
        values=self.pressures[np.isfinite(self.pressures)]
        self.pressure_limits=(min(0.,float(np.floor(np.percentile(values,.5)))),max(1.,float(np.ceil(np.percentile(values,99.5)*2)/2)))
        self.skin.cell_data['pressure_kPa']=self.pressures[0]
        self.skin_actor=self.scene.pl.add_mesh(self.skin,scalars='pressure_kPa',preference='cell',cmap='turbo',
            clim=self.pressure_limits,nan_color='#64717d',show_scalar_bar=False,lighting=False,reset_camera=False)
        self.skin_actor.SetVisibility(False)
        self.pressure_metadata=dict(source=str(pressure_path.with_suffix('.json')),color_limits_kPa=list(self.pressure_limits),
            color_range='Fixed throughout the clip; 0.5 and 99.5 percentiles rounded outward. Values beyond limits saturate.',
            temporal_interpolation='Linear only when both samples have water support; otherwise nearest snapshot, gray if unsupported.')
        self.cards=[('NET VERTICAL','N','net_vertical_N',None,1),('PITCH MOMENT','N m','moment_body_Nm',1,1),
                    ('PITCH','deg','euler_xyz_deg',1,1),('BOW ACCEL','g','bow_accel_z_g',None,1)]

    def value(self,key,t,col=None):
        q=self.data[key] if col is None else self.data[key][:,col]
        return float(np.interp(t,self.data['time_s'],q))

    def chart(self,draw,xy,card,t):
        x,y=xy;title,unit,key,col,factor=card
        v=(self.data[key] if col is None else self.data[key][:,col])*factor;tt=self.data['time_s']
        use=(tt>=4)&(tt<=16.001);tt=tt[use];v=v[use];lo=float(v.min());hi=float(v.max());pad=max((hi-lo)*.1,.01)
        lo-=pad;hi+=pad; py=lambda q:y+150-(q-lo)/(hi-lo)*80
        text(draw,(x,y),title,24,'#b9d0de',True);text(draw,(x+365,y),unit,23,'#829dad')
        text(draw,(x,y+32),f'{self.value(key,t,col)*factor:+.2f}',32,'#64cbe7',True)
        if lo<=0<=hi:draw.line((x,py(0),x+450,py(0)),fill='#2b4051',width=1)
        draw.line((x,y+150,x+450,y+150),fill='#2b4051',width=1)
        keep=tt<=t;points=list(zip((x+(tt[keep]-4)/12*450).tolist(),py(v[keep]).tolist()))
        if len(points)>1:draw.line(points,fill='#64cbe7',width=3)
        cursor=x+(t-4)/12*450;draw.line((cursor,y+76,cursor,y+150),fill='#8fa5b4',width=1)
        text(draw,(x,y+154),'4 s',17,'#6c8494');text(draw,(x+418,y+154),'16 s',17,'#6c8494')

    def auxiliary_views(self,t):
        j=int(np.clip(np.searchsorted(self.pressure_times,t)-1,0,len(self.pressure_times)-2))
        f=float(np.clip((t-self.pressure_times[j])/(self.pressure_times[j+1]-self.pressure_times[j]),0,1))
        a,b=self.pressures[j],self.pressures[j+1]
        p=(a if f<.5 else b).copy();wet=np.isfinite(a)&np.isfinite(b);p[wet]=(1-f)*a[wet]+f*b[wet]
        pose=self.scene.pose;cg=self.scene.cg
        self.skin.points=self.skin_reference@pose[:3,:3].T+pose[:3,3]
        self.skin.cell_data['pressure_kPa']=p
        self.skin_actor.SetVisibility(True);self.arrow.SetVisibility(False)
        parcel_visible=self.scene.parcels.GetVisibility() if self.scene.parcels is not None else False
        if self.scene.parcels is not None:self.scene.parcels.SetVisibility(False)
        for actor,_ in self.scene.actors:actor.GetProperty().SetOpacity(.15)
        self.scene.water.GetProperty().SetOpacity(.16)
        side=self.scene.capture([(cg[0],-20,.55),(cg[0],0,.55),(0,0,1)],1.18,(490,230))
        front=self.scene.capture([(cg[0]+20,0,.55),(cg[0],0,.55),(0,0,1)],1.02,(490,230))
        self.skin_actor.SetVisibility(False);self.arrow.SetVisibility(True)
        if self.scene.parcels is not None:self.scene.parcels.SetVisibility(parcel_visible)
        for actor,_ in self.scene.actors:actor.GetProperty().SetOpacity(1.)
        self.scene.water.GetProperty().SetOpacity(1.)
        self.scene.pl.window_size=(1080,720)
        return side,front

    def frame(self,t):
        f=self.value('net_vertical_N',t);cg=np.array([self.value('center',t,i) for i in range(3)])
        length=float(np.clip(np.sqrt(abs(f)/1000),.2,1.3));sign=1 if f>=0 else -1
        self.arrow.mapper.SetInputData(pv.Arrow(start=(cg[0],-1.0,.7),direction=(0,0,sign),scale=length))
        self.arrow.GetProperty().SetColor(*(np.array([100,203,231]) if f>=0 else np.array([239,162,94]))/255)
        image=Image.new('RGB',(W,H),BG);image.paste(self.scene.render(t),(0,145))
        side,front=self.auxiliary_views(t);image.paste(side,(40,915));image.paste(front,(550,915));d=ImageDraw.Draw(image)
        text(d,(60,42),'GUNNERUS',48,bold=True);text(d,(770,56),'1:10 MODEL',26,'#9ab0c0')
        text(d,(60,112),'FOCUSED HEAD WAVES',29,'#64cbe7',True);text(d,(800,112),f't {t:05.2f} s',25)
        text(d,(60,858),f'U {self.value("speed_body_m_s",t):+.2f} m/s',27)
        text(d,(590,858),f'Heave {self.value("heave_m",t)*1000:+.0f} mm',27)
        text(d,(60,909),'SIDE',23,'#b9d0de',True);text(d,(570,909),'FRONT',23,'#b9d0de',True)
        text(d,(60,1160),'FLUID-SIDE PRESSURE',22,'#b9d0de',True)
        gradient=(colormaps['turbo'](np.linspace(0,1,230))[:,:3]*255).astype('uint8')
        image.paste(Image.fromarray(np.repeat(gradient[None,:,:],13,axis=0)),(615,1169))
        low,high=self.pressure_limits
        text(d,(535,1158),f'≤{low:g}',22,'#b9d0de');text(d,(860,1158),f'{high:g}+ kPa',22,'#b9d0de')
        for xy,card in zip([(60,1220),(570,1220),(60,1410),(570,1410)],self.cards):self.chart(d,xy,card,t)
        text(d,(60,1630),'Arrow: net vertical force  /  Model units',23,'#91a9ba')
        text(d,(60,1670),'Pressure: nearby water  /  Gray: no sample',22,'#7b93a4')
        text(d,(60,1710),'Pitch: bow-down positive  /  Assumed thrust',22,'#7b93a4')
        text(d,(60,1780),'Visual model: IHB NTNU / vessel.js | MIT',21,'#627f92')
        return image


def encode(path,frames,count):
    proc=subprocess.Popen(['ffmpeg','-v','error','-y','-f','rawvideo','-pix_fmt','rgb24','-s',f'{W}x{H}','-r',str(FPS),'-i','-',
        '-an','-c:v','libx264','-preset','fast','-crf','19','-pix_fmt','yuv420p','-movflags','+faststart',str(path)],stdin=subprocess.PIPE)
    try:
        for i,frame in enumerate(frames):
            proc.stdin.write(np.asarray(frame).tobytes())
            if i%60==0:print(f'Encoded {i}/{count}',flush=True)
    finally:proc.stdin.close()
    assert proc.wait()==0


def main(intro_only,video_name='gunnerus_focused_waves_story_portrait',show_dilute_water=False,surface_level=.5):
    OUT.mkdir(exist_ok=True); a,b=intro_frames()
    def intro():
        for frame in [b,a]:
            for _ in range(60):yield frame
    if intro_only:
        encode(OUT/f'{PREFIX}_intro_preview.mp4',intro(),120);return
    add_bow_monitor();prepare(NAME,clip_x=7.3,clip_y=5.01,z_low=-1.3,z_high=3.)
    result=Result(show_dilute_water=show_dilute_water,surface_level=surface_level);result.frame(8.5).save(OUT/f'{PREFIX}_result_preview.png')
    def full():
        yield from intro()
        for t in np.linspace(4.,16.,18*FPS):yield result.frame(float(t))
    path=OUT/f'{video_name}.mp4';encode(path,full(),660)
    path.with_suffix('.json').write_text(json.dumps(dict(case=NAME,frames=660,fps=FPS,video_seconds=22,dimensions=[W,H],
        segments=[dict(kind='domain setup',video_s=[0,2]),dict(kind='SPH schematic',video_s=[2,4]),dict(kind='actual simulation replay',video_s=[4,22],model_time_s=[4,16])],
        metrics=['net vertical force','body pitch moment','pitch','foredeck vertical kinematic acceleration','speed','heave'],
        layout='Large close perspective view; small world side/front orthographic views with fluid-side pressure; synchronized native load/motion histories.',
        pressure=result.pressure_metadata,
        dilute_water_parcels=show_dilute_water,surface_level=surface_level,
        rendering='Matched actual particle IDs with temporal interpolation; rigid pose SLERP; no invented wave or extra downward force.',
        limitations=['Assumed CG thrust, not actuator disk.','Sealed simplified deckhouse; no flooding, air entrapment or structural failure model.']),indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--intro-only',action='store_true');ap.add_argument('--case',default=NAME)
    ap.add_argument('--prefix',default=PREFIX);ap.add_argument('--video-name',default='gunnerus_focused_waves_story_portrait')
    ap.add_argument('--dilute-water',action='store_true');ap.add_argument('--surface-level',type=float,default=.5)
    a=ap.parse_args();NAME=a.case;FOLDER=ROOT/'sph'/NAME;PREFIX=a.prefix
    main(a.intro_only,a.video_name,a.dilute_water,a.surface_level)
