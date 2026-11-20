"""Render native coupled body transforms and reconstructed SPH liquid surfaces.

Unlike render_intro_motion.py this script contains no keyframed body or water
motion. The camera observes the prescribed accelerating support frame.
"""
from pathlib import Path
import argparse,json,subprocess
import numpy as np
import pyvista as pv
from PIL import Image,ImageDraw
from scipy.integrate import cumulative_trapezoid
from analyze_coupled import convert
from render_braking_results import surface
from render_intro_motion import Scene
from render_minimal_intro import W,H,FPS,BG,WHITE,BLUE,ORANGE,text

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'media'
CG0=np.array([.8,0,-.9])
LABELS={'A_unpartitioned':'NO BAFFLES','C_separate_compartments':'COMPARTMENTS',
        'B_perforated_baffles':'PERFORATED','D_cross_baffles':'CROSS BAFFLES'}


class CoupledScene(Scene):
    def __init__(self,size,case,thickness=.08):
        super().__init__(size,truck=False,case=case)
        # First six actors are the shell and its outline. Match plate thickness
        # to the actual voxel envelope, including the stronger four-layer case.
        for j,actor in enumerate(self.actors[6:-1]):
            mesh=actor.mapper.dataset.copy(deep=True)
            if case=='D_cross_baffles' and j==2:mesh.points[:,1]*=thickness/.06
            else:
                cx=(-1.,1.)[j];mesh.points[:,0]=cx+(mesh.points[:,0]-cx)*thickness/.06
            actor.mapper.SetInputData(mesh)
        self.water.GetProperty().SetOpacity(.88)
        self.add(pv.Box(bounds=(-3.6,4.7,-.8,.8,-1.35,-1.12)),'#718491')
        self.add(pv.Box(bounds=(3.25,4.65,-.86,.86,-1.12,.65)),'#a2b3c0')
        self.add(pv.Box(bounds=(3.45,4.58,-.885,-.87,-.1,.48)),'#193a50')
        for x in (-2.1,-.8,3.8):
            for y in (-1.,1.):
                self.add(pv.Cylinder(center=(x,y,-1.6),direction=(0,1,0),radius=.5,height=.25,resolution=32),'#263746')
                self.add(pv.Cylinder(center=(x,y+np.sign(y)*.13,-1.6),direction=(0,1,0),radius=.25,height=.02,resolution=32),'#8ba1b1')
        self.brake_lights=[]
        for y in (-.65,.65):
            self.brake_lights.append(self.add(pv.Box(bounds=(-3.65,-3.60,y-.13,y+.13,-1.25,-1.03)),'#63342f',lighting=False))
        # Reference support surface, fixed in the acceleration frame.
        self.p.add_mesh(pv.Plane(center=(0,0,-2.104),direction=(0,0,1),i_size=64,j_size=64),color='#0b1420',ambient=1)
        for y in (-2.8,2.8):
            self.p.add_mesh(pv.Line((-15,y,-2.099),(15,y,-2.099)),color='#2b4052',line_width=2,lighting=False)
        self.road_cues=[]
        dashes=[]
        for x in range(-60,61,6):
            for y in (-2.5,2.5):
                dashes.append(pv.Line((x,y,-2.098),(x+2.0,y,-2.098)))
        self.road_cues.append(self.p.add_mesh(pv.merge(dashes),color='#8d9eaa',line_width=4,lighting=False))
        # Sparse roadside posts are visual reference points. They scroll by the
        # prescribed forward speed, not an independently solved road trajectory.
        posts=[];caps=[]
        for x in range(-60,61,6):
            for y in (-4.,4.):
                posts.append(pv.Box(bounds=(x-.06,x+.06,y-.06,y+.06,-2.1,-1.48)))
                caps.append(pv.Box(bounds=(x-.065,x+.065,y-.065,y+.065,-1.65,-1.54)))
        self.road_cues.append(self.p.add_mesh(pv.merge(posts),color='#8295a4',ambient=.8))
        self.road_cues.append(self.p.add_mesh(pv.merge(caps),color='#b2bfca',lighting=False))
        self.p.camera_position=[(12,-24,11),(.5,-.8,-.7),(0,0,1)]
        self.p.camera.parallel_projection=True
        self.p.camera.parallel_scale=5.8

    def actual_frame(self,rotation,center,water,distance=0.,braking=False):
        transform=np.eye(4);transform[:3,:3]=rotation;transform[:3,3]=center-rotation@CG0
        for actor in self.actors:actor.user_matrix=transform
        self.water.mapper.SetInputData(water)
        shift=np.array([center[0]-CG0[0],center[1]-CG0[1],0.])
        self.p.camera.position=np.array([12.,-24.,11.])+shift
        self.p.camera.focal_point=np.array([.5,-.8,-.7])+shift
        for cue in self.road_cues:cue.position=(-float(distance)%6.,0.,0.)
        for lamp in self.brake_lights:
            lamp.GetProperty().SetColor(*(np.array([242,75,50]) if braking else np.array([99,52,47]))/255)
        self.p.render()
        return Image.fromarray(self.p.screenshot(return_img=True))


def prepare(folder):
    pose=np.load(folder/'body_pose.npz')
    dest,rr=convert(folder,folder/'output',label='Fluid',fluid=True)
    cache=folder/'coupled_surfaces';cache.mkdir(exist_ok=True)
    for i in range(len(rr)):
        path=cache/f'water_{i:04d}.vtp'
        if path.exists():continue
        points=pv.read(dest/f'Fluid_{i:04d}.vtk').points.astype(float)
        local=(points-pose['center'][i])@pose['rotation'][i]+CG0
        surface(local).save(path)
    return pose,cache


def render(names,preview=False):
    folders=[ROOT/'sph'/n for n in names]
    analysis=[json.loads((f/'coupled_analysis.json').read_text()) for f in folders]
    setup=[json.loads((f/'setup.json').read_text()) for f in folders]
    acceleration=np.loadtxt(folders[0]/'input_acceleration.csv',delimiter=',',skiprows=1)
    drive_time=acceleration[:,0]
    drive_speed=np.maximum(0.,setup[0]['brake']['initial_speed_m_s']+cumulative_trapezoid(acceleration[:,1],drive_time,initial=0.))
    drive_distance=cumulative_trapezoid(drive_speed,drive_time,initial=0.)
    for s in setup[1:]:
        assert s['lateral_g']==setup[0]['lateral_g'] and s['fluid_mass_kg']==setup[0]['fluid_mass_kg']
    records=[prepare(f) for f in folders]
    layouts=[n.split('_coupled_')[0] for n in names]
    grid=len(names)>1
    scenes=[CoupledScene((510,500) if grid else (1080,1100),layout,s.get('numerical_plate_thickness_m',.08)) for layout,s in zip(layouts,setup)]
    count=min(len(p['time_s']) for p,_ in records)
    chosen=list(range(0,count,2))
    stem='tanker_coupled_comparison_roadmotion_portrait' if grid else 'tanker_'+names[0]+'_roadmotion_portrait'
    def frame(i,display_time=None):
        im=Image.new('RGB',(W,H),BG);d=ImageDraw.Draw(im)
        text(d,(65,150),'SAME INPUT.' if grid else 'BRAKE + SWERVE',65,bold=True)
        text(d,(65,235),'DIFFERENT RESPONSE.' if grid else LABELS[layouts[0]],52,bold=True)
        text(d,(65,345),f'PEAK: BRAKE 0.40g + LATERAL {setup[0]["lateral_g"]:.2f}g',27,fill='#a3bbcc')
        t=records[0][0]['time_s'][i] if display_time is None else display_time
        text(d,(65,410),f'{t:.1f} s',32,fill=BLUE)
        speed=float(np.interp(t,drive_time,drive_speed));distance=float(np.interp(t,drive_time,drive_distance))
        text(d,(720,390),'INPUT SPEED',17,fill='#748ea3')
        text(d,(720,417),f'{speed*3.6:.0f} km/h',32,fill=BLUE)
        for k,(scene,(p,cache),layout) in enumerate(zip(scenes,records,layouts)):
            shot=scene.actual_frame(p['rotation'][i],p['center'][i],pv.read(cache/f'water_{i:04d}.vtp'),distance,t>=setup[0]['brake']['brake_start_s'])
            x,y=(30+510*(k%2),560+570*(k//2)) if grid else (0,530)
            im.paste(shot,(x,y));d=ImageDraw.Draw(im)
            if grid:text(d,(x+20,y+15),LABELS[layout],26,fill=WHITE if k==0 else ORANGE,bold=True)
            roll=float(p['euler_deg'][i,0])
            roll_time=analysis[k].get('first_75deg_time_s')
            has_rolled=analysis[k]['overturned'] and roll_time is not None and p['time_s'][i]>=roll_time+.2
            flag='ROLLED' if has_rolled else 'ROLLING' if p['tilt_deg'][i]>=75 else 'WHEEL LIFT' if p['side_lift'][i] else ''
            text(d,(x+20,y+(470 if grid else 1020)),f'{abs(roll):.0f} deg  {flag}',26 if grid else 38,fill='#ec9b68' if flag else '#a3bbcc')
        text(d,(65,1785),'SPH + RIGID BODY  /  FIRST RESOLUTION',23,fill='#748ea3')
        return im
    if preview:
        for t in (0,3,4,5,7.9):
            i=min(chosen,key=lambda j:abs(records[0][0]['time_s'][j]-t))
            frame(i).save(OUT/f'{stem}_t{t:.1f}.png')
    else:
        path=OUT/(stem+'.mp4')
        proc=subprocess.Popen(['ffmpeg','-v','error','-y','-f','rawvideo','-pix_fmt','rgb24','-s',f'{W}x{H}',
                               '-r',str(FPS),'-i','-','-an','-c:v','libx264','-preset','fast','-crf','19','-pix_fmt','yuv420p',
                               '-movflags','+faststart',str(path)],stdin=subprocess.PIPE)
        frame_count=round(float(records[0][0]['time_s'][count-1])*FPS)+1
        for k in range(frame_count):
            t=k/FPS
            # Hold the nearest saved 20 Hz physical field. Only the prescribed
            # road cue runs continuously at 30 Hz; no water motion is invented.
            i=min(count-1,int(np.argmin(abs(records[0][0]['time_s']-t))))
            proc.stdin.write(frame(i,t).tobytes())
            if k%60==0:print(f'render {k}/{frame_count}',flush=True)
        proc.stdin.close()
        if proc.wait():raise RuntimeError('Encoder failed')
        meta=dict(video=path.name,cases=names,frames=frame_count,fps=FPS,seconds=frame_count/FPS,
                  body_pose='Rigid fit to native boundary-particle positions. No authored rotation or translation.',
                  water='Actual particle occupancy surface reconstructed in instantaneous body coordinates.',
                  temporal_render='Nearest saved 20 Hz fluid/body field; prescribed road-motion cue at 30 Hz.',
                  reference_frame='Body motion in nonrotating prescribed accelerating support frame; not a solved steering trajectory.',
                  environment='Road dashes and roadside posts scroll by the prescribed initial 50 km/h and integrated braking history. Displayed speed is this prescribed reference speed; it does not include solved relative sliding. Brake lights are a visual cue.',
                  vehicle='Hypothetical rigid vehicle; fixed wheels, frictional road contact, no suspension or calibrated tires.',
                  source_summary=analysis)
        (OUT/(stem+'.json')).write_text(json.dumps(meta,indent=2))
        print(json.dumps({'video':str(path),'seconds':meta['seconds']}),flush=True)
    for scene in scenes:scene.close()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('cases',nargs='+');p.add_argument('--preview',action='store_true')
    a=p.parse_args();render(a.cases,a.preview)
