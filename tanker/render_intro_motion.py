"""Authored risk illustration followed by four parametric study geometries.

Not a dynamics solution: vehicle motion, roll and planar liquid surface are
keyframed solely for the explicitly labelled concept intro.
"""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess

import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw

from render_minimal_intro import W, H, FPS, BG, WHITE, BLUE, ORANGE, text

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'media'
GEO = ROOT/'geometry'/'concept_v1'
GROUND = -2.04
CASES = [
    ('A_unpartitioned','NO BAFFLES'),
    ('C_separate_compartments','COMPARTMENTS'),
    ('B_perforated_baffles','PERFORATED BAFFLES'),
    ('D_cross_baffles','CROSS BAFFLES'),
]


def smooth(a,b,t):
    q=np.clip((t-a)/(b-a),0.,1.)
    return q*q*(3-2*q)


def rotation(roll,yaw=0):
    c,s=np.cos(roll),np.sin(roll)
    rx=np.array([[1,0,0],[0,c,-s],[0,s,c]])
    c,s=np.cos(yaw),np.sin(yaw)
    return np.array([[c,-s,0],[s,c,0],[0,0,1]])@rx


def cylinder():
    return pv.Cylinder(center=(0,0,0),direction=(1,0,0),radius=1,height=6,resolution=128).triangulate().clean()


def plate(x,holes):
    if not holes:
        return pv.Cylinder(center=(x,0,0),direction=(1,0,0),radius=1.005,height=.06,resolution=128).triangulate().clean()
    n=105
    q=np.linspace(-1.01,1.01,n)
    yy,zz=np.meshgrid(q,q,indexing='ij')
    points=np.column_stack((np.full(n*n,x-.03),yy.ravel(),zz.ravel()))
    faces=[]
    for i in range(n-1):
        for j in range(n-1):
            y,z=(q[i]+q[i+1])/2,(q[j]+q[j+1])/2
            if y*y+z*z>1.005**2:
                continue
            if any((y-cy)**2+(z-cz)**2<.24**2 for cy in (-.42,.42) for cz in (-.40,.40)):
                continue
            k=i*n+j
            faces.extend([4,k,k+n,k+n+1,k+1])
    return pv.PolyData(points,np.asarray(faces)).clean().triangulate().extrude((.06,0,0),capping=True).clean()


def parts(case):
    if case=='A_unpartitioned':
        return []
    result=[plate(x,case!='C_separate_compartments') for x in (-1,1)]
    if case=='D_cross_baffles':
        result.append(pv.Box(bounds=(-3.0,3.0,-.03,.03,-1.005,1.005)).triangulate())
    return result


def export_geometry():
    GEO.mkdir(parents=True,exist_ok=True)
    cylinder().save(GEO/'tank_envelope.stl')
    records=[]
    for cid,label in CASES:
        plates=parts(cid)
        paths=[]
        for i,mesh in enumerate(plates):
            path=GEO/f'{cid}_plate_{i+1}.stl'
            mesh.save(path)
            paths.append({'file':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'open_edges':mesh.n_open_edges})
        records.append({'id':cid,'label':label,'plates':paths})
    data={'status':'Parametric proposed study geometry, not yet particle-qualified or a runnable solver setup',
          'units':'m','tank_length':6.,'tank_radius':1.,'plate_positions_x':[-1.,1.],
          'plate_thickness':.06,'hole_radius':.24,'hole_centres_yz':[[-.42,-.4],[-.42,.4],[.42,-.4],[.42,.4]],
          'note':'Plate thickness and faceting are provisional numerical geometry. No exact manufacturer CAD. Compartments require separate vented gas spaces. Cross-plate intersections require particle-generation checks.',
          'cases':records}
    (GEO/'geometry.json').write_text(json.dumps(data,indent=2),encoding='utf-8')
    return data


class Scene:
    def __init__(self,size,truck=False,case='A_unpartitioned'):
        self.p=pv.Plotter(off_screen=True,window_size=size)
        self.p.set_background(BG)
        self.p.enable_anti_aliasing('ssaa')
        self.p.enable_depth_peeling(number_of_peels=8,occlusion_ratio=0)
        self.actors=[]
        self.points=[]
        self.wheels=[]
        self.truck=truck
        self.barrel=cylinder()
        self.add(self.barrel,'#aec0cd',opacity=.055)
        angle=np.linspace(0,2*np.pi,140)
        for x in (-3,3):
            line=pv.lines_from_points(np.column_stack((np.full_like(angle,x),np.cos(angle),np.sin(angle))),close=True)
            self.add(line,'#96adbe',line_width=3,lighting=False)
        for th in (np.pi/2,3*np.pi/2,np.pi):
            self.add(pv.Line((-3,np.cos(th),np.sin(th)),(3,np.cos(th),np.sin(th))),'#587186',line_width=2,lighting=False)
        for mesh in parts(case):
            self.add(mesh,ORANGE,ambient=.7)
        self.water=self.add(self.water_mesh((0,0,1)), '#2784b1',opacity=.58,ambient=.6,specular=.15)
        if truck:
            self.add(pv.Box(bounds=(-3.6,4.7,-.8,.8,-1.35,-1.12)),'#718491')
            self.add(pv.Box(bounds=(3.25,4.65,-.86,.86,-1.12,.65)),'#a2b3c0')
            self.add(pv.Box(bounds=(3.45,4.58,-.885,-.87,-.1,.48)),'#193a50')
            self.add(pv.Box(bounds=(4.66,4.67,-.72,.72,-.1,.48)),'#193a50')
            for x in (-2.1,-.8,3.8):
                for y in (-.92,.92):
                    self.add(pv.Cylinder(center=(x,y,-1.55),direction=(0,1,0),radius=.48,height=.22,resolution=40),'#263746')
                    self.add(pv.Cylinder(center=(x,y+np.sign(y)*.12,-1.55),direction=(0,1,0),radius=.25,height=.03,resolution=32),'#8ba1b1')
            self.lights=[]
            for y in (-.65,.65):
                self.lights.append(self.add(pv.Box(bounds=(-3.65,-3.60,y-.13,y+.13,-1.25,-1.03)),'#63342f'))
            self.p.add_mesh(pv.Plane(center=(0,0,GROUND-.015),direction=(0,0,1),i_size=36,j_size=7),color='#0b1420',ambient=1)
            for y in (-2.1,2.1):
                self.p.add_mesh(pv.Line((-18,y,GROUND),(18,y,GROUND)),color='#2b4052',line_width=2,lighting=False)
            for x in range(-18,19,3):
                self.p.add_mesh(pv.Line((x,-1.8,GROUND+.005),(x+1.25,-1.8,GROUND+.005)),color='#4a6072',line_width=3,lighting=False)
            self.p.camera_position=[(11,-23,10),(.5,-.8,-.4),(0,0,1)]
            self.p.camera.parallel_projection=True
            self.p.camera.parallel_scale=5.1
        else:
            self.p.camera_position=[(9,-20,8),(0,0,0),(0,0,1)]
            self.p.camera.parallel_projection=True
            self.p.camera.parallel_scale=3.1 if size[0]<700 else 3.1
        self.body_points=np.concatenate(self.points)

    def add(self,mesh,color,**kwargs):
        opts={'ambient':.6}
        opts.update(kwargs)
        actor=self.p.add_mesh(mesh,color=color,reset_camera=False,**opts)
        self.actors.append(actor)
        self.points.append(mesh.points[::max(1,len(mesh.points)//400)])
        return actor

    def water_mesh(self,up):
        # The symmetric tank cut through its centre stays exactly half filled.
        # This plane is AUTHOR-ANIMATED, not an SPH free-surface solution.
        return self.barrel.clip_closed_surface(normal=-np.asarray(up,float),origin=(0,0,0))

    def frame(self,t=None):
        if self.truck:
            # Schematic approach, smoothly decreasing speed, then a full stop.
            stop=2.7
            if t<1.6:
                distance=2.5*t
            elif t<stop:
                u=t-1.6
                distance=4+2.5*(u-u*u/(2*1.1))
            else:
                distance=5.375
            dx=-5.375+distance
            yaw=np.deg2rad(11)*smooth(1.7,2.7,t)
            roll=np.deg2rad(88)*smooth(3.4,5.05,t)
            r=rotation(roll,yaw)
            # Keep the rigid schematic in contact with the road.
            rotated=self.body_points@r.T
            dz=GROUND-rotated[:,2].min()
            dy=.45*smooth(1.7,2.7,t)-1.9*np.sin(roll)
            mat=np.eye(4);mat[:3,:3]=r;mat[:3,3]=(dx,dy,dz)
            for actor in self.actors:
                actor.user_matrix=mat
            sx=.56*smooth(1.65,2.9,t)*(1-smooth(3.25,5.3,t))
            sy=-.55*smooth(2.2,3.55,t)*(1-smooth(4.6,5.7,t))
            normal=r.T@np.array([0.,0.,1.])+np.array([-sx,-sy,0.])
            self.water.mapper.SetInputData(self.water_mesh(normal))
            for lamp in self.lights:
                lamp.GetProperty().SetColor(*(np.array([235,92,72]) if t>1.6 else np.array([99,52,47]))/255)
        self.p.render()
        return Image.fromarray(self.p.screenshot(return_img=True))

    def close(self):
        self.p.close()


def case_slides():
    slides=[]
    small=[]
    for i,(cid,label) in enumerate(CASES):
        scene=Scene((1080,880),case=cid)
        shot=scene.frame();scene.close()
        im=Image.new('RGB',(W,H),BG);im.paste(shot,(0,500));d=ImageDraw.Draw(im)
        text(d,(70,180),'WHAT CHANGES INSIDE?',57,bold=True)
        text(d,(70,1390),f'0{i+1}',38,fill=BLUE)
        text(d,(70,1460),label,52,bold=True,fill=WHITE if i==0 else ORANGE)
        text(d,(70,1770),'PROPOSED TEST GEOMETRY',21,fill='#60798d')
        slides.append(im)
        scene=Scene((510,440),case=cid)
        small.append(scene.frame());scene.close()
        im.save(OUT/f'intro_v3_case_{i+1}.png')
    grid=Image.new('RGB',(W,H),BG);d=ImageDraw.Draw(grid)
    text(d,(70,180),'FOUR DESIGNS.',62,bold=True)
    text(d,(70,265),'ONE TEST.',62,bold=True)
    for i,shot in enumerate(small):
        x=30+(i%2)*510;y=540+(i//2)*590
        grid.paste(shot,(x,y));d=ImageDraw.Draw(grid)
        text(d,(x+20,y-35),f'0{i+1}',28,fill=BLUE)
        label=['NO BAFFLES','COMPARTMENTS','PERFORATED','CROSS BAFFLES'][i]
        text(d,(x+20,y+410),label,28,fill=WHITE if i==0 else ORANGE,bold=True)
    text(d,(70,1770),'SAME TANK. SAME BRAKING.',23,fill='#879dad')
    grid.save(OUT/'intro_v3_four_geometries.png')
    return slides,grid


def action_frame(scene,t):
    im=Image.new('RGB',(W,H),BG);im.paste(scene.frame(t),(0,465));d=ImageDraw.Draw(im)
    label='DRIVE' if t<1.6 else 'BRAKE + SWERVE' if t<3.15 else 'SHIFTING LOAD' if t<4.3 else 'ROLLOVER RISK'
    text(d,(70,180),label,64,bold=True)
    text(d,(70,1770),'CONCEPT ANIMATION',22,fill='#7a93a7')
    return im


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--preview',action='store_true');args=parser.parse_args()
    OUT.mkdir(exist_ok=True)
    geo=export_geometry()
    slides,grid=case_slides()
    scene=Scene((1080,1050),truck=True)
    if args.preview:
        for i,t in enumerate((.7,2.6,3.6,4.5,5.6)):
            action_frame(scene,t).save(OUT/f'intro_v3_motion_preview_{i}.png')
        scene.close();print('Preview frames and four geometry views written.');return
    count=390
    video=OUT/'tanker_intro_v3_motion_portrait.mp4'
    proc=subprocess.Popen(['ffmpeg','-v','error','-y','-f','rawvideo','-pix_fmt','rgb24','-s',f'{W}x{H}','-r',str(FPS),'-i','-',
                           '-an','-c:v','libx264','-preset','fast','-crf','19','-pix_fmt','yuv420p','-movflags','+faststart',str(video)],stdin=subprocess.PIPE)
    try:
        for k in range(count):
            t=k/FPS
            if t<6:
                im=action_frame(scene,t)
            elif t<11.2:
                im=slides[min(3,int((t-6)/1.3))]
            else:
                im=grid
            proc.stdin.write(im.tobytes())
            if k%60==0:
                print(f'frame {k}/{count}',flush=True)
        proc.stdin.close()
        if proc.wait()!=0:raise RuntimeError('Encoder failed')
    finally:
        scene.close()
    report={'status':'Concept animation completed; no SPH or vehicle-dynamics solution',
            'video':video.name,'resolution':[W,H],'fps':FPS,'frames':count,'duration_s':count/FPS,
            'story':[{'time_s':[0,1.6],'label':'DRIVE'},{'time_s':[1.6,3.15],'label':'BRAKE + SWERVE'},
                     {'time_s':[3.15,4.3],'label':'SHIFTING LOAD'},{'time_s':[4.3,6.0],'label':'ROLLOVER RISK'},
                     {'time_s':[6,11.2],'label':'Four proposed study geometries'},{'time_s':[11.2,13],'label':'FOUR DESIGNS. ONE TEST.'}],
            'disclosure':'Truck pose and liquid plane are authored concept animation. No predicted rollover threshold, forces or safety improvement. Intro includes a swerve; the first SPH comparison is prescribed straight braking.',
            'geometry_manifest':'../geometry/concept_v1/geometry.json',
            'technical_reference':'https://www.fmcsa.dot.gov/sites/fmcsa.dot.gov/files/docs/2005%20CDL%20DRIVER%20MANUAL%20FINAL%20July%202010.pdf',
            'source_note':'Driver manual section 8 discusses surge and rapid steering while braking as a rollover risk. This does not validate the authored animation.',
            'checks':{'all_exported_plates_closed':all(p['open_edges']==0 for c in geo['cases'] for p in c['plates']),
                      'half_volume_liquid':'Centre-plane clip of centrally symmetric cylinder; illustrative, not an integrated free-surface solution'}}
    (OUT/'intro_v3_motion.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report),flush=True)


if __name__=='__main__':
    main()
