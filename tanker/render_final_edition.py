"""Final minimal portrait and landscape edits, using the qualified v8 data.

Only presentation timing/composition changes. Native body poses, particle fields,
loads, and numerical setups are unchanged. Intro runs the entire8s record at2x.
"""
from pathlib import Path
from functools import lru_cache
from io import BytesIO
import argparse,json,subprocess
import numpy as np
import pyvista as pv
from PIL import Image,ImageDraw
from scipy.integrate import cumulative_trapezoid
from matplotlib import colormaps
from matplotlib import rc_context
from matplotlib.mathtext import math_to_image
from matplotlib.font_manager import FontProperties
from render_minimal_intro import BG,WHITE,BLUE,ORANGE,text,arrow
from render_intro_motion import Scene,CASES
from render_braking_results import surface
from render_coupled_results import CoupledScene,prepare,LABELS
from render_load_comparison import LoadScene,CLIM

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'media';FPS=30
NAMES=[name+'_coupled_lat690_v8' for name,_ in CASES]


@lru_cache(None)
def equation(value,width,color=WHITE):
    buffer=BytesIO()
    with rc_context({'savefig.transparent':True,'savefig.facecolor':'none','figure.facecolor':'none'}):
        math_to_image('$'+value+'$',buffer,prop=FontProperties(size=29),dpi=180,format='png',color=color)
    im=Image.open(buffer).convert('RGBA')
    # Retain true glyph alpha edges without a rectangular background.
    box=im.getbbox()
    if box:im=im.crop(box)
    if im.width>width:im=im.resize((width,round(im.height*width/im.width)),Image.Resampling.LANCZOS)
    return im


def paste_equation(im,value,centre_y,x,width,color=WHITE):
    eq=equation(value,width,color);im.paste(eq,(x+(width-eq.width)//2,round(centre_y-eq.height/2)),eq)


class Edition:
    def __init__(self,orientation):
        self.orientation=orientation;self.land=orientation=='landscape'
        self.w,self.h=(1920,1080) if self.land else (1080,1920)
        self.folders=[ROOT/'sph'/name for name in NAMES]
        self.records=[prepare(f) for f in self.folders]
        self.setup=[json.loads((f/'setup.json').read_text()) for f in self.folders]
        self.summary=json.loads((ROOT/'comparison_summary.json').read_text())
        assert all(s['fluid_particles']==137812 for s in self.setup)
        acceleration=np.loadtxt(self.folders[0]/'input_acceleration.csv',delimiter=',',skiprows=1)
        self.time=acceleration[:,0]
        self.speed=np.maximum(0,50/3.6+cumulative_trapezoid(acceleration[:,1],self.time,initial=0))
        self.distance=cumulative_trapezoid(self.speed,self.time,initial=0)
        self.sections=[]

    def canvas(self):return Image.new('RGB',(self.w,self.h),BG)

    def write_clip(self,name,duration,frame,previews=()):
        path=OUT/f'tanker_final_{self.orientation}_{name}.mp4';n=round(duration*FPS)
        proc=subprocess.Popen(['ffmpeg','-v','error','-y','-f','rawvideo','-pix_fmt','rgb24','-s',f'{self.w}x{self.h}',
            '-r',str(FPS),'-i','-','-an','-c:v','libx264','-preset','fast','-crf','19',
            '-pix_fmt','yuv420p','-movflags','+faststart',str(path)],stdin=subprocess.PIPE)
        try:
            for k in range(n):
                im=frame(k/(n-1) if n>1 else 0.)
                assert im.size==(self.w,self.h)
                if k in previews:im.save(OUT/f'tanker_final_{self.orientation}_{name}_{k:03d}.png')
                proc.stdin.write(im.tobytes())
                if k%90==0:print(f'{self.orientation} {name} {k}/{n}',flush=True)
        finally:proc.stdin.close()
        if proc.wait():raise RuntimeError(f'Encoding {path.name}')
        self.sections.append(dict(name=name,file=path.name,seconds=n/FPS))

    def actual(self,scene,case,t):
        pose,cache=self.records[case];index=int(np.argmin(abs(pose['time_s']-t)))
        image=scene.actual_frame(pose['rotation'][index],pose['center'][index],
            pv.read(cache/f'water_{index:04d}.vtp'),np.interp(t,self.time,self.distance),t>=1.5)
        return image,float(pose['tilt_deg'][index])

    def baseline(self,intro=False):
        if self.land:
            size,xy,scale=((1450,1080),(470,0),4.65) if intro else ((1920,900),(0,160),3.2)
        else:size,xy,scale=(1080,1240),(0,425),5.8
        scene=CoupledScene(size,CASES[0][0]);scene.p.camera.parallel_scale=scale
        def frame(q):
            t=8*q;im=self.canvas();shot,tilt=self.actual(scene,0,t);im.paste(shot,xy);d=ImageDraw.Draw(im)
            if self.land and intro:
                text(d,(65,190),'BRAKE.',65,bold=True);text(d,(65,280),'SWERVE.',65,bold=True)
                text(d,(65,435),'8.8 t',84,bold=True,fill=BLUE);text(d,(65,535),'OF WATER',36,fill=BLUE)
                text(d,(65,700),'50 → 0 km/h',36,fill='#9aafbe')
                text(d,(65,900),'2×',32,fill='#9aafbe')
            elif self.land:
                text(d,(65,60),'NO BAFFLES',56,bold=True)
                text(d,(1650,65),f'{t:.1f} s',40,fill=BLUE)
                text(d,(65,980),f'{tilt:.0f}°',40,fill='#9aafbe')
            else:
                text(d,(65,150),'BRAKE + SWERVE' if intro else 'NO BAFFLES',65,bold=True)
                if intro:
                    text(d,(65,245),'8.8 t OF WATER',58,bold=True,fill=BLUE)
                    text(d,(65,345),'50 → 0 km/h',34,fill='#9aafbe')
                    text(d,(920,345),'2×',30,fill='#9aafbe')
                else:
                    text(d,(65,285),f'{t:.1f} s',36,fill=BLUE)
                    text(d,(65,1710),f'{tilt:.0f}°',40,fill='#9aafbe')
            return im
        duration=4 if intro else 6
        self.write_clip('intro' if intro else 'baseline',duration,frame,(0,round(.5*(duration*30-1)),duration*30-1))
        scene.close()

    def geometry(self):
        im=self.canvas();d=ImageDraw.Draw(im)
        if self.land:text(d,(65,75),'FOUR DESIGNS. SAME INPUT.',62,bold=True)
        else:
            text(d,(65,150),'FOUR DESIGNS.',65,bold=True);text(d,(65,235),'SAME INPUT.',65,bold=True)
        for k,(layout,_) in enumerate(CASES):
            suffix='brake_040_v1' if k==0 else 'thick_040_v2'
            source=ROOT/'sph'/f'{layout}_{suffix}'
            size=(460,450) if self.land else (510,440)
            scene=Scene(size,case=layout)
            for j,actor in enumerate(scene.actors[6:-1]):
                mesh=actor.mapper.dataset.copy(deep=True)
                if layout=='D_cross_baffles' and j==2:mesh.points[:,1]*=.16/.06
                else:
                    cx=(-1.,1.)[j];mesh.points[:,0]=cx+(mesh.points[:,0]-cx)*.16/.06
                actor.mapper.SetInputData(mesh)
            scene.water.mapper.SetInputData(surface(pv.read(source/'fluid.vtk').points))
            scene.water.GetProperty().SetOpacity(.78)
            if self.land:scene.p.camera.parallel_scale=3.3
            x,y=(30+470*k,290) if self.land else (30+510*(k%2),520+570*(k//2))
            im.paste(scene.frame(),(x,y));scene.close();d=ImageDraw.Draw(im)
            text(d,(x+20,y+430),LABELS[layout],27,bold=True,fill=WHITE if k==0 else ORANGE)
        im.save(OUT/f'tanker_final_{self.orientation}_geometry.png')
        self.write_clip('geometry',2.5 if self.land else 3,lambda q:im)

    def method(self):
        source=ROOT/'sph/A_unpartitioned_brake_040_v1';im=self.canvas();d=ImageDraw.Draw(im)
        text(d,(65,75 if self.land else 140),'137,812',80 if self.land else 94,bold=True,fill=BLUE)
        text(d,(430,120),'FLUID PARTICLES',32) if self.land else text(d,(65,255),'FLUID PARTICLES',36)
        scene=Scene((1130,720) if self.land else (1080,690));scene.water.SetVisibility(False)
        rng=np.random.default_rng(81)
        for label,color,count,size in [('fluid',BLUE,6000,4),('walls',ORANGE,1600,3)]:
            points=pv.read(source/f'{label}.vtk').points
            scene.p.add_points(points[rng.choice(len(points),count,replace=False)],color=color,
                point_size=size,render_points_as_spheres=False,lighting=False)
        im.paste(scene.frame(),(0,250 if self.land else 400));scene.close();d=ImageDraw.Draw(im)
        # Matches the local solver's symmetric pressure term; remaining fluid
        # viscosity/body accelerations are grouped, not silently omitted.
        fluid=r'\dot{\mathbf{v}}_i=-\sum_j m_j\frac{p_i+p_j}{\rho_i\rho_j}\nabla_i W_{ij}+\mathbf{a}_{\mathrm{visc}}+\mathbf{g}-\mathbf{a}_0'
        translation=r'm\ddot{\mathbf{x}}=\sum\mathbf{F}'
        rotation=r'\mathbf{I}\dot{\boldsymbol{\omega}}+\boldsymbol{\omega}\times(\mathbf{I}\boldsymbol{\omega})=\sum\boldsymbol{\tau}'
        if self.land:
            text(d,(1190,270),'FLUID',30,bold=True,fill=BLUE)
            paste_equation(im,fluid,380,1170,690)
            d=ImageDraw.Draw(im);arrow(d,(1490,475),(1490,570),BLUE,7);arrow(d,(1560,570),(1560,475),ORANGE,7)
            text(d,(1190,650),'RIGID BODY',30,bold=True,fill=ORANGE)
            paste_equation(im,translation,750,1170,690)
            paste_equation(im,rotation,880,1170,690)
        else:
            text(d,(65,1120),'FLUID',30,bold=True,fill=BLUE)
            paste_equation(im,fluid,1230,65,950)
            d=ImageDraw.Draw(im);arrow(d,(500,1340),(500,1410),BLUE,7);arrow(d,(575,1410),(575,1340),ORANGE,7)
            text(d,(65,1480),'RIGID BODY',30,bold=True,fill=ORANGE)
            paste_equation(im,translation,1600,65,950)
            paste_equation(im,rotation,1730,65,950)
        im.save(OUT/f'tanker_final_{self.orientation}_method.png')
        self.write_clip('method',2.5,lambda q:im)

    def panel_layout(self,load=False):
        if self.land:return (930,365),[(30+960*(k%2),145+400*(k//2)) for k in range(4)],3.7 if load else 3.05
        return (510,540),[(30+510*(k%2),430+600*(k//2)) for k in range(4)],5.8

    def comparison(self):
        size,positions,scale=self.panel_layout()
        scenes=[CoupledScene(size,c[0],s.get('numerical_plate_thickness_m',.08)) for c,s in zip(CASES,self.setup)]
        for scene in scenes:scene.p.camera.parallel_scale=scale
        def frame(q):
            t=8*q;im=self.canvas();d=ImageDraw.Draw(im)
            text(d,(65,55 if self.land else 150),'SAME INPUT.',56 if self.land else 65,bold=True)
            text(d,(1680,65) if self.land else (65,265),f'{t:.1f} s',36,fill=BLUE)
            for k,(scene,(x,y)) in enumerate(zip(scenes,positions)):
                shot,tilt=self.actual(scene,k,t);im.paste(shot,(x,y));d=ImageDraw.Draw(im)
                text(d,(x+20,y+10),LABELS[CASES[k][0]],28,bold=True,fill=WHITE if k==0 else ORANGE)
                text(d,(x+20,y+size[1]-45),f'{tilt:.0f}°',31,fill='#9aafbe')
            return im
        self.write_clip('comparison',8,frame,(0,90,120,239))
        for scene in scenes:scene.close()

    def loads(self):
        size,positions,scale=self.panel_layout(True)
        scenes=[LoadScene(f,size) for f in self.folders]
        for scene in scenes:scene.p.camera.parallel_scale=scale
        def frame(q):
            t=8*q;im=self.canvas();d=ImageDraw.Draw(im)
            text(d,(65,55 if self.land else 150),'FORCES + MOMENTS',56 if self.land else 58,bold=True)
            text(d,(1680,65) if self.land else (65,265),f'{t:.1f} s',36,fill='#9aafbe')
            for k,(scene,(x,y),(pose,cache)) in enumerate(zip(scenes,positions,self.records)):
                i=int(np.argmin(abs(pose['time_s']-t)))
                shot,f,m=scene.load_frame(i,pose,pv.read(cache/f'water_{i:04d}.vtp'),t)
                im.paste(shot,(x,y));d=ImageDraw.Draw(im)
                text(d,(x+20,y+10),LABELS[CASES[k][0]],28,bold=True)
                if self.land:
                    text(d,(x+20,y+size[1]-38),f'|F| {f:.0f} kN',30,fill=BLUE)
                    text(d,(x+490,y+size[1]-38),f'|M| {m:.0f} kNm',30,fill=ORANGE)
                else:
                    text(d,(x+20,y+size[1]-58),f'|F| {f:.0f} kN',29,fill=BLUE)
                    text(d,(x+20,y+size[1]-15),f'|M| {m:.0f} kNm',29,fill=ORANGE)
            bx,by,bw=(360,970,1200) if self.land else (130,1695,820)
            bar=np.uint8(colormaps['jet'](np.linspace(0,1,bw))[None,:,:3]*255)
            im.paste(Image.fromarray(bar).resize((bw,18)),(bx,by));d=ImageDraw.Draw(im)
            for value in [0,10,100,1000,CLIM]:
                x=bx+bw*np.log10(1+value)/np.log10(1+CLIM)
                text(d,(x,by+28),f'{value:.0f}',22,fill='#9aafbe',anchor='mt')
            text(d,(bx,by+70),'N / PARTICLE  ·  LOG',23,fill='#9aafbe')
            return im
        self.write_clip('loads',8,frame,(0,90,120,210))
        for scene in scenes:scene.close()

    def conclusion(self):
        records=self.summary['cases'];im=self.canvas();d=ImageDraw.Draw(im)
        maxmoment=max(r['peak_roll_moment_kNm'] for r in records)
        reduction=100*(1-records[3]['peak_roll_moment_kNm']/records[0]['peak_roll_moment_kNm'])
        text(d,(65,65 if self.land else 150),'CROSS BAFFLES',62 if self.land else 65,bold=True)
        text(d,(65,155 if self.land else 245),'STAYED UPRIGHT',43 if self.land else 57,bold=True,fill=BLUE)
        text(d,(65,285 if self.land else 480),'PEAK ROLL MOMENT [kNm]',27,fill=ORANGE)
        text(d,(840,285 if self.land else 480),'75° AT',27,fill='#9aafbe')
        for k,r in enumerate(records):
            y=(380+145*k) if self.land else (600+205*k)
            text(d,(65,y),r['label'],31,bold=True)
            color=BLUE if k==3 else ORANGE;bw=490*r['peak_roll_moment_kNm']/maxmoment
            d.rounded_rectangle((65,y+55,65+bw,y+85),radius=5,fill=color)
            text(d,(590,y+45),f"{r['peak_roll_moment_kNm']:.1f}",36,bold=True,fill=color)
            value=f"{r['first_75deg_s']:.2f} s" if r['first_75deg_s'] is not None else 'NONE'
            text(d,(840,y+45),value,33,bold=True,fill=color)
        if self.land:
            scene=CoupledScene((800,490),CASES[3][0],.16);scene.p.camera.parallel_scale=3.7
            shot,_=self.actual(scene,3,3.5);scene.close();im.paste(shot,(1100,220));d=ImageDraw.Draw(im)
            text(d,(1270,740),f'−{reduction:.0f}%',108,bold=True,fill=BLUE)
            text(d,(1230,875),'PEAK ROLL MOMENT',33,fill=BLUE)
            text(d,(65,1005),'1.5–3.5 s',26,fill='#9aafbe')
        else:
            text(d,(65,1470),f'−{reduction:.0f}%',110,bold=True,fill=BLUE)
            text(d,(65,1600),'PEAK ROLL MOMENT',37,fill=BLUE)
            text(d,(65,1750),'1.5–3.5 s',26,fill='#9aafbe')
        im.save(OUT/f'tanker_final_{self.orientation}_conclusion.png')
        self.write_clip('conclusion',5,lambda q:im)

    def assemble(self):
        listing=OUT/f'final_{self.orientation}_concat.txt'
        listing.write_text(''.join(f"file '{s['file']}'\n" for s in self.sections),encoding='utf-8')
        title='Tanker sloshing: four designs, one input'
        metadata=[';FFMETADATA1','title='+title];t=0.
        for s in self.sections:
            s['start_s']=t;s['end_s']=t+s['seconds']
            metadata.extend(['[CHAPTER]','TIMEBASE=1/1000',f'START={round(t*1000)}',f'END={round(s["end_s"]*1000)}',f'title={s["name"]}'])
            t=s['end_s']
        meta=OUT/f'final_{self.orientation}_chapters.txt';meta.write_text('\n'.join(metadata),encoding='utf-8')
        video=OUT/f'tanker_sloshing_final_{self.orientation}.mp4'
        subprocess.run(['ffmpeg','-v','error','-y','-f','concat','-safe','0','-i',listing.name,'-i',meta.name,
            '-map','0:v:0','-map_metadata','1','-map_chapters','1','-c','copy','-movflags','+faststart',video.name],cwd=OUT,check=True)
        subprocess.run(['ffmpeg','-v','error','-i',str(video),'-f','null','-'],check=True)
        result=dict(video=video.name,width=self.w,height=self.h,fps=FPS,seconds=t,sections=self.sections,
            input_data='Qualified v8 data, unchanged. ../comparison_summary.json',
            pacing='Intro entire0-8s at2x; single replay0-8s in6s; comparison and loads0-8s in8s.',
            method='137812 real fluid particles;6000 displayed initial fluid points and1600 wall points. Fluid symmetric pressure-gradient equation plus viscosity, gravity, reference-frame acceleration. Rigid linear/angular momentum balances.',
            formula_source='work/DualSPHysics-source/src/source/JSphGpu_ker.cu:pressure term nearline594. Rigid-body sums include water/contact/gravity/frame forces; torque equation in body axes.',
            visual_changes='Removed method banners, repeated setup descriptions and caveat text; kept case labels, physical time, units, contour scaling and comparison window.',
            load_reference='Moment about dry-vehicleCG;0.4s filtered resultants. Per-boundary-particleN is not pressure. Common log10(1+load)0-5000N jet.',
            limitations=self.summary['limitations'],verification='Full encoded video decoded without error.')
        (OUT/f'tanker_sloshing_final_{self.orientation}.json').write_text(json.dumps(result,indent=2))
        print(json.dumps({'video':str(video),'seconds':t,'size':[self.w,self.h]}),flush=True)

    def run(self,preview=False):
        if preview:
            # Short preliminary render validates both mathematical and cinematic layouts.
            original=self.write_clip
            def still(name,duration,frame,previews=()):
                for q in ([0,.5,.9] if name in ['intro','baseline','comparison','loads'] else [0]):
                    frame(q).save(OUT/f'preview_final_{self.orientation}_{name}_{q:.1f}.png')
            self.write_clip=still
        self.baseline(True);self.geometry();self.method();self.baseline();self.comparison();self.loads();self.conclusion()
        if not preview:self.assemble()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('orientation',choices=['portrait','landscape']);p.add_argument('--preview',action='store_true')
    args=p.parse_args();Edition(args.orientation).run(args.preview)
