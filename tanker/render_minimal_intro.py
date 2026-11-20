"""Minimal concept intro. Original schematic geometry; no fluid results."""
from pathlib import Path
import json,os
import subprocess
from functools import lru_cache

import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'media'
W, H, FPS = 1080, 1920, 30
BG = '#080d16'
WHITE = '#e5edf3'
BLUE = '#4eabd4'
ORANGE = '#e9ad61'


@lru_cache(None)
def font(size, bold=False):
    return ImageFont.truetype(str(Path(os.environ.get('TANKER_FONT_DIR','C:/Windows/Fonts'))/('arialbd.ttf' if bold else 'arial.ttf')), size)


def text(draw, xy, value, size, fill=WHITE, bold=False, anchor=None):
    draw.text(xy, value, font=font(size, bold), fill=fill, anchor=anchor)


def arrow(draw, start, end, color, width=7):
    a, b = np.array(start, float), np.array(end, float)
    direction = (b-a) / np.linalg.norm(b-a)
    side = np.array([-direction[1], direction[0]])
    draw.line((tuple(a), tuple(b)), fill=color, width=width)
    draw.polygon([tuple(b), tuple(b-direction*23+side*12), tuple(b-direction*23-side*12)], fill=color)


def perforated_plate(x):
    # A visual plate with explicit open passages. Not a solver mesh.
    n = 75
    q = np.linspace(-1.0, 1.0, n)
    yy, zz = np.meshgrid(q, q, indexing='ij')
    points = np.column_stack((np.full(n*n, x), yy.ravel(), zz.ravel()))
    faces = []
    for i in range(n-1):
        for j in range(n-1):
            y, z = (q[i]+q[i+1])/2, (q[j]+q[j+1])/2
            if y*y+z*z > .97**2:
                continue
            if any((y-cy)**2+(z-cz)**2 < .24**2 for cy in [-.42,.42] for cz in [-.40,.40]):
                continue
            k = i*n+j
            faces.extend([4, k, k+n, k+n+1, k+1])
    return pv.PolyData(points, np.array(faces)).clean().triangulate().extrude((.055,0,0), capping=True)


def tank_image(baffles=False, truck=False, size=(1080,650)):
    p = pv.Plotter(off_screen=True, window_size=size)
    p.set_background(BG)
    p.enable_anti_aliasing('ssaa')
    p.enable_depth_peeling(number_of_peels=12, occlusion_ratio=0)
    barrel = pv.Cylinder(center=(0,0,0), direction=(1,0,0), radius=1, height=6, resolution=128).triangulate()
    p.add_mesh(barrel, color='#b6cad8', opacity=.055, smooth_shading=True, ambient=.7)
    water = pv.Cylinder(center=(0,0,0), direction=(1,0,0), radius=.965, height=5.95, resolution=128).triangulate().clean()
    water = water.clip_closed_surface(normal=(0,0,1), origin=(0,0,0))
    # clip_closed_surface keeps the positive side, so reflect to get the lower half.
    if water.center[2] > 0:
        water.points[:,2] *= -1
        water = water.compute_normals(auto_orient_normals=True)
    p.add_mesh(water, color='#247eac', opacity=.60, smooth_shading=True, ambient=.55, specular=.2)
    th = np.linspace(0,2*np.pi,150)
    for x in (-3,3):
        rim = np.column_stack((np.full(th.shape,x),np.cos(th),np.sin(th)))
        p.add_mesh(pv.lines_from_points(rim, close=True), color='#a0b6c7', line_width=3, lighting=False)
    for theta in (np.pi/2, 3*np.pi/2, np.pi):
        p.add_mesh(pv.Line((-3,np.cos(theta),np.sin(theta)),(3,np.cos(theta),np.sin(theta))), color='#60798d', line_width=2, lighting=False)
    if baffles:
        for x in (-1,1):
            p.add_mesh(perforated_plate(x), color=ORANGE, ambient=.65, diffuse=.4)
    if truck:
        p.add_mesh(pv.Box(bounds=(-3.6,4.7,-.8,.8,-1.35,-1.12)),color='#718491',ambient=.6)
        p.add_mesh(pv.Box(bounds=(3.25,4.65,-.86,.86,-1.12,.65)),color='#8d9ca7',ambient=.7)
        p.add_mesh(pv.Box(bounds=(3.45,4.58,-.88,-.87,-.1,.48)),color='#193a50',ambient=.9)
        for x in (-2.1,-.8,3.8):
            for y in (-.92,.92):
                p.add_mesh(pv.Cylinder(center=(x,y,-1.55),direction=(0,1,0),radius=.48,height=.22,resolution=48),color='#24303c',ambient=.65)
                p.add_mesh(pv.Cylinder(center=(x,y+np.sign(y)*.12,-1.55),direction=(0,1,0),radius=.25,height=.025,resolution=40),color='#8396a4',ambient=.65)
        p.add_mesh(pv.Line((-4.4,-1.3,-2.06),(5.4,-1.3,-2.06)),color='#2f4354',line_width=2,lighting=False)
        target = (.6,0,-.3)
        p.camera_position = [(10,-20,8),target,(0,0,1)]
        p.camera.parallel_projection = True
        p.camera.parallel_scale = 3.3
    else:
        p.camera_position = [(9,-20,8),(0,0,0),(0,0,1)]
        p.camera.parallel_projection = True
        p.camera.parallel_scale = 2.25
    result = Image.fromarray(p.screenshot(return_img=True))
    p.close()
    return result


def make_slides():
    OUT.mkdir(exist_ok=True)
    truck = tank_image(truck=True,size=(1080,760))
    open_tank = tank_image(size=(1080,560))
    baffled = tank_image(baffles=True,size=(1080,560))

    first = Image.new('RGB',(W,H),BG)
    first.paste(truck,(0,530))
    d = ImageDraw.Draw(first)
    text(d,(70,185),'SUDDEN BRAKING',64,bold=True)
    text(d,(70,340),'DECELERATION',28,fill='#a3b7c6')
    arrow(d,(790,385),(260,385),ORANGE)
    text(d,(70,1330),'LIQUID SURGE',28,fill=BLUE)
    arrow(d,(260,1380),(790,1380),BLUE)
    text(d,(70,1770),'CONCEPT',21,fill='#60798d')

    second = Image.new('RGB',(W,H),BG)
    second.paste(open_tank,(0,400))
    second.paste(baffled,(0,1050))
    d = ImageDraw.Draw(second)
    text(d,(70,185),'CAN BAFFLES HELP?',64,bold=True)
    text(d,(70,365),'NO BAFFLES',30,fill='#9fb3c3',bold=True)
    text(d,(70,1015),'WITH BAFFLES',30,fill=ORANGE,bold=True)
    text(d,(70,1770),'CONCEPT',21,fill='#60798d')
    first.save(OUT/'minimal_intro_v2_braking.png')
    second.save(OUT/'minimal_intro_v2_compare.png')
    return first, second


def encode(slides):
    frames = 225
    path = OUT/'tanker_minimal_intro_v2_portrait.mp4'
    cmd=['ffmpeg','-v','error','-y','-f','rawvideo','-pix_fmt','rgb24','-s',f'{W}x{H}','-r',str(FPS),'-i','-',
         '-an','-c:v','libx264','-preset','fast','-crf','19','-pix_fmt','yuv420p','-movflags','+faststart',str(path)]
    proc=subprocess.Popen(cmd,stdin=subprocess.PIPE)
    for k in range(frames):
        t=k/FPS
        if t<3.1:
            im=slides[0]
        elif t<3.4:
            im=Image.blend(slides[0],slides[1],(t-3.1)/.3)
        else:
            im=slides[1]
        proc.stdin.write(im.tobytes())
    proc.stdin.close()
    if proc.wait()!=0:
        raise RuntimeError('Video encoding failed')
    data={
        'status':'Editorial concept only. No SPH calculation started.',
        'video':str(path.name),'resolution':[W,H],'fps':FPS,'frames':frames,'duration_s':frames/FPS,
        'source':'Original procedural schematic geometry; no generated truck photographs used.',
        'physical_scope':'Prescribed straight braking. Arrows represent braking acceleration and initial relative liquid surge direction, not measured velocity or force.',
        'fluid_visual':'Static horizontal initial fill; no synthetic transient results.',
        'story':[{'time_s':[0,3.1],'text':'SUDDEN BRAKING','purpose':'Opposite braking and relative surge directions.'},
                 {'time_s':[3.4,7.5],'text':'CAN BAFFLES HELP?','purpose':'Same tank and initial fill, with only two perforated internal plates added.'}],
        'note':'The full study still includes four layouts. Intro shows two to explain the question, not a ranking.'
    }
    (OUT/'minimal_intro_v2.json').write_text(json.dumps(data,indent=2),encoding='utf-8')
    print(json.dumps(data))


if __name__=='__main__':
    encode(make_slides())
