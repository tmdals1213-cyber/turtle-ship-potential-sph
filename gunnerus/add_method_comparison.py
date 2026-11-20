"""Insert a conceptual FVM / IBM / SPH comparison into the completed film.

This edits presentation only. The CFD replay comes from the existing video.
"""
from pathlib import Path
import os
import hashlib
import json
import re
import subprocess

import numpy as np
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'media'
SOURCE = OUT / 'gunnerus_stronger_waves_story_portrait.mp4'
TARGET = OUT / 'gunnerus_stronger_waves_methods_portrait.mp4'
BG = '#080d16'
BLUE = '#64cbe7'
GOLD = '#e3b75e'


def text(draw, pos, value, size, color='#e6eef3', bold=False):
    font = ImageFont.truetype((os.environ.get('SHIP_FONT_DIR', 'C:/Windows/Fonts').rstrip('/\\') + '/') + ('segoeuib.ttf' if bold else 'segoeui.ttf'), size)
    draw.text(pos, value, font=font, fill=color)


def diagram(kind):
    w, h = 936, 240
    picture = Image.new('RGB', (w, h), '#101b28')
    d = ImageDraw.Draw(picture)
    x = np.linspace(0, w - 1, w)
    surface = 110 + 24 * np.sin(x / 120) - 42 * np.exp(-((x - 710) / 100)**2)
    water = [(int(a), int(b)) for a, b in zip(x, surface)] + [(w - 1, h - 1), (0, h - 1)]
    if kind != 'sph':
        d.polygon(water, fill='#12405c')
        d.line(water[:-2], fill=BLUE, width=3)
    angle = np.linspace(0, 2 * np.pi, 80, endpoint=False)
    inner = np.column_stack((165 * np.cos(angle), np.where(np.sin(angle) > 0, 52, 17) * np.sin(angle)))
    a = np.radians(-12)
    rotation = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    center = np.array([w / 2, 108.])
    inner = inner @ rotation.T + center
    hull = [tuple(v) for v in inner]
    solid = Image.new('1', (w, h), 0)
    ImageDraw.Draw(solid).polygon(hull, fill=1)
    mask = np.asarray(solid)

    if kind == 'fvm':
        # An illustrative O-grid: each inner edge conforms to the hull polygon.
        rays = inner - center
        outer = []
        for rx, ry in rays:
            sx = ((w - 1 - center[0]) if rx > 0 else -center[0]) / rx if abs(rx) > 1e-8 else np.inf
            sy = ((h - 1 - center[1]) if ry > 0 else -center[1]) / ry if abs(ry) > 1e-8 else np.inf
            outer.append(center + min(sx, sy) * np.array([rx, ry]))
        outer = np.array(outer)
        for f in np.linspace(0, 1, 10)**1.35:
            ring = inner * (1 - f) + outer * f
            d.line([tuple(v) for v in np.vstack((ring, ring[0]))], fill='#395166', width=1)
        for v, end in zip(inner, outer):
            d.line([tuple(v), tuple(end)], fill='#395166', width=1)
    elif kind == 'ibm':
        for xx in range(0, w, 27):
            d.line((xx, 0, xx, h), fill='#395166', width=1)
        for yy in range(0, h, 27):
            d.line((0, yy, w, yy), fill='#395166', width=1)
        # Highlight the cells cut by the immersed surface.
        cells = {(int(px // 27) * 27, int(py // 27) * 27) for px, py in inner}
        for xx, yy in cells:
            d.rectangle((xx, yy, xx + 27, yy + 27), outline='#ad8443', width=2)
    else:
        # Conceptual parcels, explicitly separated from the actual CFD replay.
        for yy in range(8, h, 15):
            for xx in range(8, w, 15):
                if yy > surface[xx] and not mask[yy, xx]:
                    d.ellipse((xx - 3, yy - 3, xx + 3, yy + 3), fill=BLUE)
        for xx, yy in [(730, 40), (756, 52), (775, 32), (795, 48)]:
            d.ellipse((xx - 4, yy - 4, xx + 4, yy + 4), fill=BLUE)
        for xx, yy in [(115, 184), (740, 178)]:
            d.line((xx, yy, xx - 62, yy - 10), fill='#d8edf6', width=3)
            d.polygon([(xx - 62, yy - 10), (xx - 48, yy - 15), (xx - 51, yy)], fill='#d8edf6')
    d.polygon(hull, fill='#675334')
    d.line(hull + [hull[0]], fill=GOLD, width=3)
    if kind == 'sph':
        for px, py in inner[::2]:
            d.ellipse((px - 3, py - 3, px + 3, py + 3), fill=GOLD)
    d.rectangle((0, 0, w - 1, h - 1), outline='#24394a', width=2)
    return picture


def main():
    frame = Image.new('RGB', (1080, 1920), BG)
    draw = ImageDraw.Draw(frame)
    text(draw, (60, 58), 'WHY SPH HERE?', 58, bold=True)
    text(draw, (60, 150), 'Breaking waves + large body motion', 29, '#a8c4d4')
    rows = [('fvm', 'BODY-FITTED FVM', 'Conforming mesh + VOF'),
            ('ibm', 'FVM + IBM', 'Fixed grid + immersed hull + VOF'),
            ('sph', 'SPH', 'Particles follow the water')]
    for i, (kind, title, caption) in enumerate(rows):
        top = 285 + i * 435
        text(draw, (72, top), title, 38, BLUE if kind == 'sph' else '#e6eef3', True)
        text(draw, (72, top + 61), caption, 27, '#a8c4d4')
        frame.paste(diagram(kind), (72, top + 115))
    text(draw, (60, 1655), 'IBM handles solid boundaries.', 25, '#91a9ba')
    text(draw, (60, 1700), 'VOF captures the water-air interface.', 25, '#91a9ba')
    text(draw, (60, 1810), 'Conceptual comparison  /  Not a solver benchmark', 22, '#627f92')
    slide = OUT / 'gunnerus_method_comparison.png'
    frame.save(slide)
    # Decode source frames exactly once. Preserve all original 660 frames in order,
    # adding a 120-frame slide after the 60-frame opening setup.
    filters = '[0:v]split=2[p][q];[p]trim=end_frame=60,setpts=PTS-STARTPTS[a];[q]trim=start_frame=60,setpts=PTS-STARTPTS[b];[1:v]trim=end_frame=120,setpts=PTS-STARTPTS,format=yuv420p[c];[a][c][b]concat=n=3:v=1:a=0[v]'
    subprocess.run(['ffmpeg', '-v', 'error', '-y', '-i', str(SOURCE), '-loop', '1', '-framerate', '30', '-i', str(slide),
                    '-filter_complex', filters, '-map', '[v]', '-frames:v', '780', '-an', '-c:v', 'libx264',
                    '-preset', 'fast', '-crf', '18', '-pix_fmt', 'yuv420p', '-movflags', '+faststart', str(TARGET)], check=True)
    result = subprocess.run(['ffmpeg', '-v', 'error', '-nostats', '-progress', 'pipe:1', '-i', str(TARGET), '-f', 'null', '-'],
                            capture_output=True, text=True, check=True)
    frames = [int(n) for n in re.findall(r'^frame=(\d+)', result.stdout, re.M)]
    assert frames[-1] == 780
    metadata = json.loads(SOURCE.with_suffix('.json').read_text())
    metadata.update(frames=780, video_seconds=26, source_video=SOURCE.name,
                    segments=[{'kind':'domain setup','video_s':[0,2]},
                              {'kind':'conceptual method comparison','video_s':[2,6]},
                              {'kind':'SPH schematic','video_s':[6,8]},
                              {'kind':'actual simulation replay','video_s':[8,26],'model_time_s':[4,16]}],
                    comparison_note='Illustrative discretizations, not three CFD results. IBM can be combined with FVM. No solver ranking or runtime comparison is implied.',
                    sources=['https://www.openfoam.com/documentation/user-guide/a-reference/a.1-standard-solvers',
                             'https://engineering.jhu.edu/fsag/wp-content/uploads/2014/06/jcp-vol174.pdf',
                             'https://github.com/DualSPHysics/DualSPHysics/wiki/3.-SPH-formulation'],
                    decoded_frames=frames[-1],sha256=hashlib.sha256(TARGET.read_bytes()).hexdigest(),
                    final_visual_review='pending')
    TARGET.with_suffix('.json').write_text(json.dumps(metadata,indent=2))
    for t in [1.9,2.1,5.9,6.1,8.1,14.5,25.9]:
        subprocess.run(['ffmpeg','-v','error','-y','-ss',str(t),'-i',str(TARGET),'-frames:v','1','-update','1',
                        str(OUT/f'gunnerus_methods_check_{t:04.1f}.png')],check=True)
    print(json.dumps({'video':str(TARGET),'frames':frames[-1],'sha256':metadata['sha256']}),flush=True)


if __name__ == '__main__':
    main()
