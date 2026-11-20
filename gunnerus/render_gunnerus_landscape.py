"""Landscape edit of the completed Gunnerus run; no new simulation."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import subprocess

import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw
from matplotlib import colormaps

import render_gunnerus_focus_story as story
from render_gunnerus_comparison import text, BG, ROOT, OUT
from add_method_comparison import diagram

CASE = 'gunnerus_focus_A115_x5.5_D300_050_16s'
NAME = 'gunnerus_stronger_waves_methods_landscape'
SIZE = (1920, 1080)
FPS = 30
story.NAME = CASE
story.FOLDER = ROOT / 'sph' / CASE


def intro_frames():
    setup = Image.new('RGB', SIZE, BG)
    d = ImageDraw.Draw(setup)
    text(d, (60, 35), 'GUNNERUS', 55, bold=True)
    text(d, (60, 110), 'FOCUSED HEAD WAVES', 30, '#64cbe7', True)
    basin = Image.open(OUT / 'gunnerus_stronger_setup.png').crop((0, 330, 1080, 1030))
    setup.paste(basin.resize((1170, 758), Image.Resampling.LANCZOS), (0, 190))
    text(d, (1220, 210), '360 x 100 m', 49, bold=True)
    text(d, (1220, 278), '30 m depth', 35, '#a8c4d4')
    text(d, (1220, 337), 'Full-scale equivalent', 25, '#91a9ba')
    text(d, (1220, 450), 'H ~ 13-15 m', 49, '#64cbe7', True)
    text(d, (1220, 519), 'Measured in the empty tank', 27, '#91a9ba')
    text(d, (1220, 562), 'Crest-to-trough / equivalent', 25, '#91a9ba')
    text(d, (1220, 695), 'FREE 6 DOF', 38, bold=True)
    text(d, (1220, 757), 'Initial speed 8 kn', 29, '#a8c4d4')
    text(d, (1220, 805), 'Assumed thrust', 29, '#a8c4d4')
    text(d, (60, 1015), '1:10 simulation  /  9.99 M particles  /  RTX 5080', 25, '#738d9e')

    comparison = Image.new('RGB', SIZE, BG)
    d = ImageDraw.Draw(comparison)
    text(d, (60, 35), 'WHY SPH HERE?', 51, bold=True)
    text(d, (720, 56), 'Breaking waves + large body motion', 31, '#a8c4d4')
    for i, (kind, title, caption) in enumerate([
        ('fvm', 'BODY-FITTED FVM', 'Conforming mesh + VOF'),
        ('ibm', 'FVM + IBM', 'Fixed grid + immersed hull + VOF'),
        ('sph', 'SPH', 'Particles follow the water')]):
        y = 155 + i * 283
        text(d, (70, y + 65), title, 40, '#64cbe7' if kind == 'sph' else '#e6eef3', True)
        text(d, (70, y + 127), caption, 29, '#a8c4d4')
        comparison.paste(diagram(kind).resize((1014, 260), Image.Resampling.LANCZOS), (820, y))
    text(d, (60, 1030), 'IBM: solid boundary  /  VOF: water-air interface  /  Conceptual comparison', 23, '#738d9e')

    method = Image.new('RGB', SIZE, BG)
    d = ImageDraw.Draw(method)
    text(d, (60, 35), 'SPH', 60, bold=True)
    text(d, (220, 62), 'PARTICLES TO SHIP MOTION', 33, '#a8c4d4')
    text(d, (60, 200), 'INITIAL PARTICLE SLICE', 28, '#91a9ba', True)
    particles = Image.open(OUT / 'gunnerus_stronger_method.png').crop((0, 500, 1080, 980))
    method.paste(particles.resize((1242, 552), Image.Resampling.LANCZOS), (0, 265))
    text(d, (60, 845), 'WATER', 30, '#64cbe7', True)
    text(d, (330, 845), 'HULL + DECKHOUSE', 30, '#e3b75e', True)
    text(d, (1300, 265), '9.99 M PARTICLES', 39, '#64cbe7', True)
    text(d, (1300, 345), 'NEIGHBOR FORCES', 34, bold=True)
    text(d, (1300, 401), 'Fluid loads + gravity', 29, '#a8c4d4')
    story.arrow(d, (1535, 490), (1535, 590), '#64cbe7', width=7)
    text(d, (1300, 660), 'FREE 6 DOF', 44, bold=True)
    text(d, (1300, 727), 'Translation + rotation', 29, '#a8c4d4')
    text(d, (60, 1015), 'DualSPHysics  /  GPU SPH  /  mDBC hull', 25, '#738d9e')
    for title, image in [('setup', setup), ('comparison', comparison), ('method', method)]:
        image.save(OUT / f'gunnerus_landscape_{title}.png')
    return setup, comparison, method


class LandscapeResult(story.Result):
    def __init__(self):
        super().__init__(show_dilute_water=True, surface_level=.38)

    def chart(self, draw, xy, card, t):
        x, y = xy
        title, unit, key, col, factor = card
        v = (self.data[key] if col is None else self.data[key][:, col]) * factor
        tt = self.data['time_s']; use = (tt >= 4) & (tt <= 16.001)
        tt, v = tt[use], v[use]
        low, high = float(v.min()), float(v.max())
        pad = max((high - low) * .1, .01); low -= pad; high += pad
        py = lambda q: y + 135 - (q - low) / (high - low) * 55
        text(draw, (x, y), title, 22, '#b9d0de', True)
        text(draw, (x + 310, y), unit, 21, '#829dad')
        text(draw, (x, y + 28), f'{self.value(key, t, col) * factor:+.2f}', 30, '#64cbe7', True)
        if low <= 0 <= high:
            draw.line((x, py(0), x + 390, py(0)), fill='#2b4051', width=1)
        draw.line((x, y + 135, x + 390, y + 135), fill='#2b4051', width=1)
        keep = tt <= t
        points = list(zip((x + (tt[keep] - 4) / 12 * 390).tolist(), py(v[keep]).tolist()))
        if len(points) > 1:
            draw.line(points, fill='#64cbe7', width=2)
        cursor = x + (t - 4) / 12 * 390
        draw.line((cursor, y + 78, cursor, y + 135), fill='#8fa5b4', width=1)
        text(draw, (x, y + 139), '4 s', 16, '#6c8494')
        text(draw, (x + 354, y + 139), '16 s', 16, '#6c8494')

    def frame(self, t):
        force = self.value('net_vertical_N', t)
        cg = np.array([self.value('center', t, i) for i in range(3)])
        length = float(np.clip(np.sqrt(abs(force) / 1000), .2, 1.3))
        sign = 1 if force >= 0 else -1
        self.arrow.mapper.SetInputData(pv.Arrow(start=(cg[0], -1., .7), direction=(0, 0, sign), scale=length))
        self.arrow.GetProperty().SetColor(*(np.array([100, 203, 231]) if force >= 0 else np.array([239, 162, 94])) / 255)
        self.scene.pl.window_size = (1320, 670)
        close = self.scene.render(t)
        side, front = self.auxiliary_views(t)
        image = Image.new('RGB', SIZE, BG)
        image.paste(close, (0, 130))
        image.paste(side, (1380, 155)); image.paste(front, (1380, 435))
        d = ImageDraw.Draw(image)
        text(d, (50, 28), 'GUNNERUS', 47, bold=True)
        text(d, (50, 90), 'FOCUSED HEAD WAVES', 26, '#64cbe7', True)
        text(d, (1010, 38), f't {t:05.2f} s', 29)
        text(d, (1480, 38), '1:10 MODEL', 28, '#9ab0c0')
        text(d, (1380, 120), 'SIDE', 24, '#b9d0de', True)
        text(d, (1380, 400), 'FRONT', 24, '#b9d0de', True)
        text(d, (1380, 707), 'FLUID-SIDE PRESSURE', 23, '#b9d0de', True)
        colors = (colormaps['turbo'](np.linspace(0, 1, 270))[:, :3] * 255).astype('uint8')
        image.paste(Image.fromarray(np.repeat(colors[None, :, :], 14, axis=0)), (1475, 755))
        low, high = self.pressure_limits
        text(d, (1380, 745), f'≤{low:g}', 24, '#b9d0de')
        text(d, (1770, 745), f'{high:g}+ kPa', 24, '#b9d0de')
        text(d, (50, 806), f'U {self.value("speed_body_m_s", t):+.2f} m/s', 27)
        text(d, (650, 806), f'Heave {self.value("heave_m", t) * 1000:+.0f} mm', 27)
        for x, card in zip([50, 525, 1000, 1475], self.cards):
            self.chart(d, (x, 866), card, t)
        text(d, (50, 1042), 'Model units  /  Pitch: bow-down positive  /  Assumed thrust', 20, '#738d9e')
        text(d, (1220, 1042), 'Visual model: IHB NTNU / vessel.js | MIT', 20, '#627f92')
        return image


def main(preview_only):
    intro = intro_frames()
    result = LandscapeResult()
    if preview_only:
        for t in [8.34, 12.3]:
            result.frame(t).save(OUT / f'gunnerus_landscape_result_{t:g}.png')
        return
    target = OUT / f'{NAME}.mp4'
    proc = subprocess.Popen(['ffmpeg', '-v', 'error', '-y', '-f', 'rawvideo', '-pix_fmt', 'rgb24', '-s', '1920x1080', '-r', str(FPS),
                             '-i', '-', '-an', '-c:v', 'libx264', '-preset', 'fast', '-crf', '19', '-pix_fmt', 'yuv420p',
                             '-movflags', '+faststart', str(target)], stdin=subprocess.PIPE)
    try:
        for frame, seconds in zip(intro, [2, 4, 2]):
            for _ in range(seconds * FPS):
                proc.stdin.write(np.asarray(frame).tobytes())
        for i, t in enumerate(np.linspace(4., 16., 18 * FPS)):
            proc.stdin.write(np.asarray(result.frame(float(t))).tobytes())
            if i % 60 == 0:
                print(f'Landscape result frames {i}/540', flush=True)
    finally:
        proc.stdin.close()
    assert proc.wait() == 0
    check = subprocess.run(['ffmpeg', '-v', 'error', '-nostats', '-progress', 'pipe:1', '-i', str(target), '-f', 'null', '-'], capture_output=True, text=True, check=True)
    count = [int(n) for n in re.findall(r'^frame=(\d+)', check.stdout, re.M)][-1]
    assert count == 780
    metadata = json.loads((OUT / 'gunnerus_stronger_waves_methods_portrait.json').read_text())
    metadata.update(dimensions=list(SIZE), sha256=hashlib.sha256(target.read_bytes()).hexdigest(), source_video=None,
                    rendering='Re-rendered from the same actual cached particle IDs, rigid poses and native monitors; no new simulation.',
                    layout='Large perspective left, side and front pressure right, four synchronized native diagnostic traces below.',
                    decoded_frames=count, final_visual_review='pending')
    target.with_suffix('.json').write_text(json.dumps(metadata, indent=2))
    print(f'Complete: {target}', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--preview-only', action='store_true')
    main(parser.parse_args().preview_only)
