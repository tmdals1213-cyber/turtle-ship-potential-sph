# SPDX-License-Identifier: GPL-3.0-or-later
from project_paths import ASSET_DIR, SOLVER_BIN, FONT_DIR
"""Cache measured SPH water on a fixed world crop for the portrait film.

Only display reconstruction changes. There is no imposed ship translation and
no analytical wave filling outside the measured data.
"""
from pathlib import Path
import json
import numpy as np
import pyvista as pv
from scipy.ndimage import gaussian_filter

ROOT = Path(__file__).resolve().parent
GEO = ASSET_DIR
cache = GEO / 'free6_history_40_portrait_fixed'
cache.mkdir(exist_ok=True)
source = GEO / 'free6_history_40_render'
dx = .04
origin = np.array([-2.4, 0., 0.])
axes = [np.arange(-3.2, 3.2 + dx/2, dx), np.arange(-1.2, 1.2 + dx/2, dx), np.arange(-.4, .72 + dx/2, dx)]
extended = [np.r_[v[0]-dx*np.arange(3, 0, -1), v, v[-1]+dx*np.arange(1, 4)] for v in axes]
edges = [np.r_[v-dx/2, v[-1]+dx/2] for v in extended]
for path in sorted(source.glob('frame_*.npz')):
    target = cache/path.name
    if target.exists():
        continue
    previous = np.load(path)
    part = int(previous['part'])
    mesh = pv.read(ROOT/f'sph/free6_history_40/output/vtks/PartVtk_{part:04d}.vtk')
    points = mesh.points[mesh['Type'] == 3].astype(float) - origin
    hist = np.histogramdd(points, bins=edges)[0]
    occ = gaussian_filter(hist.astype('f4'), sigma=1.1, mode='constant')[3:-3, 3:-3, 3:-3]
    np.savez_compressed(target, occupancy=occ.astype('f2'), x=axes[0], y=axes[1], z=axes[2],
        origin=origin, R=previous['R'], translation=previous['translation'], cg=previous['cg'],
        time=previous['time'], part=part, visualization_version=1)
    if part % 20 == 0:
        print(f'Fixed-world surface {part}', flush=True)
report = dict(source_case='free6_history_40', camera_frame='fixed laboratory',
    fixed_world_origin_m=origin.tolist(), crop_world_x_m=[-5.6,.8],
    surface='Gaussian-smoothed occupancy from actual saved water particles',
    motion='Original rigid poses; no artificial forward translation',
    cached_frames=len(list(cache.glob('frame_*.npz'))))
(ROOT/'results/portrait_fixed_scene.json').write_text(json.dumps(report, indent=2))
print(json.dumps(report))
