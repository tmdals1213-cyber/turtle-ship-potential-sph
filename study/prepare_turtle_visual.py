# SPDX-License-Identifier: GPL-3.0-or-later
from project_paths import ASSET_DIR, SOLVER_BIN, FONT_DIR
"""Prepare a locally licensed visual skin and a separately declared simple hull.

The downloaded source and all converted meshes stay in work/, outside the
distributable project. Run export_turtle_blender.py first.
"""
from pathlib import Path
import json,sys,hashlib
import numpy as np
import trimesh
from scipy.spatial import ConvexHull
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from common import ROOT,config,physics
p=ASSET_DIR
a=np.load(p/'turtle_source_mesh.npz');v=a['vertices'].copy();f=a['faces'];ids=a['object_ids'];names=a['names']
# Uniform scale; chosen body length, waterline and mass are visualization/demo
# assumptions, not historical ship measurements or a reconstruction validation.
scale=3/14.6;origin=np.array([-1.70,.0198,-4.90]);v=(v-origin)*scale
m=trimesh.Trimesh(v,f,process=False)
m.export(p/'turtle_visual.stl')
# The closed topside overlaps suppress any hidden water sheet inside the hull.
# Cut a waterline opening from the outer hull intersections (not the oars).
hull_ids=[i for i,n in enumerate(names) if n in ['Cube.000','Cube.003','Cube.004','Cube.006','Cube.008','Cube.009']]
tri=v[f[np.isin(ids,hull_ids)]];cross=[]
for j in range(3):
    aa,bb=tri[:,j],tri[:,(j+1)%3];ok=aa[:,2]*bb[:,2]<0
    q=aa[ok]+(bb[ok]-aa[ok])*(-aa[ok,2]/(bb[ok,2]-aa[ok,2]))[:,None];cross.extend(q[:,:2])
xy=np.asarray(cross);poly=xy[ConvexHull(xy).vertices]
# Simple muted material groups. No invented historical material claims.
centres=v[f].mean(1);colors=np.tile([151,106,65],(len(f),1)).astype(np.uint8)
colors[centres[:,2]<.16]=[104,72,45]
colors[centres[:,2]>.55]=[186,155,103]
colors[ids==list(names).index('Object01')]=[181,146,86]
np.savez_compressed(p/'turtle_visual.npz',vertices=v,faces=f,colors=colors,waterline_polygon=poly)
c=config();c['hull'].update(length_m=3,beam_m=1.0,draft_m=.28,freeboard_m=.35,center_of_mass_m=[0,0,.08])
c['description']='Illustrative broad Wigley proxy carrying a turtle-ship visual skin; not turtle-ship hydrodynamic prediction.'
(ROOT/'case_turtle_proxy.json').write_text(json.dumps(c,indent=2))
meta=dict(source_page='https://wifewell.tistory.com/entry/거북선-3D파일FBXOBJ-및-블렌더-무료-나눔',author='예은이네 / 와이프잘모시기',source_archive_sha256=(hashlib.sha256((p/'turtle_blend.7z').read_bytes()).hexdigest() if (p/'turtle_blend.7z').exists() else None),
    terms='Author allows rendered educational/commercial video use; model redistribution and sale prohibited. Keep all original and converted geometry local.',
    original_objects=len(names),triangles=len(f),scale=scale,source_origin=origin.tolist(),
    local_mesh='work/turtle-visual-source/turtle_visual.stl',proxy_hydrostatics=physics(c),
    physical_scope='Simple broad Wigley hull; illustrative CG and inertia; skin follows calculated heave and pitch rigidly. No oar, sail, rudder, deck or dragon-head hydrodynamic forces.')
(ROOT/'results/turtle_visual_provenance.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(meta,ensure_ascii=False,indent=2))
