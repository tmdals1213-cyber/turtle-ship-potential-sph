# SPDX-License-Identifier: GPL-3.0-or-later
from project_paths import ASSET_DIR, SOLVER_BIN, FONT_DIR
"""Plausible art direction for the licensed render mesh; no historical material claim."""
from pathlib import Path
import json
import numpy as np
ROOT=Path(__file__).resolve().parent;geo=ASSET_DIR
src=np.load(geo/'turtle_source_mesh.npz');skin=np.load(geo/'turtle_history_visual.npz')
assert np.array_equal(src['faces'],skin['faces'])
v=(src['vertices']-[-1.70,.0198,-4.90])*3/14.6;f=src['faces'];ids=src['object_ids'];names=src['names'];centers=v[f].mean(1)
colors=np.tile([105.,74.,47.],(len(f),1));groups={};palette=np.zeros((len(names),3))
for i,name in enumerate(names):
    mask=ids==i
    if not mask.any():continue
    pts=v[f[mask]];lo=pts.min((0,1));hi=pts.max((0,1));ext=hi-lo;count=int(mask.sum())
    rgb=[105,74,47];material='warm timber'
    if name.startswith('Cube.00'):
        rgb=[100,66,42];material='weathered hull timber'
    if name=='Cube.001':rgb=[99,103,95];material='weathered timber roof'
    if name=='Cube.005':rgb=[129,87,50];material='gunwale timber'
    if name in ['Cube.006','Cube.008']:rgb=[64,47,34];material='dark lower hull'
    if name=='Circle':rgb=[119,103,77];material='cordage / spars'
    if name.startswith('Cone'):
        rgb=[67,75,76];material='dark iron fittings'
        if lo[2]>2:rgb=[145,110,57];material='aged brass finial'
    if name.startswith('Cylinder'):
        rgb=[94,65,39];material='timber spar / oar'
        if count==3260:rgb=[43,57,56];material='dark gunmetal'
    if name=='Object01':rgb=[122,55,38];material='muted red carved bow'
    if i==500:rgb=[133,38,31];material='muted red textile flag'
    colors[mask]=rgb;palette[i]=rgb;groups[str(i)]=dict(source_object=str(name),material=material)
# Restrained, coherent plank/object variation instead of random triangle noise.
plank=np.floor(centers[:,0]/.19)+7*np.floor(centers[:,2]/.075)
variation=1+.055*np.sin(plank*3.17+ids*.83)+.025*np.sin(centers[:,0]*13+centers[:,2]*42)
colors=np.clip(colors*variation[:,None],0,255).astype('u1')
vertex_ids=np.zeros(len(v),dtype=int);vertex_ids[f.ravel()]=np.repeat(ids,3)
grain=1+.025*np.sin(v[:,0]*35+2*np.sin(v[:,2]*7))+.025*np.sin(v[:,2]*17)
vertex_colors=np.clip(palette[vertex_ids]*grain[:,None],0,255).astype('u1')
np.savez_compressed(geo/'turtle_history_materials.npz',colors=colors,vertex_colors=vertex_colors)
(ROOT/'results/turtle_material_direction.json').write_text(json.dumps(dict(
    purpose='Plausible visual materials, not historical evidence. Computational shape, mass and motion are unchanged.',
    palette=['weathered brown timber','gray-brown timber roof','dark gunmetal fittings','muted red bow and flag','natural cordage'],
    source='Same locally licensed turtle mesh as v3; all mesh derivatives remain local-only.',
    material_assignments=groups),indent=2))
print('Prepared',len(colors),'face colors; calculation geometry unchanged.')
