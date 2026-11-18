# SPDX-License-Identifier: GPL-3.0-or-later
"""Run with the source-compatible Blender 2.79, with auto-execution disabled.

Local rendering asset only. The source author permits rendered videos but not
redistribution of the model. Do not bundle source or converted mesh files.
"""
import bpy, json, os, sys
import numpy as np
from mathutils import Vector

out = sys.argv[sys.argv.index('--')+1]
os.makedirs(out, exist_ok=True)
scene = bpy.context.scene
scene.frame_set(1)
records=[]
for obj in scene.objects:
    if obj.type not in ('MESH','CURVE','SURFACE','FONT'): continue
    corners=[list(obj.matrix_world * Vector(p)) for p in obj.bound_box]
    records.append(dict(name=obj.name, type=obj.type, hidden_render=obj.hide_render,
        layers=list(obj.layers), bounds=[np.min(corners,axis=0).tolist(),np.max(corners,axis=0).tolist()],
        vertices=len(obj.data.vertices) if obj.type=='MESH' else 0,
        modifiers=[m.type for m in obj.modifiers]))
with open(os.path.join(out,'objects.json'),'w',encoding='utf-8') as f:json.dump(records,f,indent=2,ensure_ascii=False)
print('OBJECTS',len(records))
for r in records:print(r['name'],r['vertices'],r['bounds'],r['hidden_render'])
if '--export' in sys.argv:
    vertices=[];faces=[];object_ids=[];names=[]
    for obj in scene.objects:
        if obj.type not in ('MESH','CURVE','SURFACE','FONT') or obj.hide_render: continue
        if not any(a and b for a,b in zip(obj.layers,scene.layers)):continue
        mesh=obj.to_mesh(scene,True,'PREVIEW')
        mesh.calc_tessface()
        off=len(vertices);oi=len(names);names.append(obj.name)
        vertices.extend([list(obj.matrix_world*v.co) for v in mesh.vertices])
        for f in mesh.tessfaces:
            ids=list(f.vertices)
            for j in range(1,len(ids)-1):faces.append([off+ids[0],off+ids[j],off+ids[j+1]]);object_ids.append(oi)
        bpy.data.meshes.remove(mesh)
    np.savez_compressed(os.path.join(out,'turtle_source_mesh.npz'),vertices=vertices,faces=faces,object_ids=object_ids,names=names)
    print('EXPORTED',len(vertices),len(faces),len(names))
