"""Add declared convex deckhouse envelopes to the sealed main hull.

Original submerged geometry and assumed vessel mass/CG/inertia are retained.
The new contact envelope is not detailed drainage, flooding or structural CAD.
"""
import hashlib,json
from pathlib import Path
import numpy as np
import trimesh

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'geometry/gunnerus_impact'


def main():
    if OUT.exists(): raise FileExistsError('Preserve existing impact geometry')
    OUT.mkdir()
    p=json.loads((ROOT/'geometry/gunnerus/physical_model.json').read_text())
    original=trimesh.load(p['hull_file'],process=True)
    scene=trimesh.load(ROOT/'geometry/candidates/gunnerus_vesseljs/Gunnerus.glb',force='scene',process=False)
    vertices={}
    for node in scene.graph.nodes_geometry:
        matrix,key=scene.graph[node]; v=trimesh.transform_points(scene.geometry[key].vertices,matrix)
        vertices[node]=np.c_[-v[:,2],-v[:,0],v[:,1]+p['visual_to_model']['vertical_shift_before_scale_m']]/10
    lower=np.vstack([v for name,v in vertices.items() if name in ['WhitePaint','VesselGlass','BluePaint','bridgeGlass']])
    lower=lower[(lower[:,0]>=-.38)&(lower[:,0]<=1.3)&(abs(lower[:,1])<=.37)&(lower[:,2]>=.33)&(lower[:,2]<=.705)]
    upper=np.vstack([v for name,v in vertices.items() if name in ['bridgeGlass','bridgeBlue','bridgeInside'] or name.startswith('vesselShared_bridge_')])
    upper=upper[(upper[:,0]>=-.36)&(upper[:,0]<=.47)&(abs(upper[:,1])<=.37)&(upper[:,2]>=.66)&(upper[:,2]<=.97)]
    envelopes=[trimesh.convex.convex_hull(v) for v in [lower,upper]]
    for name,mesh in zip(['lower_deckhouse','wheelhouse'],envelopes): mesh.export(OUT/f'{name}.stl')
    hull=trimesh.boolean.union([original]+envelopes,engine='manifold',check_volume=True)
    assert hull.is_watertight and hull.is_winding_consistent and hull.volume>original.volume
    assert len(hull.split())==1, 'Deckhouses must connect to the main hull'
    old_wet=original.slice_plane([0,0,0],[0,0,-1],cap=False)
    new_wet=hull.slice_plane([0,0,0],[0,0,-1],cap=False)
    rel=abs(new_wet.volume-old_wet.volume)/old_wet.volume; assert rel<1e-6
    file=OUT/'hull.stl'; hull.export(file); normal=OUT/'hull_normals.vtk'
    with normal.open('w') as f:
        f.write('# vtk DataFile Version 3.0\nGunnerus hull plus deckhouse contact envelopes\nASCII\nDATASET POLYDATA\n')
        f.write(f'POINTS {len(hull.vertices)} float\n'); np.savetxt(f,hull.vertices,fmt='%.9g')
        f.write(f'POLYGONS {len(hull.faces)} {4*len(hull.faces)}\n'); np.savetxt(f,np.c_[np.full(len(hull.faces),3),hull.faces],fmt='%d')
        f.write(f'CELL_DATA {len(hull.faces)}\nSCALARS Mk unsigned_short 1\nLOOKUP_TABLE default\n'); np.savetxt(f,np.full(len(hull.faces),12),fmt='%d')
    p.update(model='Gunnerus sealed main hull plus simplified deckhouse contact envelopes',hull_file=str(file.resolve()),
        normals_file=str(normal.resolve()),solver_triangles=len(hull.faces),bounds_m=hull.bounds.tolist(),closed_volume_m3=float(hull.volume),
        contact_geometry='Manifold union of existing main hull and two convex envelopes fitted to MIT visual deckhouse vertices.',
        contact_envelope_bounds_m=[x.bounds.tolist() for x in envelopes],submerged_volume_relative_change=float(rel),
        previous_main_hull_sha256=hashlib.sha256((ROOT/'geometry/gunnerus/hull.stl').read_bytes()).hexdigest(),
        mass_geometry_policy='Mass, CG and inertia are the unchanged assumed loading condition; adding upper contact geometry does not add a second structural mass.')
    p['limitations']=['Approximate sealed hull and convex deckhouse envelopes, not production CAD or internal flooding.',
        'Masts, rails, cranes and propulsion appendages remain visual-only.',
        'Drainage/freeing ports and local panel flexibility are not resolved; deck water and impact loads are qualitative diagnostics.',
        'Mass, KG and radii of gyration remain declared assumptions from the prior pilot.',
        'Froude scale does not preserve full-scale Reynolds number.']
    (OUT/'physical_model.json').write_text(json.dumps(p,indent=2))
    audit=dict(watertight=True,connected_components=1,triangles=len(hull.faces),submerged_volume_relative_change=float(rel),
        added_above_water_volume_m3=float(hull.volume-original.volume),bounds_m=hull.bounds.tolist(),mass_kg=p['mass_kg'])
    (OUT/'geometry_audit.json').write_text(json.dumps(audit,indent=2)); print(json.dumps(audit,indent=2))


if __name__=='__main__': main()
