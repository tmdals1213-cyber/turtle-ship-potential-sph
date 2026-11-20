"""Place the unchanged contact geometry at a declared initial tank station."""
import json
from pathlib import Path
import numpy as np
import trimesh
ROOT=Path(__file__).resolve().parent

def at_station(x):
    base=ROOT/'geometry/gunnerus_impact'; out=ROOT/'geometry'/f'gunnerus_impact_x{x:g}'
    if out.exists(): return json.loads((out/'physical_model.json').read_text())
    p=json.loads((base/'physical_model.json').read_text()); hull=trimesh.load(p['hull_file'],process=True)
    shift=np.array([x,0,0]); hull.apply_translation(shift); out.mkdir()
    file=out/'hull.stl'; hull.export(file); normal=out/'hull_normals.vtk'
    with normal.open('w') as f:
        f.write('# vtk DataFile Version 3.0\nGunnerus positioned contact envelope\nASCII\nDATASET POLYDATA\n')
        f.write(f'POINTS {len(hull.vertices)} float\n'); np.savetxt(f,hull.vertices,fmt='%.9g')
        f.write(f'POLYGONS {len(hull.faces)} {4*len(hull.faces)}\n'); np.savetxt(f,np.c_[np.full(len(hull.faces),3),hull.faces],fmt='%d')
        f.write(f'CELL_DATA {len(hull.faces)}\nSCALARS Mk unsigned_short 1\nLOOKUP_TABLE default\n'); np.savetxt(f,np.full(len(hull.faces),12),fmt='%d')
    p.update(hull_file=str(file.resolve()),normals_file=str(normal.resolve()),initial_station_x_m=x,bounds_m=hull.bounds.tolist())
    for key in ['center_of_mass_m','center_of_buoyancy_m']: p[key]=(np.array(p[key])+shift).tolist()
    p['contact_envelope_bounds_m']=[(np.array(b)+shift).tolist() for b in p['contact_envelope_bounds_m']]
    p['visual_to_model']['translation_model_m']=shift.tolist()
    (out/'physical_model.json').write_text(json.dumps(p,indent=2)); return p
