# SPDX-License-Identifier: GPL-3.0-or-later
"""Generate common hull, hydrostatics and Airy-wave reference. No CFD claim."""
from pathlib import Path
import sys,json
import numpy as np
import trimesh
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from common import ROOT,config,physics,wave,hull_quads,airy

def triangulated(v,f):
    tri=np.vstack([f[:,[0,1,2]],f[:,[0,2,3]]])
    # Bow/keel parameter corners collapse onto the centreplane. Do not retain
    # the two coincident internal triangles produced by the split quadrilateral.
    tri=tri[~np.all(np.abs(v[tri,1])<1e-12,axis=1)]
    mesh=trimesh.Trimesh(vertices=v,faces=tri,process=True)
    mesh.update_faces(mesh.nondegenerate_faces());mesh.remove_unreferenced_vertices();mesh.fix_normals()
    return mesh

if __name__=='__main__':
    c=config();out=ROOT/'results';geo=ROOT/'geometry'
    out.mkdir(exist_ok=True);geo.mkdir(exist_ok=True)
    v,f=hull_quads(160,50,above=True);mesh=triangulated(v,f)
    assert mesh.is_watertight and mesh.volume>0
    mesh.export(geo/'wigley_closed.stl')
    v,f=hull_quads(160,50);wet=triangulated(v,f)
    wet.export(geo/'wigley_wetted.stl')
    p=physics(c);w=wave(c=c)
    # Signed volume is well defined here since the missing cap lies at z=0.
    displacement=abs(wet.volume);error=abs(displacement-p['volume_m3'])/p['volume_m3']
    assert error<.001,error
    stats=dict(hull_watertight=bool(mesh.is_watertight),hull_triangles=len(mesh.faces),
        analytic=p,baseline_wave=w,stl_displacement_m3=displacement,relative_displacement_error=error,
        fluid_particle_estimates={str(dp):12*1.5*1.0/dp**3 for dp in [.06,.03,.015,.01]},
        limits=['Analytic surrogate hull; not a manufacturer CAD model.',
                'Mass distribution and CG are illustrative declared assumptions.',
                'Finite-depth wave-flume case, not a full-scale open-ocean prediction.',
                'Theory animation is an analytic incident-wave reference; no ship response is invented.'])
    (out/'setup_summary.json').write_text(json.dumps(stats,indent=2))
    rows=[wave(r,c=c) for r in c['waves']['wavelength_over_length_sweep']]
    (out/'wave_matrix.json').write_text(json.dumps(rows,indent=2))
    t=np.linspace(0,8*w['encounter_period_s'],801)
    eta,u,vz=airy(0,0,t,w,c)
    np.savetxt(out/'theory_at_ship_origin.csv',np.c_[t,eta,u,vz],delimiter=',',
        header='time_s,eta_m,u_ship_frame_m_s,w_m_s',comments='')
    print(json.dumps(stats,indent=2))
