"""Sample nearby water pressure with the supplied, unmodified MeasureTool.

This is a fluid-side pressure diagnostic, not wall traction or a new force
integration. Native rigid-body forces remain the source for load histories.
"""
from pathlib import Path
import argparse,csv,json,subprocess,os
import numpy as np
import pyvista as pv
from sph_segments import part_rows,output_dirs,vtk_path

ROOT=Path(__file__).resolve().parent
TOOL=Path(os.environ.get('DUALSPHYSICS_BIN', 'external/DualSPHysics/bin/windows'))/'MeasureTool_win64.exe'


def prepare_pressure(name,first=40,last=None):
    folder=ROOT/'sph'/name; out=folder/'surface_pressure'; out.mkdir(exist_ok=True)
    setup=json.loads((folder/'setup.json').read_text()); dp=setup['dp_m']
    rows=part_rows(folder)
    times={int(r['Part']):float(r['TimeStep [s]']) for r in rows}
    last=max(times) if last is None else last
    assert first in times and last in times
    target=out/f'pressure_{first:04d}_{last:04d}.npz'
    if target.exists(): return target
    mesh=pv.read(setup['physical_model']['hull_file']).triangulate().clean()
    mesh=mesh.compute_normals(cell_normals=True,point_normals=False,consistent_normals=True,auto_orient_normals=True,split_vertices=False)
    assert mesh.n_open_edges==0
    samples=mesh.cell_centers().points+mesh.cell_data['Normals']*(.5*dp)
    mesh.save(out/'contact_mesh.vtp')
    points=out/'sample_points.txt'
    points.write_text(''.join('POINTS\n'+ ' '.join(f'{x:.10g}' for x in q)+'\n' for q in samples))
    records={int(r['Part']):r for r in rows};initial_body=None
    def world_points(part):
        nonlocal initial_body
        cache=ROOT/f'render_cache/gunnerus_comparison/{name}/frame_{part:04d}.npz'
        if cache.exists():
            with np.load(cache) as data:rot=data['R'];trans=data['translation']
        else:
            if initial_body is None:
                a=pv.read(vtk_path(records[0]));mask=a['Type']==2;order=np.argsort(a['Idp'][mask])
                initial_body=(a.points[mask][order].astype(float),a['Idp'][mask][order])
            before,ids=initial_body;b=pv.read(vtk_path(records[part]));mask=b['Type']==2;order=np.argsort(b['Idp'][mask])
            assert np.array_equal(ids,b['Idp'][mask][order]);after=b.points[mask][order].astype(float)
            u,_,vt=np.linalg.svd((before-before.mean(0)).T@(after-after.mean(0)));rot=vt.T@u.T
            assert np.linalg.det(rot)>0;trans=after.mean(0)-rot@before.mean(0)
            assert np.sqrt(np.mean((before@rot.T+trans-after)**2))<dp*.02
        return samples@rot.T+trans
    with (out/'measure.log').open('w') as log:
        for segment,directory in enumerate(output_dirs(folder)):
            selected=[int(r['Part']) for r in rows if r['output_directory']==str(directory) and first<=int(r['Part'])<=last]
            if not selected:continue
            # MeasureTool's trackmk assumes a zero-based part sequence. For a
            # resumed segment, supply the measured world positions explicitly.
            batches=[[p] for p in selected] if segment else [selected]
            for batch in batches:
                current_points=points;tracking=['-trackmk:12']
                if segment:
                    current_points=out/f'points_{batch[0]:04d}.txt';tracking=[]
                    current_points.write_text(''.join('POINTS\n'+' '.join(f'{x:.10g}' for x in q)+'\n' for q in world_points(batch[0])))
                command=[str(TOOL),'-dirdata',str(directory/'data'),'-filexml','AUTO',f'-first:{min(batch)}',f'-last:{max(batch)}',
                    '-points',str(current_points),*tracking,'-onlytype:-all,fluid','-vars:-all,press,kcorr',
                    '-kclimit:0.5','-kcdummy:-999999','-savevtk',str(out/'sample'),'-threads:2']
                subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    pressure=[]; support=[]; stamps=[]
    for part in range(first,last+1):
        d=pv.read(out/f'sample_{part:04d}.vtk'); assert d.n_points==mesh.n_cells
        p=d['Press'].astype('f4'); kc=d['Kcorr'].astype('f4')
        valid=(kc>=.5)&(p>-999998)&np.isfinite(p)
        p[~valid]=np.nan; pressure.append(p/1000); support.append(kc); stamps.append(times[part])
    values=np.array(pressure); assert np.isfinite(values).any()
    np.savez_compressed(target,time_s=stamps,pressure_kPa=values,kernel_sum=support)
    meta=dict(case=name,frames=[first,last],mesh_cells=mesh.n_cells,offset_outward_m=.5*dp,
        source='DualSPHysics MeasureTool, fluid particles only, pressure interpolation with kernel correction',
        point_motion='Original segment: native trackmk:12. Resumed segments: explicit world points from rigid fit to actual hull particle IDs.',min_kernel_sum=.5,
        unsupported='NaN, displayed gray; insufficient water support is not a measured zero pressure.',
        units='Model-scale gauge pressure, kPa',pressure_range_kPa=[float(np.nanmin(values)),float(np.nanmax(values))],
        limitations=['Pressure sampled half a particle spacing outside each contact-face center; not wall traction.',
                    '0.1 s snapshots do not resolve peak slam pressure. No pressure-derived force integration.',
                    'Native solver integrated fluid forces and moments are displayed separately.'])
    target.with_suffix('.json').write_text(json.dumps(meta,indent=2))
    print(json.dumps(meta),flush=True)
    return target


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('case');ap.add_argument('--first',type=int,default=40);ap.add_argument('--last',type=int)
    a=ap.parse_args();prepare_pressure(a.case,a.first,a.last)
