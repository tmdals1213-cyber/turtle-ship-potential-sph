"""Fluid-facing analytic normal surfaces for tanker mDBC2 boundaries."""
from pathlib import Path
import numpy as np
import trimesh
from prepare_braking import el

ROOT=Path(__file__).resolve().parent


def configure(root,folder,layout,floating=False):
    tank=trimesh.creation.cylinder(radius=1.,height=6.,sections=128)
    tank.apply_transform(trimesh.transformations.rotation_matrix(np.pi/2,[0,1,0]))
    tank.faces=tank.faces[:,::-1]  # Fluid is inside the envelope.
    meshes=[tank];markers=[12 if floating else 10]
    paths=sorted((ROOT/'geometry/concept_v1').glob(f'{layout}_plate_*.stl'))
    for i,path in enumerate(paths):
        m=trimesh.load(path,process=True);m.fix_normals()
        if layout=='D_cross_baffles' and i==2:m.vertices[:,1]*=.08/.06
        else:
            centre=(-1,1)[i];m.vertices[:,0]=centre+(m.vertices[:,0]-centre)*.08/.06
        meshes.append(m);markers.append(12 if floating else 11)
    p=[];faces=[];mk=[]
    for mesh,marker in zip(meshes,markers):
        offset=len(p);p.extend(mesh.vertices.tolist());faces.extend((mesh.faces+offset).tolist());mk.extend([marker]*len(mesh.faces))
    with (folder/'fluid_facing_normals.vtk').open('w') as out:
        out.write('# vtk DataFile Version 3.0\nFluid-facing tanker surfaces\nASCII\nDATASET POLYDATA\n')
        out.write(f'POINTS {len(p)} float\n');np.savetxt(out,p,fmt='%.8f')
        out.write(f'POLYGONS {len(faces)} {4*len(faces)}\n');np.savetxt(out,np.column_stack((np.full(len(faces),3),faces)),fmt='%d')
        out.write(f'CELL_DATA {len(faces)}\nSCALARS Mk unsigned_short 1\nLOOKUP_TABLE default\n');np.savetxt(out,mk,fmt='%d')
    cd=root.find('casedef');existing=cd.find('normals')
    if existing is not None:cd.remove(existing)
    n=el(el(cd,'normals',active='true'),'norgeometry')
    el(n,'geometryfile',file='../fluid_facing_normals.vtk');el(n,'distanceh',v=2);el(n,'maxsizeh',v=2.1)
    pars=root.find('execution/parameters')
    for k,v in [('Boundary',2),('SlipMode',2),('NoPenetration',1)]:
        node=pars.find(f"parameter[@key='{k}']")
        if node is None:el(pars,'parameter',key=k,value=v)
        else:node.set('value',str(v))
    cd.find('constantsdef/speedsound').set('value','110')


def fixed_test(layout='C_separate_compartments'):
    import json,shutil,subprocess,xml.etree.ElementTree as ET
    source=ROOT/'sph'/f'{layout}_brake_040_v1';folder=ROOT/'sph'/f'{layout}_mdbc_check_v2'
    if (folder/'output').exists():raise FileExistsError(folder)
    (folder/'generated').mkdir(parents=True,exist_ok=True)
    root=ET.parse(source/'Case_Def.xml').getroot();configure(root,folder,layout)
    root.find("execution/parameters/parameter[@key='TimeMax']").set('value','6')
    ET.indent(root);ET.ElementTree(root).write(folder/'Case_Def.xml',encoding='utf-8',xml_declaration=True)
    setup=json.loads((source/'setup.json').read_text());setup.update(case=folder.name,duration_s=6,boundary='mDBC2 no-slip/no-penetration',sound_speed_m_s=110)
    (folder/'setup.json').write_text(json.dumps(setup,indent=2));shutil.copy2(source/'braking.csv',folder/'braking.csv')
    shutil.copy2(source/'fluid.vtk',folder/'fluid.vtk')
    from tanker_paths import BIN
    exe=BIN/'GenCase_win64.exe'
    with (folder/'gencase.log').open('w') as out:
        r=subprocess.run([str(exe),'Case_Def',str(folder/'generated/Case'),'-save:all'],cwd=folder,stdout=out,stderr=subprocess.STDOUT)
    if r.returncode:raise RuntimeError('GenCase normals test failed')
    shutil.copy2(source/'initial_audit.json',folder/'initial_audit.json')
    return folder


if __name__=='__main__':print(fixed_test())
