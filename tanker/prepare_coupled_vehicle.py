"""Prepare native SPH/Chrono rigid-truck surrogate with unilateral road contact.

The non-rotating support frame has a prescribed horizontal acceleration.
Inertial loads act on BOTH water and dry vehicle. Tank pose responds freely to
SPH and contact forces. This is a rigid support-platform surrogate, not a
steering/ABS/tire-handling model of an identified road vehicle.
"""
from pathlib import Path
import argparse,json,os,subprocess
import xml.etree.ElementTree as ET
import numpy as np
import pyvista as pv
import trimesh
from prepare_braking import el,profile,CASES
from mdbc_tanker import configure

ROOT=Path(__file__).resolve().parent
from tanker_paths import BIN


def write_points(path,xyz):
    with path.open('w') as out:
        out.write('# vtk DataFile Version 3.0\nTruck boundary particles\nASCII\nDATASET POLYDATA\n')
        out.write(f'POINTS {len(xyz)} float\n');np.savetxt(out,xyz,fmt='%.8f')
        out.write(f'VERTICES {len(xyz)} {2*len(xyz)}\n')
        np.savetxt(out,np.column_stack((np.ones(len(xyz),int),np.arange(len(xyz)))),fmt='%d')


def box(bounds):
    lo=np.array(bounds)[::2];hi=np.array(bounds)[1::2]
    return trimesh.creation.box(extents=hi-lo,transform=trimesh.transformations.translation_matrix((hi+lo)/2))


def collisions(folder):
    meshes=[box([-3.6,4.7,-.8,.8,-1.35,-1.12]),box([3.25,4.65,-.86,.86,-1.12,.65])]
    tank=trimesh.creation.cylinder(radius=1.02,height=6.04,sections=64)
    tank.apply_transform(trimesh.transformations.rotation_matrix(np.pi/2,[0,1,0]));meshes.append(tank)
    # Fixed wheels on a rigid chassis: no wheel spin or suspension is claimed.
    for x in (-2.1,-.8,3.8):
        for y in (-1.,1.):
            wheel=trimesh.creation.cylinder(radius=.5,height=.25,sections=32)
            wheel.apply_transform(trimesh.transformations.rotation_matrix(np.pi/2,[1,0,0]))
            wheel.apply_translation((x,y,-1.60));meshes.append(wheel)
    objects={'truck_contact':trimesh.util.concatenate(meshes),'road_contact':box([-32,32,-32,32,-2.5,-2.1])}
    for name,mesh in objects.items():
        mesh.export(folder/f'{name}.obj')
        with (folder/f'{name}.vtk').open('w') as out:
            out.write('# vtk DataFile Version 3.0\nCollision surface\nASCII\nDATASET POLYDATA\n')
            out.write(f'POINTS {len(mesh.vertices)} float\n');np.savetxt(out,mesh.vertices,fmt='%.8f')
            out.write(f'POLYGONS {len(mesh.faces)} {4*len(mesh.faces)}\n');np.savetxt(out,np.column_stack((np.full(len(mesh.faces),3),mesh.faces)),fmt='%d')


def prepare(layout,lat_g,duration=8.,version='v1',source_suffix='brake_040_v1',lateral_start=2.,sound_speed=110.,cfl=.25):
    source=ROOT/'sph'/f'{layout}_{source_suffix}'
    name=f'{layout}_coupled_lat{round(lat_g*1000):03d}_{version}'
    folder=ROOT/'sph'/name
    if folder.exists():
        if (folder/'initial_audit.json').exists():return folder
        if (folder/'output').exists():raise FileExistsError('Preserve previous simulation')
    (folder/'generated').mkdir(parents=True,exist_ok=True)
    collisions(folder)
    boundary=[pv.read(source/'walls.vtk').points]
    if (source/'plates.vtk').exists():boundary.append(pv.read(source/'plates.vtk').points)
    write_points(folder/'body_particles.vtk',np.concatenate(boundary))
    # One bookkeeping particle assigns a fixed-body marker to the contact road.
    write_points(folder/'road_marker.vtk',np.array([[-5.5,0,-2.34]]))
    root=ET.parse(source/'Case_Def.xml').getroot();cd=root.find('casedef');ge=cd.find('geometry');de=ge.find('definition')
    de.find('pointmin').attrib.update(x='-5.7',y='-1.3',z='-2.5');de.find('pointmax').attrib.update(x='4.9',y='1.3',z='1.3')
    main=ge.find('commands/mainlist');main.clear();el(main,'setshapemode').text='dp | bound';el(main,'setdrawmode',mode='full')
    for tag,mk,path in [('setmkfluid',0,source/'fluid.vtk'),('setmkbound',0,folder/'road_marker.vtk'),('setmkbound',2,folder/'body_particles.vtk')]:
        el(main,tag,mk=mk);el(el(main,'drawfilevtk',file=str(path)),'polyselec').text='points'
    prop=el(cd,'properties');link=el(prop,'links');el(link,'link',mkbound=0,property='contact')
    material=el(prop,'property',name='contact')
    for k,v in [('Young_Modulus',2e7),('PoissonRatio',.3),('Kfric',.8),('Sfric',.8),('Restitution_Coefficient',.02)]:el(material,k,value=v)
    floating=el(el(cd,'floatings'),'floating',mkbound=2,property='contact')
    mass=9000.
    el(floating,'massbody',value=mass);el(floating,'center',x=.8,y=0,z=-.9)
    el(floating,'inertia',x=mass*.9**2,y=mass*2.2**2,z=mass*2.3**2)
    el(floating,'translationDOF',x=1,y=1,z=1);el(floating,'rotationDOF',x=1,y=1,z=1)
    # Fluid added linear acceleration is -a_reference; add same physical force at dry CG.
    t,ax,brake=profile();lateral=np.zeros_like(t)
    rise=(t>=lateral_start)&(t<lateral_start+.6);hold=(t>=lateral_start+.6)&(t<4.4);fall=(t>=4.4)&(t<5.)
    lateral[rise]=lat_g*9.81*.5*(1-np.cos(np.pi*(t[rise]-lateral_start)/.6))
    lateral[hold]=lat_g*9.81
    lateral[fall]=lat_g*9.81*.5*(1+np.cos(np.pi*(t[fall]-4.4)/.6))
    ex=root.find('execution');sp=ex.find('special');series=sp.find('accinputs/accinput/acctimes');series.clear()
    forces=el(floating,'linearforce')
    for tm,a,b in zip(t,ax,lateral):
        el(series,'timevalue',time=f'{tm:.9f}',linx=f'{a:.9f}',liny=f'{-b:.9f}')
        el(forces,'force',time=f'{tm:.9f}',x=f'{mass*a:.8f}',y=f'{-mass*b:.8f}',z=0)
    chrono=el(sp,'chrono');el(chrono,'savedata',value=.01);el(chrono,'ompthreads',value=2)
    col=el(chrono,'collision',activate='true');el(col,'distancedp',value=.05);el(col,'contactmethod',value=0)
    el(chrono,'bodyfixed',id='road',mkbound=0,modelfile='../road_contact.vtk')
    el(chrono,'bodyfloating',id='truck',mkbound=2,modelfile='../truck_contact.vtk')
    pars=ex.find('parameters')
    for k,v in [('RigidAlgorithm',3),('FtPause',0),('FtSaveAce',1),('TimeMax',duration),('TimeOut',.05)]:
        node=pars.find(f"parameter[@key='{k}']")
        if node is None:el(pars,'parameter',key=k,value=v)
        else:node.set('value',str(v))
    dom=pars.find('simulationdomain');dom.find('posmin').attrib.update(x='-30',y='-30',z='-3');dom.find('posmax').attrib.update(x='30',y='30',z='8')
    configure(root,folder,layout,floating=True)
    # The current mDBC2 no-penetration force branch explicitly excludes floating
    # neighbors (JSphGpu_ker.cu). Use DBC for the moving tank, and independently
    # qualify compartment isolation; the fixed-tank mDBC test is insufficient.
    pars.find("parameter[@key='Boundary']").set('value','1')
    pars.find("parameter[@key='NoPenetration']").set('value','0')
    cd.remove(cd.find('normals'))
    cd.find('constantsdef/speedsound').set('value',str(sound_speed))
    cd.find('constantsdef/cflnumber').set('value',str(cfl))
    ET.indent(root);ET.ElementTree(root).write(folder/'Case_Def.xml',encoding='utf-8',xml_declaration=True)
    source_setup=json.loads((source/'setup.json').read_text())
    setup={**source_setup,'case':name,'duration_s':duration,'coupled':True,'lateral_g':lat_g,
           'lateral_profile':{'start_s':lateral_start,'ramp_s':.6,'release_start_s':4.4,'end_s':5.},
           'boundary':'DBC for floating boundaries; compartment isolation must be checked','sound_speed_m_s':sound_speed,'cfl_number':cfl,
           'dry_mass_kg':mass,'dry_cg_m':[.8,0,-.9],'dry_inertia_kg_m2':[mass*.9**2,mass*2.2**2,mass*2.3**2],
           'contact_mu':.8,'road_z_m':-2.1,'contact_model':'Chrono NSC with six fixed wheel cylinders, chassis/cab/tank contact meshes',
           'support_frame':'Non-rotating accelerating support-platform surrogate. Effective tire forces arise from friction/contact; no wheel spin, steering geometry or suspension.',
           'coupling':'Native two-way SPH to Chrono: fluid forces alter rigid-body pose; moving tank boundaries feed back on the fluid.',
           'input_resultant_peak_g':float(np.hypot(ax,lateral).max()/9.81),
           'rollover_policy':'Report actual motion under declared assumptions. No keyframed body motion, no required outcome.',
           'frame':'Nonrotating accelerating support frame; rigid tank pose is solved relative to the support surface.',
           'limitations':['Single liquid with atmospheric gas pressure; no air compression or cavitation.',
                          'First spatial resolution; boundary checks do not establish resolution convergence.',
                          'Rigid fixed wheels with frictional contact, no calibrated tires or suspension.',
                          'Same dry mass and inertia in all designs; additional baffle structural mass not separately represented.'],
           'vehicle_parameters_status':'Hypothetical generic vehicle; not calibrated production vehicle data.'}
    (folder/'setup.json').write_text(json.dumps(setup,indent=2))
    np.savetxt(folder/'input_acceleration.csv',np.column_stack((t,-ax,lateral)),delimiter=',',header='time_s,reference_ax_m_s2,reference_ay_m_s2',comments='')
    env=os.environ.copy();env['OMP_NUM_THREADS']='4'
    with (folder/'gencase.log').open('w') as log:
        p=subprocess.run([str(BIN/'GenCase_win64.exe'),'Case_Def',str(folder/'generated/Case'),'-save:all'],cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT)
    if p.returncode:raise RuntimeError(f'GenCase failed: {folder}')
    mesh=pv.read(folder/'generated/Case_All.vtk');counts={str(i):int((mesh['Type']==i).sum()) for i in np.unique(mesh['Type'])}
    assert counts.get('3')==source_setup['fluid_particles'] and counts.get('2')==len(np.concatenate(boundary))
    (folder/'initial_audit.json').write_text(json.dumps({'passed':True,'counts_by_type':counts,'fluid_mass_kg':source_setup['fluid_mass_kg']},indent=2))
    print(json.dumps({'case':name,'counts':counts,'max_resultant_g':setup['input_resultant_peak_g']}),flush=True)
    return folder


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--layout',default='A_unpartitioned',choices=CASES);ap.add_argument('--lateral-g',type=float,default=.55);ap.add_argument('--duration',type=float,default=8);ap.add_argument('--version',default='v1')
    ap.add_argument('--source-suffix',default='brake_040_v1')
    ap.add_argument('--lateral-start',type=float,default=2.)
    ap.add_argument('--sound-speed',type=float,default=110.);ap.add_argument('--cfl',type=float,default=.25)
    a=ap.parse_args();prepare(a.layout,a.lateral_g,a.duration,a.version,a.source_suffix,a.lateral_start,a.sound_speed,a.cfl)
