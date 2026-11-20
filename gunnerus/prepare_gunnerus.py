"""Prepare declared Gunnerus flotation / idealized-force startup inputs.

The MIT visual mesh is an approximation, with openings sealed. The external
force is world-fixed at the CG: no resolved propeller or fluid momentum source.
"""
from pathlib import Path
import argparse, hashlib, json, os, xml.etree.ElementTree as ET
import numpy as np
import trimesh
from xml_helpers import el, box

ROOT = Path(__file__).resolve().parent
GEO = ROOT/'geometry/gunnerus'
SCALE, RHO, DRAFT, KG = 10., 1025., 2.7, 4.0

def physical_model():
    GEO.mkdir(exist_ok=True)
    source=ROOT/'geometry/candidates/gunnerus_vesseljs/hull_closed_candidate.stl'
    m=trimesh.load(source,process=True); m.fix_normals()
    v=m.vertices.copy()
    # glTF: bow -Z, starboard +X, up +Y. This is a proper rotation.
    m.vertices=np.c_[-v[:,2],-v[:,0],v[:,1]]
    keel=float(m.bounds[0,2]); shift=-keel-DRAFT
    m.apply_translation([0,0,shift]);m.apply_scale(1/SCALE)
    assert m.is_watertight and m.is_winding_consistent and m.volume>0
    wet=m.slice_plane([0,0,0],[0,0,-1],cap=False)
    assert wet.volume>0
    cb=wet.center_mass;cg=np.array([cb[0],cb[1],(KG-DRAFT)/SCALE])
    section=m.section(plane_origin=[0,0,0],plane_normal=[0,0,1])
    area=ixx=0.
    for pts in section.discrete:
        q=pts[:,:2]-cg[:2];x,y=q.T;xn,yn=np.roll(q,-1,axis=0).T
        cross=x*yn-xn*y
        area+=abs(cross.sum()/2)
        ixx+=abs(((y*y+y*yn+yn*yn)*cross).sum()/12)
    mass=RHO*wet.volume;beam=float(np.ptp(m.bounds[:,1]));lpp=33.90/SCALE
    inertia=mass*np.array([.35*beam,.25*lpp,.25*lpp])**2
    gm=float(cb[2]+ixx/wet.volume-cg[2]);assert gm>0
    hull=GEO/'hull.stl';m.export(hull)
    normal=GEO/'hull_normals.vtk'
    with normal.open('w',encoding='utf-8') as f:
        f.write('# vtk DataFile Version 3.0\nGunnerus closed visual hull\nASCII\nDATASET POLYDATA\n')
        f.write(f'POINTS {len(m.vertices)} float\n');np.savetxt(f,m.vertices,fmt='%.9g')
        f.write(f'POLYGONS {len(m.faces)} {4*len(m.faces)}\n')
        np.savetxt(f,np.c_[np.full(len(m.faces),3),m.faces],fmt='%d')
        f.write(f'CELL_DATA {len(m.faces)}\nSCALARS Mk unsigned_short 1\nLOOKUP_TABLE default\n')
        np.savetxt(f,np.full(len(m.faces),12),fmt='%d')
    p=dict(model='Gunnerus visual-hull demonstration',scale_ratio=SCALE,rho_kg_m3=RHO,
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        hull_file=str(hull),normals_file=str(normal),closed_watertight=True,solver_triangles=len(m.faces),
        mass_kg=float(mass),full_scale_equivalent_mass_t=float(mass*SCALE**3/1000),
        mass_basis='Derived from closed visual-hull displacement at the declared 2.7 m draft; not measured vessel mass.',
        Lpp_m=lpp,beam_reference_m=9.6/SCALE,draft_reference_m=DRAFT/SCALE,
        center_of_mass_m=cg.tolist(),center_of_buoyancy_m=cb.tolist(),inertia_kg_m2=inertia.tolist(),
        assumed_full_scale_KG_m=KG,GM_roll_m=gm,waterplane_area_m2=float(area),
        displacement_m3=float(wet.volume),closed_volume_m3=float(m.volume),bounds_m=m.bounds.tolist(),
        visual_to_model=dict(axes='x=-glTF.z, y=-glTF.x, z=glTF.y',scale=1/SCALE,vertical_shift_before_scale_m=shift),
        limitations=['Approximate visual hull, not production CAD; windows/openings sealed and flooding excluded.',
            'LCG aligned with geometric buoyancy centre for initial zero trim; KG and radii of gyration are declared assumptions.',
            'Propellers, rudders, bow thruster and above-hull superstructure do not act on the fluid in these startup cases.',
            'Froude model scale does not also preserve full-scale Reynolds number.',
            'The deck/bulwark drainage geometry still needs separate qualification before a green-water load claim.'])
    (GEO/'physical_model.json').write_text(json.dumps(p,indent=2),encoding='utf-8')
    return p

def prepare(kind,p,dp=.05,duration=2.):
    name=f'gunnerus_{kind}_{round(dp*1000):03d}_{duration:g}s'
    folder=ROOT/'sph'/name
    if folder.exists():raise FileExistsError(f'Preserving existing case: {folder}')
    folder.mkdir()
    depth=1.25;x0,x1,y0,y1,z1=-7.,7.,-2.,2.,1.
    speed=8*.514444/SCALE**.5 if kind=='thrust' else 0.
    thrust=30. if kind=='thrust' else 0.
    root=ET.Element('case');cd=el(root,'casedef');co=el(cd,'constantsdef')
    el(co,'gravity',x=0,y=0,z=-9.81)
    for tag,value,extra in [('rhop0',RHO,{}),('rhopgradient',2,{}),('hswl',depth,{'auto':'false'}),
        ('gamma',7,{}),('speedsystem',max(speed,1.3),{'auto':'false'}),('coefsound',20,{}),
        ('speedsound',60,{'auto':'false'}),('coefh',1.2,{}),('cflnumber',.2,{})]:el(co,tag,value=value,**extra)
    el(cd,'mkconfig',boundcount=240,fluidcount=9)
    de=el(el(cd,'geometry'),'definition',dp=dp)
    # Half-spacing placement relative to nominal free surface.
    el(de,'pointmin',x=x0-5*dp,y=y0-5*dp,z=-depth-5*dp-dp/2)
    el(de,'pointmax',x=x1+5*dp,y=y1+5*dp,z=z1+dp)
    cmd=el(el(cd.find('geometry'),'commands'),'mainlist');el(cmd,'setshapemode').text='dp | bound'
    el(cmd,'setdrawmode',mode='full');el(cmd,'setmkfluid',mk=0)
    box(cmd,(x0,y0,-depth),(x1-x0,y1-y0,depth-dp/2))
    el(cmd,'setmkbound',mk=0)
    box(cmd,(x0-3*dp,y0-3*dp,-depth-3*dp),(x1-x0+6*dp,y1-y0+6*dp,3*dp))
    for yy in [y0-3*dp,y1]:box(cmd,(x0-3*dp,yy,-depth),(x1-x0+6*dp,3*dp,z1+depth))
    for xx in [x0-3*dp,x1]:box(cmd,(xx,y0,-depth),(3*dp,y1-y0,z1+depth))
    el(cmd,'setshapemode').text='actual | bound';el(cmd,'setmkbound',mk=2)
    stl=el(cmd,'drawfilestl',file=p['hull_file'],autofill='true',advanced='true');el(stl,'depth',depthmin=.001*dp)
    el(el(cmd,'fillpoint',x=0,y=0,z=-.1),'modefill').text='void | fluid'
    body=el(el(cd,'floatings'),'floating',mkbound=2)
    el(body,'massbody',value=p['mass_kg']);el(body,'center',**dict(zip('xyz',p['center_of_mass_m'])))
    el(body,'inertia',**dict(zip('xyz',p['inertia_kg_m2'])))
    el(body,'translationDOF',x=1,y=1,z=1);el(body,'rotationDOF',x=1,y=1,z=1)
    if speed:el(body,'linearvelini',x=speed,y=0,z=0)
    if thrust:
        schedule=el(body,'linearforce')
        for t,f in [(0,0),(.3,0),(1.3,thrust),(duration+1,thrust)]:el(schedule,'force',time=t,x=f,y=0,z=0)
    ng=el(el(cd,'normals',active='true'),'norgeometry')
    el(ng,'geometryfile',file=os.path.relpath(p['normals_file'],folder/'generated').replace('\\','/'))
    el(ng,'distanceh',v=2);el(ng,'maxsizeh',v=2.1)
    pars=el(el(root,'execution'),'parameters')
    params=dict(Boundary=2,SlipMode=1,NoPenetration=0,StepAlgorithm=2,Kernel=2,
        ViscoTreatment=1,Visco=.01,DensityDT=2,DensityDTvalue=.1,Shifting=0,RigidAlgorithm=0,
        FtPause=0,FtSaveAce=1,TimeMax=duration,TimeOut=.05,CoefDtMin=.01,DtIni=0,DtMin=0,
        RhopOutMin=700,RhopOutMax=1300,MinFluidStop=.99)
    for key,value in params.items():el(pars,'parameter',key=key,value=value)
    domain=el(pars,'simulationdomain');el(domain,'posmin',x=x0-5*dp,y=y0-5*dp,z=-depth-5*dp-dp/2)
    el(domain,'posmax',x=x1+5*dp,y=y1+5*dp,z=z1+dp)
    ET.indent(root);ET.ElementTree(root).write(folder/'Case_Def.xml',encoding='utf-8',xml_declaration=True)
    setup=dict(status='prepared_not_run',case=name,kind=kind,dp_m=dp,duration_s=duration,scale_ratio=SCALE,
        domain_bounds_m=[x0,x1,y0,y1,-depth,z1],water_depth_m=depth,physical_model=p,
        initial_velocity_m_s=[speed,0,0],all_six_dofs_free=True,
        propulsion=dict(mode='world-fixed external force at CG',plateau_force_N=[thrust,0,0],
            full_scale_Froude_force_equivalent_kN=thrust*SCALE**3/1000,ramp_interval_s=[.3,1.3],
            basis='Assumed startup diagnostic force, not obtained from propeller RPM or measured open-water curves.',
            acts_on_fluid=False,resolved_propeller=False,body_following=False),
        gpu_wall_limit_s=None,limitations=p['limitations']+[
            'Startup only: no equilibrium, calibrated resistance or propeller-efficiency claim.',
            'Finite closed calm-water tank; no incident wave and no outlet.',
            'Fixed laboratory force direction is an idealized straight-run approximation; no thrust-line moment.',
            'Native force saved every 0.05 model seconds; this is not sufficient for a resolved slam-peak claim.'])
    (folder/'setup.json').write_text(json.dumps(setup,indent=2),encoding='utf-8');return name

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--kind',choices=['float','thrust','both'],default='both')
    a=ap.parse_args();p=physical_model()
    cases=[prepare(k,p) for k in (['float','thrust'] if a.kind=='both' else [a.kind])]
    print(json.dumps(dict(cases=cases,mass_kg=p['mass_kg'],GM_roll_m=p['GM_roll_m']),indent=2))
