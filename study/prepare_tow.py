# SPDX-License-Identifier: GPL-3.0-or-later
from project_paths import ASSET_DIR, SOLVER_BIN, FONT_DIR
"""Prepare a piston tank. --history --free6 selects the published turtle case.
Other switches retain earlier proxy/testing variants; see the reference XML.
"""
import argparse,json,sys,os,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
import trimesh
from prepare_sph import el,box
from prepare import triangulated
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from common import ROOT,physics,wave,hull_quads
ap=argparse.ArgumentParser();ap.add_argument('--name',default='tow_proxy_30');ap.add_argument('--dp',type=float,default=.03);ap.add_argument('--duration',type=float,default=18);ap.add_argument('--visco',type=float,default=.01)
ap.add_argument('--fitted',action='store_true');ap.add_argument('--height',type=float,default=.06);ap.add_argument('--coefh',type=float,default=1);ap.add_argument('--slip-mode',type=int,default=2,choices=[1,2,3])
ap.add_argument('--history',action='store_true');ap.add_argument('--free6',action='store_true');ap.add_argument('--mixed-walls',action='store_true');ap.add_argument('--initial-speed',type=float,default=.5)
a=ap.parse_args();c=json.loads((ROOT/('case_turtle_fitted.json' if a.fitted else 'case_turtle_proxy.json')).read_text(encoding='utf-8'));p=json.loads((ROOT/'results/turtle_fitted_hydrostatics.json').read_text(encoding='utf-8')) if a.fitted else physics(c);w=wave(height=a.height,c=c);dp=a.dp
if a.history:
    a.fitted=True;c=json.loads((ROOT/'case_turtle_history.json').read_text(encoding='utf-8'));p=json.loads((ROOT/'results/turtle_history_model.json').read_text(encoding='utf-8'));w=wave(height=a.height,c=c)
out=ROOT/'sph'/a.name;out.mkdir(parents=True,exist_ok=True);offset=np.array([-3.,0.,0.]);cg=np.array(c['hull']['center_of_mass_m'])+offset
if a.fitted:
    hm=trimesh.load(ASSET_DIR/('turtle_history_closed.stl' if a.history else 'turtle_fitted_closed.stl'),process=True);hm.apply_translation(offset)
    hullpath=ASSET_DIR/f'{a.name}_hull.stl'
else:
    v,f=hull_quads(160,48,above=True,c=c);hm=triangulated(v+offset,f);hullpath=out/'hull.stl'
assert hm.is_watertight;hm.export(hullpath)
root=ET.Element('case');cd=el(root,'casedef');con=el(cd,'constantsdef');el(con,'gravity',x=0,y=0,z=-9.81)
for tag,val,at in [('rhop0',c['environment']['rho_kg_m3'],{}),('rhopgradient',2,{}),('hswl',1,{'auto':'false'}),('gamma',7,{}),('speedsystem',.5,{'auto':'false'}),('coefsound',20,{}),('speedsound',50,{'auto':'false'}),('coefh',a.coefh,{}),('cflnumber',.2,{})]:el(con,tag,value=val,**at)
el(cd,'mkconfig',boundcount=240,fluidcount=9)
x0,x1,y0,y1,z0,z1=-9.,9.,-1.2,1.2,-1.,.7
geom=el(cd,'geometry');de=el(geom,'definition',dp=dp);el(de,'pointmin',x=x0-4*dp,y=y0-4*dp,z=z0-4*dp);el(de,'pointmax',x=x1+.3,y=y1+4*dp,z=z1)
main=el(el(geom,'commands'),'mainlist');el(main,'setshapemode').text='dp | bound' if a.mixed_walls else 'actual | bound';el(main,'setdrawmode',mode='full')
el(main,'setmkfluid',mk=0);box(main,(x0,y0,z0),(x1-x0,y1-y0,-z0))
el(main,'setmkbound',mk=0)
box(main,(x0-3*dp,y0-3*dp,z0-3*dp),(x1-x0+6*dp,y1-y0+6*dp,3*dp))
box(main,(x0-3*dp,y0-3*dp,z0),(x1-x0+6*dp,3*dp,z1-z0))
box(main,(x0-3*dp,y1,z0),(x1-x0+6*dp,3*dp,z1-z0))
box(main,(x0-3*dp,y0,z0),(3*dp,y1-y0,z1-z0))
# Extend the piston through the floor and side-wall thickness. The former
# flush-ended piston exposed a leakage path at its three-dimensional corners.
el(main,'setmkbound',mk=1);box(main,(x1,y0-3*dp,z0-3*dp),(4*dp,y1-y0+6*dp,z1-z0+3*dp))
if a.mixed_walls:el(main,'setshapemode').text='actual | bound'
el(main,'setmkbound',mk=2);draw=el(main,'drawfilestl',file=str(hullpath.resolve()),autofill='true',advanced='true');el(draw,'depth',depthmin=.001*dp);el(el(main,'fillpoint',x=-3,y=0,z=-.08),'modefill').text='void | fluid'
fl=el(el(cd,'floatings'),'floating',mkbound=2);el(fl,'massbody',value=p['mass_kg']);el(fl,'center',**dict(zip('xyz',cg)));el(fl,'inertia',**dict(zip('xyz',p['inertia_kg_m2'])))
el(fl,'translationDOF',x=1,y=int(a.free6),z=1);el(fl,'rotationDOF',x=int(a.free6),y=1,z=int(a.free6))
if a.free6:
    el(fl,'linearvelini',x=a.initial_speed,y=0,z=0);el(fl,'angularvelini',x=0,y=0,z=0)
else:
    lv=el(fl,'linearvel');el(lv,'vel',time=0,x=0,y=0,z='none');el(lv,'vel',time=2,x=.5,y=0,z='none');el(lv,'vel',time=max(3,a.duration+1),x=.5,y=0,z='none')
mo=el(el(cd,'motion'),'objreal',ref=1);el(mo,'begin',mov=1,start=0);mv=el(mo,'mvrectfile',id=1,duration=a.duration+1);el(mv,'file',name='piston_motion.dat',fields=2,fieldtime=0,fieldx=1)
# Fluid-facing mDBC geometry, including the translating piston and the hull.
pts=[];faces=[];mks=[]
def quad(q,normal,mk=10):
    q=np.array(q,float);off=len(pts)
    if np.dot(np.cross(q[1]-q[0],q[2]-q[0]),normal)<0:q=q[::-1]
    pts.extend(q.tolist());faces.extend([[off,off+1,off+2],[off,off+2,off+3]]);mks.extend([mk,mk])
zz=z0+dp/2
quad([(x0-3*dp,y0-3*dp,zz),(x1+3*dp,y0-3*dp,zz),(x1+3*dp,y1+3*dp,zz),(x0-3*dp,y1+3*dp,zz)],[0,0,1])
for yy,n in [(y0+dp/2,[0,1,0]),(y1-dp/2,[0,-1,0])]:quad([(x0-3*dp,yy,z0-3*dp),(x1+3*dp,yy,z0-3*dp),(x1+3*dp,yy,z1),(x0-3*dp,yy,z1)],n)
for xx,n,mk in [(x0+dp/2,[1,0,0],10),(x1-dp/2,[-1,0,0],11)]:quad([(xx,y0-3*dp,z0-3*dp),(xx,y1+3*dp,z0-3*dp),(xx,y1+3*dp,z1),(xx,y0-3*dp,z1)],n,mk)
if a.mixed_walls:pts=[];faces=[];mks=[] # Zero-normal tank particles retain DBC; mDBC1 correction acts on hull normals only.
off=len(pts);pts.extend(hm.vertices.tolist());faces.extend((hm.faces+off).tolist());mks.extend([12]*len(hm.faces))
normalpath=ASSET_DIR/f'{a.name}_normals.vtk' if a.fitted else out/'boundary_normals.vtk'
with normalpath.open('w') as vf:
    vf.write('# vtk DataFile Version 3.0\nFluid-facing surfaces\nASCII\nDATASET POLYDATA\n');vf.write(f'POINTS {len(pts)} float\n');np.savetxt(vf,pts,fmt='%.9g')
    vf.write(f'POLYGONS {len(faces)} {len(faces)*4}\n');np.savetxt(vf,np.c_[np.full(len(faces),3),faces],fmt='%d')
    vf.write(f'CELL_DATA {len(faces)}\nSCALARS Mk unsigned_short 1\nLOOKUP_TABLE default\n');np.savetxt(vf,mks,fmt='%d')
ng=el(el(cd,'normals',active='true'),'norgeometry');el(ng,'geometryfile',file=os.path.relpath(normalpath,out/'generated').replace('\\','/'));el(ng,'distanceh',v=2);el(ng,'maxsizeh',v=2.1)
ex=el(root,'execution');sp=el(ex,'special');gauges=el(sp,'gauges');default=el(gauges,'default');el(default,'computedt',value=.01);el(default,'output',value='true');el(default,'outputdt',value=.01)
for xg in [6.,4.5,3.,0.]:
    gu=el(gauges,'swl',name='swl_'+str(xg).replace('.','p'));el(gu,'pointdp',coefdp=.25);el(gu,'point0',x=xg,y=.9,z=-.25);el(gu,'point2',x=xg,y=.9,z=.25)
damp=el(el(sp,'damping'),'dampingzone');el(damp,'limitmin',x=-5.8,y=0,z=0);el(damp,'limitmax',x=-9,y=0,z=0);el(damp,'overlimit',value=.2);el(damp,'redumax',value=10);el(damp,'factorxyz',x=1,y=1,z=1)
pars=el(ex,'parameters')
params={'Boundary':2,'SlipMode':a.slip_mode,'NoPenetration':0 if a.slip_mode==1 else 1,'StepAlgorithm':2,'Kernel':2,'ViscoTreatment':1,'Visco':a.visco,'DensityDT':2,'DensityDTvalue':.1,'Shifting':0,'RigidAlgorithm':0,'FtPause':0,'TimeMax':a.duration,'TimeOut':min(.15,a.duration/5),'CoefDtMin':.01,'DtIni':0,'DtMin':0,'RhopOutMin':700,'RhopOutMax':1300,'MinFluidStop':.99}
for k,val in params.items():el(pars,'parameter',key=k,value=val)
dom=el(pars,'simulationdomain');el(dom,'posmin',x=x0-5*dp,y=y0-5*dp,z=z0-5*dp);el(dom,'posmax',x=x1+.6,y=y1+5*dp,z=z1)
ET.indent(root);ET.ElementTree(root).write(out/'Case_Def.xml',encoding='utf-8',xml_declaration=True)
t=np.arange(0,a.duration+1,.002);ramp=.5*(1-np.cos(np.pi*np.minimum(t/(3*w['period_s']),1)));kh=w['k_rad_m']
transfer=4*np.sinh(kh)**2/(np.sinh(2*kh)+2*kh);stroke=w['height_m']/transfer
np.savetxt(out/'piston_motion.dat',np.c_[t,-.5*stroke*ramp*np.sin(w['omega_rad_s']*t)],fmt='%.10g')
note=dict(status='prepared_not_validated',frame='laboratory',hull='broad_Wigley_proxy',visual_skin='turtle_ship_local_render_only',dp_m=dp,duration_s=a.duration,domain_m=[18,2.4,1.7],
    mass_kg=p['mass_kg'],center_of_mass_m=cg.tolist(),inertia_kg_m2=p['inertia_kg_m2'],surge='0 to 0.5 m/s linear ramp over 2 s; imposed x, free z/pitch',
    boundary='mDBC, no collisions expected or modeled',particle_mass_policy='GenCase default fluid interaction mass; independent physical rigid-body mass',artificial_viscosity_alpha=a.visco,
    target_H_m=.06,wavelength_m=3,limitations=['Coarse 30 mm spacing; only two spacings per target wave height.','Finite tank width and reflections; not an open-ocean performance prediction.','No paddle/sail/propeller forces or detailed turtle-body hydrodynamics.','Run a short pilot only after buoyancy and wave-maker checks.'])
(out/'setup.json').write_text(json.dumps(note,indent=2));print(out)
note.update(hull='fitted_turtle_envelope' if a.fitted else 'broad_Wigley_proxy',target_H_m=a.height,coefh=a.coefh,slip_mode=a.slip_mode,closed_hull_volume_m3=float(hm.volume),hull_file=str(hullpath.resolve()),limitations=['Initial 3-D commissioning diagnostic; not yet a validated ship response.','Finite tank and coarse spacing; reflections and resolution must be checked.','No oar, sail or rudder forces; assumed mass distribution.'])
(out/'setup.json').write_text(json.dumps(note,indent=2))
if a.history or a.free6:
    note.update(hull='history_informed_turtle_envelope' if a.history else note['hull'],rho_kg_m3=c['environment']['rho_kg_m3'],
        surge='Initial speed only; all six rigid-body DOFs free' if a.free6 else note['surge'],free_dofs=c['motion']['free_dofs'],
        initial_velocity_m_s=[a.initial_speed,0,0] if a.free6 else [0,0,0],initial_angular_velocity_rad_s=[0,0,0],
        prescribed_body_velocity=False if a.free6 else True,
        boundary='DBC tank and piston; mDBC1 hull normals only' if a.mixed_walls else note['boundary'],
        scale_ratio=10 if a.history else None,source_model='results/turtle_history_model.json' if a.history else None,
        limitations=['Coarse demonstration, not a validated seakeeping or resistance prediction.','Impulsive initial body velocity in initially still water creates a startup transient.','Finite tank, artificial viscosity, unresolved near-wall flow and no collisions/flooding model.','Head-on waves and symmetry can keep sway, roll and yaw small even though those DOFs are free.'])
    (out/'setup.json').write_text(json.dumps(note,indent=2))
