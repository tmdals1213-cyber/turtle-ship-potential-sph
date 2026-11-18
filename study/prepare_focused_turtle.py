# SPDX-License-Identifier: GPL-3.0-or-later
from project_paths import ASSET_DIR, SOLVER_BIN, FONT_DIR
"""Focused-wave commissioning and free-six-DOF turtle impact case.

First-order NewWave/JONSWAP piston input; the SPH wave is measured separately.
No analytic surface, prescribed ship motion or artificial water jet is used.
"""
from pathlib import Path
import argparse, json, os, xml.etree.ElementTree as ET
import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parent
GEO = ASSET_DIR

def el(parent, tag, **attrs):
    return ET.SubElement(parent, tag, {k:str(v) for k,v in attrs.items()})

def box(parent, point, size):
    b=el(parent,'drawbox'); el(b,'boxfill').text='solid'
    el(b,'point',**dict(zip('xyz',point))); el(b,'size',**dict(zip('xyz',size)))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--name',required=True)
    ap.add_argument('--empty-2d',action='store_true')
    ap.add_argument('--amplitude',type=float,default=.5)
    ap.add_argument('--dp',type=float,default=.04)
    ap.add_argument('--duration',type=float,default=18.)
    ap.add_argument('--focus-x',type=float,default=1.)
    ap.add_argument('--focus-time',type=float,default=12.)
    ap.add_argument('--peak-period',type=float,default=1.88)
    ap.add_argument('--piston-x',type=float,default=12.)
    a=ap.parse_args()
    folder=ROOT/'sph'/a.name
    if (folder/'output').exists():
        raise FileExistsError('Preserve existing solver output; choose a new case name.')
    folder.mkdir(parents=True,exist_ok=True)
    dp=a.dp; g=9.81; depth=1.2; x1=a.piston_x; x0=x1-24.; z0=-depth; z1=1.6
    y0,y1=(0.,0.) if a.empty_2d else (-2.4,2.4)
    model=json.loads((ROOT/'results/turtle_history_model.json').read_text(encoding='utf-8'))
    config=json.loads((ROOT/'case_turtle_history.json').read_text(encoding='utf-8'))
    f=np.linspace(.30,.95,81); omega=2*np.pi*f; fp=1/a.peak_period
    sigma=np.where(f<=fp,.07,.09)
    spectrum=f**-5*np.exp(-1.25*(fp/f)**4)*3.3**np.exp(-.5*((f-fp)/(sigma*fp))**2)
    amp=a.amplitude*spectrum/spectrum.sum()
    k=np.array([brentq(lambda q:q*np.tanh(q*depth)-w*w/g,1e-5,20) for w in omega])
    transfer=4*np.sinh(k*depth)**2/(np.sinh(2*k*depth)+2*k*depth)
    # s=x1-x is positive in the direction of incident propagation (-world x).
    # eta=sum a*cos(k*(s-sf)-omega*(t-tf)); piston displacement is its sine.
    phase=k*(x1-a.focus_x)-omega*a.focus_time
    time=np.arange(0,a.duration+1+.001,.002)
    ramp=.5*(1-np.cos(np.pi*np.minimum(time/2,1)))
    end=.5*(1-np.cos(np.pi*np.minimum(np.maximum(a.duration+1-time,0)/2,1)))
    displacement=-np.sum(amp/transfer*np.sin(time[:,None]*omega+phase),axis=1)*ramp*end
    np.savetxt(folder/'piston_motion.dat',np.c_[time,displacement],fmt='%.10g')
    np.savez(folder/'wave_components.npz',f=f,omega=omega,k=k,amplitude=amp,
             transfer=transfer,phase=phase,time=time,piston_x=displacement)
    if np.max(np.abs(displacement))>.85:
        raise ValueError('Piston motion exceeds the extended-floor allowance.')
    root=ET.Element('case'); cd=el(root,'casedef'); con=el(cd,'constantsdef')
    el(con,'gravity',x=0,y=0,z=-g)
    for tag,v,attrs in [('rhop0',1025,{}),('rhopgradient',2,{}),('hswl',depth,{'auto':'false'}),
        ('gamma',7,{}),('speedsystem',3,{'auto':'false'}),('coefsound',20,{}),
        ('speedsound',60,{'auto':'false'}),('coefh',1.2,{}),('cflnumber',.2,{})]:
        el(con,tag,value=v,**attrs)
    el(cd,'mkconfig',boundcount=240,fluidcount=9)
    geom=el(cd,'geometry'); definition=el(geom,'definition',dp=dp)
    el(definition,'pointmin',x=x0-5*dp,y=0 if a.empty_2d else y0-5*dp,z=z0-5*dp)
    el(definition,'pointmax',x=x1+1.,y=0 if a.empty_2d else y1+5*dp,z=z1)
    commands=el(el(geom,'commands'),'mainlist')
    el(commands,'setshapemode').text='dp | bound'; el(commands,'setdrawmode',mode='full')
    el(commands,'setmkfluid',mk=0); box(commands,(x0,y0,z0),(x1-x0,y1-y0,depth))
    # The stationary floor and side walls extend behind the complete piston stroke.
    yy0,yy1=(0.,0.) if a.empty_2d else (y0-3*dp,y1+3*dp)
    el(commands,'setmkbound',mk=1)
    box(commands,(x1,yy0,z0-3*dp),(4*dp,yy1-yy0,z1-z0+3*dp))
    # Static walls overwrite shared points. Otherwise the piston steals floor
    # particles and opens a hole at its original footprint when it retreats.
    el(commands,'setmkbound',mk=0)
    box(commands,(x0-3*dp,yy0,z0-3*dp),(x1+1-x0+3*dp,yy1-yy0,3*dp))
    box(commands,(x0-3*dp,y0,z0),(3*dp,y1-y0,z1-z0))
    if not a.empty_2d:
        for yy in [y0-3*dp,y1]:
            box(commands,(x0-3*dp,yy,z0),(x1+1-x0+3*dp,3*dp,z1-z0))
    if not a.empty_2d:
        import trimesh
        hm=trimesh.load(GEO/'turtle_history_closed.stl',process=True)
        assert hm.is_watertight
        el(commands,'setshapemode').text='actual | bound'; el(commands,'setmkbound',mk=2)
        draw=el(commands,'drawfilestl',file=str((GEO/'turtle_history_closed.stl').resolve()),autofill='true',advanced='true')
        el(draw,'depth',depthmin=.001*dp)
        el(el(commands,'fillpoint',x=0,y=0,z=-.08),'modefill').text='void | fluid'
        fl=el(el(cd,'floatings'),'floating',mkbound=2)
        el(fl,'massbody',value=model['mass_kg'])
        el(fl,'center',**dict(zip('xyz',model['center_of_mass_m'])))
        el(fl,'inertia',**dict(zip('xyz',model['inertia_kg_m2'])))
        el(fl,'translationDOF',x=1,y=1,z=1); el(fl,'rotationDOF',x=1,y=1,z=1)
        normals=GEO/f'{a.name}_normals.vtk'
        with normals.open('w') as stream:
            stream.write('# vtk DataFile Version 3.0\nHull fluid-facing normals\nASCII\nDATASET POLYDATA\n')
            stream.write(f'POINTS {len(hm.vertices)} float\n'); np.savetxt(stream,hm.vertices,fmt='%.9g')
            stream.write(f'POLYGONS {len(hm.faces)} {4*len(hm.faces)}\n')
            np.savetxt(stream,np.c_[np.full(len(hm.faces),3),hm.faces],fmt='%d')
            stream.write(f'CELL_DATA {len(hm.faces)}\nSCALARS Mk unsigned_short 1\nLOOKUP_TABLE default\n')
            np.savetxt(stream,np.full(len(hm.faces),12),fmt='%d')
        ng=el(el(cd,'normals',active='true'),'norgeometry')
        el(ng,'geometryfile',file=os.path.relpath(normals,folder/'generated').replace('\\','/'))
        el(ng,'distanceh',v=2); el(ng,'maxsizeh',v=2.1)
    mo=el(el(cd,'motion'),'objreal',ref=1); el(mo,'begin',mov=1,start=0)
    el(el(mo,'mvrectfile',id=1,duration=a.duration+1),'file',name='piston_motion.dat',fields=2,fieldtime=0,fieldx=1)
    ex=el(root,'execution'); special=el(ex,'special'); gauges=el(special,'gauges')
    default=el(gauges,'default'); el(default,'computedt',value=.01)
    el(default,'output',value='true'); el(default,'outputdt',value=.01)
    for xx in [8.,6.,4.,3.,2.,1.,0.,-1.,-2.]:
        if xx>=x1-.8:continue
        gauge=el(gauges,'swl',name='swl_'+str(xx).replace('-','m').replace('.','p'))
        el(gauge,'pointdp',coefdp=.25)
        el(gauge,'point0',x=xx,y=0 if a.empty_2d else 1.,z=-.7)
        el(gauge,'point2',x=xx,y=0 if a.empty_2d else 1.,z=1.3)
    damping=el(el(special,'damping'),'dampingzone')
    el(damping,'limitmin',x=x0+5,y=0,z=0);el(damping,'limitmax',x=x0,y=0,z=0)
    el(damping,'overlimit',value=.2);el(damping,'redumax',value=10);el(damping,'factorxyz',x=1,y=1,z=1)
    params={'Boundary':1 if a.empty_2d else 2,'SlipMode':1,'NoPenetration':0,
        'StepAlgorithm':2,'Kernel':2,'ViscoTreatment':1,'Visco':.01,'DensityDT':2,
        'DensityDTvalue':.1,'Shifting':0,'RigidAlgorithm':0,'FtPause':0,
        'TimeMax':a.duration,'TimeOut':.08,'CoefDtMin':.01,'DtIni':0,'DtMin':0,
        'RhopOutMin':700,'RhopOutMax':1300,'MinFluidStop':.99}
    pars=el(ex,'parameters')
    for key,value in params.items():el(pars,'parameter',key=key,value=value)
    dom=el(pars,'simulationdomain')
    el(dom,'posmin',x=x0-5*dp,y=0 if a.empty_2d else y0-5*dp,z=z0-5*dp)
    el(dom,'posmax',x=x1+1.,y=0 if a.empty_2d else y1+5*dp,z=z1)
    ET.indent(root);ET.ElementTree(root).write(folder/'Case_Def.xml',encoding='utf-8',xml_declaration=True)
    note=dict(status='prepared_not_wave_calibrated',purpose='Empty-tank focused-wave commissioning' if a.empty_2d else 'Turtle ship extreme-wave intro',
        dimension=2 if a.empty_2d else 3,dp_m=dp,duration_s=a.duration,domain_m=[24,y1-y0,z1-z0],water_depth_m=depth,
        wave_type='First-order crest-focused NewWave using JONSWAP gamma=3.3; actual nonlinear wave measured from SPH',
        target_linear_crest_amplitude_m=a.amplitude,peak_period_s=a.peak_period,frequency_band_hz=[float(f[0]),float(f[-1])],
        focus_x_m=a.focus_x,focus_time_s=a.focus_time,piston_x_m=x1,wave_direction='negative world x',
        piston_displacement_range_m=[float(displacement.min()),float(displacement.max())],
        piston_max_speed_m_s=float(abs(np.gradient(displacement,time)).max()),
        boundary='DBC tank and piston; mDBC1 hull normals only' if not a.empty_2d else 'DBC',
        scale_ratio=10,rho_kg_m3=1025,initial_velocity_m_s=[0,0,0],prescribed_body_velocity=False,
        free_dofs=[] if a.empty_2d else ['Surge','Sway','Heave','Roll','Pitch','Yaw'],
        sources=['https://github.com/DualSPHysics/DualSPHysics/wiki/3.-SPH-formulation',
          'https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/green-water-loads-on-prismatic-obstacles/DF5D3EB421E38AEF11D70DFD967BEC12'],
        limitations=['Target linear crest amplitude is not measured wave height or significant wave height.',
          'First-order piston forcing; nonlinear shift of focus/breaking must be measured in the empty tank.',
          'Finite-width tank with downstream damping; reflected waves are not separated, especially in later events.',
          'Coarse demonstration: no validated peak impact pressures, air entrapment, structural damage or internal flooding.',
          'Upper closed hull interacts with water; mast, oars and decorative appendages remain display-only.',
          'Existing 40 mm hull hydrostatic bias is not removed by adjusting the physical mass.'])
    if not a.empty_2d:
        note.update(mass_kg=model['mass_kg'],center_of_mass_m=model['center_of_mass_m'],inertia_kg_m2=model['inertia_kg_m2'],
                    closed_hull_volume_m3=float(hm.volume),hull_file=str((GEO/'turtle_history_closed.stl').resolve()))
    (folder/'setup.json').write_text(json.dumps(note,indent=2),encoding='utf-8')
    print(json.dumps(note,indent=2))

if __name__=='__main__':main()
