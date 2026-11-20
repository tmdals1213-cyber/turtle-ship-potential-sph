"""Matched Gunnerus calm/head-sea pilots and an independent 2-D wavemaker check.

Reuses the commissioned closed hull, six free DOFs and assumed CG force.
No actuator disk or propeller wake is implemented by this input generator.
"""
from pathlib import Path
import argparse, json, os, xml.etree.ElementTree as ET
import numpy as np
from scipy.optimize import brentq
from prepare_gunnerus import ROOT, SCALE, RHO, el, box


def prepare(kind, *, name=None, dp=.05, depth=1.25, duration=10., domain=None, physical_model=None):
    empty = kind == 'empty'
    name = name or f'gunnerus_{kind}_050_10s'
    folder = ROOT/'sph'/name
    if folder.exists():
        raise FileExistsError(f'Preserving existing case: {folder}')
    folder.mkdir()
    p = physical_model or json.loads((ROOT/'geometry/gunnerus/physical_model.json').read_text())
    domain = domain or [-12.,12.,-3.,3.,-depth,1.6]
    x0, x1, z1 = domain[0],domain[1],domain[5]
    y0, y1 = (0., 0.) if empty else (domain[2],domain[3])
    height, period = .4, 2.2
    omega = 2*np.pi/period
    k = brentq(lambda q: 9.81*q*np.tanh(q*depth)-omega**2, .001, 20.)
    transfer = 4*np.sinh(k*depth)**2/(np.sinh(2*k*depth)+2*k*depth)
    t = np.arange(0, duration+1.001, .002)
    ramp_time = 2*period
    ramp = .5*(1-np.cos(np.pi*np.minimum(t/ramp_time, 1)))
    displacement = -height/(2*transfer)*ramp*np.sin(omega*t)
    if kind == 'calm':
        displacement[:] = 0
    np.savetxt(folder/'piston_motion.dat', np.c_[t,displacement], fmt='%.10g')
    root = ET.Element('case'); cd = el(root,'casedef'); co = el(cd,'constantsdef')
    el(co,'gravity',x=0,y=0,z=-9.81)
    for tag,value,attrs in [('rhop0',RHO,{}),('rhopgradient',2,{}),('hswl',depth,{'auto':'false'}),
        ('gamma',7,{}),('speedsystem',2,{'auto':'false'}),('coefsound',20,{}),
        ('speedsound',60,{'auto':'false'}),('coefh',1.2,{}),('cflnumber',.2,{})]:
        el(co,tag,value=value,**attrs)
    el(cd,'mkconfig',boundcount=240,fluidcount=9)
    geom = el(cd,'geometry'); de = el(geom,'definition',dp=dp)
    bounds_min = dict(x=x0-5*dp,y=0 if empty else y0-5*dp,z=-depth-5*dp-dp/2)
    bounds_max = dict(x=x1+1,y=0 if empty else y1+5*dp,z=z1+dp)
    el(de,'pointmin',**bounds_min); el(de,'pointmax',**bounds_max)
    cmd = el(el(geom,'commands'),'mainlist')
    el(cmd,'setshapemode').text = 'dp | bound'; el(cmd,'setdrawmode',mode='full')
    el(cmd,'setmkfluid',mk=0); box(cmd,(x0,y0,-depth),(x1-x0,y1-y0,depth-dp/2))
    yy0, yy1 = (0.,0.) if empty else (y0-3*dp,y1+3*dp)
    el(cmd,'setmkbound',mk=1)
    box(cmd,(x1,yy0,-depth-3*dp),(4*dp,yy1-yy0,z1+depth+3*dp))
    # Fixed walls overwrite intersections, retaining a continuous floor when
    # the piston moves back. The floor extends behind the complete stroke.
    el(cmd,'setmkbound',mk=0)
    box(cmd,(x0-3*dp,yy0,-depth-3*dp),(x1+1-x0+3*dp,yy1-yy0,3*dp))
    box(cmd,(x0-3*dp,y0,-depth),(3*dp,y1-y0,z1+depth))
    if not empty:
        for yy in [y0-3*dp,y1]:
            box(cmd,(x0-3*dp,yy,-depth),(x1+1-x0+3*dp,3*dp,z1+depth))
        el(cmd,'setshapemode').text = 'actual | bound'; el(cmd,'setmkbound',mk=2)
        stl = el(cmd,'drawfilestl',file=p['hull_file'],autofill='true',advanced='true')
        el(stl,'depth',depthmin=.001*dp)
        el(el(cmd,'fillpoint',x=0,y=0,z=-.1),'modefill').text = 'void | fluid'
        body = el(el(cd,'floatings'),'floating',mkbound=2)
        el(body,'massbody',value=p['mass_kg'])
        el(body,'center',**dict(zip('xyz',p['center_of_mass_m'])))
        el(body,'inertia',**dict(zip('xyz',p['inertia_kg_m2'])))
        el(body,'translationDOF',x=1,y=1,z=1); el(body,'rotationDOF',x=1,y=1,z=1)
        el(body,'linearvelini',x=8*.514444/np.sqrt(SCALE),y=0,z=0)
        force = el(body,'linearforce')
        for tt,ff in [(0,0),(.3,0),(1.3,30),(duration+1,30)]:
            el(force,'force',time=tt,x=ff,y=0,z=0)
        ng = el(el(cd,'normals',active='true'),'norgeometry')
        el(ng,'geometryfile',file=os.path.relpath(p['normals_file'],folder/'generated').replace('\\','/'))
        el(ng,'distanceh',v=2); el(ng,'maxsizeh',v=2.1)
    mo = el(el(cd,'motion'),'objreal',ref=1); el(mo,'begin',mov=1,start=0)
    el(el(mo,'mvrectfile',id=1,duration=duration+1),'file',name='piston_motion.dat',fields=2,fieldtime=0,fieldx=1)
    ex = el(root,'execution'); special = el(ex,'special'); gauges = el(special,'gauges')
    default = el(gauges,'default'); el(default,'computedt',value=.01)
    el(default,'output',value='true'); el(default,'outputdt',value=.01)
    for xx in [10.,8.,6.,4.,2.,0.]:
        gu = el(gauges,'swl',name=f'x{int(xx)}')
        el(gu,'pointdp',coefdp=.25)
        el(gu,'point0',x=xx,y=0 if empty else 2.,z=-.6)
        el(gu,'point2',x=xx,y=0 if empty else 2.,z=1.1)
    damp = el(el(special,'damping'),'dampingzone')
    el(damp,'limitmin',x=-7,y=0,z=0); el(damp,'limitmax',x=x0,y=0,z=0)
    el(damp,'overlimit',value=.2); el(damp,'redumax',value=10)
    el(damp,'factorxyz',x=1,y=1,z=1)
    pars = el(ex,'parameters')
    values = dict(Boundary=1 if empty else 2,SlipMode=1,NoPenetration=0,
        StepAlgorithm=2,Kernel=2,ViscoTreatment=1,Visco=.01,DensityDT=2,DensityDTvalue=.1,
        Shifting=0,RigidAlgorithm=0,FtPause=0,FtSaveAce=1,TimeMax=duration,TimeOut=.1,
        TimeOutExtra=.01,CoefDtMin=.01,DtIni=0,DtMin=0,RhopOutMin=700,RhopOutMax=1300,MinFluidStop=.99)
    for key,value in values.items(): el(pars,'parameter',key=key,value=value)
    dom = el(pars,'simulationdomain'); el(dom,'posmin',**bounds_min); el(dom,'posmax',**bounds_max)
    ET.indent(root); ET.ElementTree(root).write(folder/'Case_Def.xml',encoding='utf-8',xml_declaration=True)
    setup = dict(status='prepared_not_run',case=name,kind=kind,dp_m=dp,duration_s=duration,
        scale_ratio=SCALE,dimension=2 if empty else 3,physical_model=None if empty else p,
        domain_bounds_m=[x0,x1,y0,y1,-depth,z1],water_depth_m=depth,
        target_H_m=0. if kind=='calm' else height,intrinsic_period_s=period,
        wavelength_m=float(2*np.pi/k),wave_direction='negative world x',
        wave_input='First-order Biesel piston, two-period cosine ramp; actual wave measured independently.',
        wave_ramp_s=ramp_time,piston_displacement_range_m=[float(displacement.min()),float(displacement.max())],
        initial_velocity_m_s=[0 if empty else 8*.514444/np.sqrt(SCALE),0,0],all_six_dofs_free=not empty,
        propulsion=dict(mode='none' if empty else 'world-fixed external force at CG',plateau_force_N=[0 if empty else 30,0,0],
            acts_on_fluid=False,body_following=False,resolved_propeller=False,basis='Assumed diagnostic force, not a calibrated propulsion operating point'),
        output=dict(particles_interval_s=.1,native_force_motion_extra_interval_s=.01),gpu_wall_limit_s=None,
        limitations=([] if empty else p['limitations'])+[
            'Coarse pilot; not validated resistance, propulsion efficiency or impact pressure.',
            'Finite tank with downstream damping; later reflections and finite depth remain part of this setup.',
            'Matched force and initial speed do not maintain equal speeds, so differences are not equal-speed added resistance.',
            'Assumed thrust does not add propeller momentum to the water; no actuator disk is present.',
            '10-second transient includes wavemaker ramp and ship startup.'])
    (folder/'setup.json').write_text(json.dumps(setup,indent=2),encoding='utf-8')
    return name


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--kind',choices=['empty','calm','waves','all'],default='all')
    args = ap.parse_args()
    print(json.dumps([prepare(k) for k in (['empty','calm','waves'] if args.kind=='all' else [args.kind])]))
