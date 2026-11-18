# SPDX-License-Identifier: GPL-3.0-or-later
"""Generate a reproducible, deliberately coarse GPU startup case and wave-test inputs.

The inlet is a Z/time grid (upstream GPU supports one X column). This is a
preparation/pilot configuration, not yet an accepted seakeeping calculation.
"""
from pathlib import Path
import argparse,json,sys,shutil,xml.etree.ElementTree as ET
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from common import ROOT,config,physics,wave,airy

def el(p,tag,**attrs):
    return ET.SubElement(p,tag,{k:str(v) for k,v in attrs.items()})
def box(p,point,size,fill='solid'):
    b=el(p,'drawbox');el(b,'boxfill').text=fill
    el(b,'point',**dict(zip('xyz',point)));el(b,'size',**dict(zip('xyz',size)))

def make(name,dp,duration,hull=True,rest=False,mdbc=False):
    c=config()
    if rest:
        c['motion']['forward_speed_m_s']=0;c['waves']['baseline_height_m']=0
        c['sph']['domain_m']=dict(xmin=-2.1,xmax=2.1,ymin=-.6,ymax=.6,zmin=-1,zmax=.35)
    p=physics(c);w=wave(c=c);e=c['environment'];s=c['sph'];d=s['domain_m']
    out=ROOT/'sph'/name;out.mkdir(parents=True,exist_ok=True)
    xmin,xmax,ymin,ymax,zmin,zmax=[d[k] for k in ['xmin','xmax','ymin','ymax','zmin','zmax']]
    U=c['motion']['forward_speed_m_s'];buffer=6*dp;maxwave=.15
    root=ET.Element('case');cd=el(root,'casedef');con=el(cd,'constantsdef')
    el(con,'gravity',x=0,y=0,z=-e['g_m_s2'])
    for tag,value,extra in [('rhop0',e['rho_kg_m3'],{}),('rhopgradient',2,{}),('hswl',-zmin,{'auto':'false'}),('gamma',7,{}),('speedsystem',U,{'auto':'false'}),('coefsound',20,{}),('speedsound',s['sound_speed_m_s'],{'auto':'false'}),('coefh',1,{}),('cflnumber',s['cfl'],{})]:el(con,tag,value=value,**extra)
    el(cd,'mkconfig',boundcount=240,fluidcount=9)
    geom=el(cd,'geometry');de=el(geom,'definition',dp=dp)
    el(de,'pointmin',x=xmin-buffer-4*dp,y=ymin-4*dp,z=zmin-4*dp)
    el(de,'pointmax',x=xmax+buffer+4*dp,y=ymax+4*dp,z=zmax)
    main=el(el(geom,'commands'),'mainlist');el(main,'setshapemode').text='actual | bound' if mdbc else 'dp | bound'
    el(main,'setdrawmode',mode='full');el(main,'setmkfluid',mk=0)
    box(main,(xmin,ymin,zmin),(xmax-xmin,ymax-ymin,-zmin))
    el(main,'setmkbound',mk=0)
    length=xmax-xmin+2*buffer
    box(main,(xmin-buffer,ymin-3*dp,zmin-3*dp),(length,ymax-ymin+6*dp,3*dp))
    box(main,(xmin-buffer,ymin-3*dp,zmin),(length,3*dp,zmax-zmin))
    box(main,(xmin-buffer,ymax,zmin),(length,3*dp,zmax-zmin))
    if rest:
        box(main,(xmin-3*dp,ymin,zmin),(3*dp,ymax-ymin,zmax-zmin))
        box(main,(xmax,ymin,zmin),(3*dp,ymax-ymin,zmax-zmin))
    if hull:
        shutil.copyfile(ROOT/'geometry/wigley_closed.stl',out/'hull.stl')
        el(main,'setmkbound',mk=1)
        draw=el(main,'drawfilestl',file='hull.stl',autofill='true',**({'advanced':'true'} if mdbc else {}))
        if mdbc:el(draw,'depth',depthmin=.001*dp)
        # STL autofill closes rasterized boundary gaps; it does not replace
        # pre-existing fluid enclosed by the shell. Fill the closed interior
        # with the same floating marker before assigning the physical mass.
        inside=el(main,'fillpoint',x=0,y=0,z=-.05)
        el(inside,'modefill').text='void | fluid'
        f=el(el(cd,'floatings'),'floating',mkbound=1)
        # Physical rigid-body mass is independent of SPH quadrature mass.
        # Leave masspart at GenCase's fluid-particle-mass default, per upstream.
        el(f,'massbody',value=p['mass_kg'])
        el(f,'center',**dict(zip('xyz',c['hull']['center_of_mass_m'])))
        el(f,'inertia',**dict(zip('xyz',p['inertia_kg_m2'])))
        el(f,'translationDOF',x=0,y=0,z=1);el(f,'rotationDOF',x=0,y=1,z=0)
    if mdbc:
        # Supply only the fluid-facing wall surfaces. Including the outer sides
        # of thick solid wall boxes gives wrong ghost-node directions.
        import trimesh
        points=[];faces=[]
        def quad(q,normal):
            off=len(points);q=np.asarray(q,float)
            if np.dot(np.cross(q[1]-q[0],q[2]-q[0]),normal)<0:q=q[::-1]
            points.extend(q.tolist());faces.extend([[off,off+1,off+2],[off,off+2,off+3]])
        lo,hi=xmin+dp/2,xmax-dp/2;near,far=ymin+dp/2,ymax-dp/2;floor=zmin+dp/2
        quad([(xmin-buffer,ymin-3*dp,floor),(xmax+buffer,ymin-3*dp,floor),(xmax+buffer,ymax+3*dp,floor),(xmin-buffer,ymax+3*dp,floor)],[0,0,1])
        for yy,n in [(near,[0,1,0]),(far,[0,-1,0])]:quad([(xmin-buffer,yy,zmin-3*dp),(xmax+buffer,yy,zmin-3*dp),(xmax+buffer,yy,zmax),(xmin-buffer,yy,zmax)],n)
        if rest:
            for xx,n in [(lo,[1,0,0]),(hi,[-1,0,0])]:quad([(xx,ymin-3*dp,zmin-3*dp),(xx,ymax+3*dp,zmin-3*dp),(xx,ymax+3*dp,zmax),(xx,ymin-3*dp,zmax)],n)
        mk=[10]*len(faces)
        if hull:
            hm=trimesh.load(out/'hull.stl');off=len(points);points.extend(hm.vertices.tolist());faces.extend((hm.faces+off).tolist());mk.extend([11]*len(hm.faces))
        with (out/'boundary_normals.vtk').open('w') as vf:
            vf.write('# vtk DataFile Version 3.0\nFluid-facing surfaces\nASCII\nDATASET POLYDATA\n')
            vf.write(f'POINTS {len(points)} float\n');np.savetxt(vf,points,fmt='%.9g')
            vf.write(f'POLYGONS {len(faces)} {4*len(faces)}\n');np.savetxt(vf,np.c_[np.full(len(faces),3),faces],fmt='%d')
            vf.write(f'CELL_DATA {len(faces)}\nSCALARS Mk unsigned_short 1\nLOOKUP_TABLE default\n');np.savetxt(vf,mk,fmt='%d')
        ng=el(el(cd,'normals',active='true'),'norgeometry')
        el(ng,'geometryfile',file='../boundary_normals.vtk');el(ng,'distanceh',v=2);el(ng,'maxsizeh',v=2.1)
    el(el(cd,'initials'),'velocity',mkfluid=0,x=-U,y=0,z=0)
    ex=el(root,'execution');sp=el(ex,'special')
    gauges=el(sp,'gauges');default=el(gauges,'default')
    el(default,'computedt',value=.01);el(default,'output',value='true');el(default,'outputdt',value=.01)
    for xg in [4.5,3.75,3.0,0.0,-3.0]:
        g=el(gauges,'swl',name='swl_'+str(xg).replace('-','m').replace('.','p'))
        el(g,'pointdp',coefdp=.25);el(g,'point0',x=xg,y=.45 if hull else 0,z=-.3);el(g,'point2',x=xg,y=.45 if hull else 0,z=.3)
    io=el(sp,'inout')
    el(io,'memoryresize',size0=2,size=4);el(io,'extrapolatemode',value=1)
    for inlet,x,direction in [(True,xmax,-1),(False,xmin,1)]:
        zone=el(io,'inoutzone');el(zone,'refilling',value=2);el(zone,'inputtreatment',value=2 if inlet else 1);el(zone,'layers',value=6)
        b=el(el(zone,'zone3d'),'box')
        el(b,'point',x=x,y=ymin+dp,z=zmin+dp);el(b,'size',x=0,y=ymax-ymin-2*dp,z=maxwave-zmin-dp)
        el(b,'direction',x=direction,y=0,z=0)
        vel=el(zone,'imposevelocity',mode=3 if inlet else 2)
        if inlet:
            el(vel,'gridveldata',file='inlet_velocity.csv');el(vel,'gridposzero',x=xmax,z=zmin)
        el(zone,'imposerhop',mode=1 if inlet else 2)
        zs=el(zone,'imposezsurf',mode=1 if inlet else 2);el(zs,'savevtk',value='false');el(zs,'remove',value='false')
        if inlet:el(zs,'zsurffile',file='inlet_surface.csv')
        else:el(zs,'zsurf0',value=0);el(zs,'zsurfmin',value=-.2);el(zs,'zsurffit',valuedp=.5)
    params={'SavePosDouble':0,'StepAlgorithm':2,'Kernel':2,'ViscoTreatment':1,'Visco':.01,'ViscoBoundFactor':1,'DensityDT':2,'DensityDTvalue':.1,'Shifting':0,'RigidAlgorithm':1,'FtPause':0,'CoefDtMin':.05,'DtIni':0,'DtMin':0,'DtFixed':0,'DtAllParticles':0,'TimeMax':duration,'TimeOut':min(.05,duration/5),'PartsOutMax':1,'RhopOutMin':700,'RhopOutMax':1300}
    if rest:sp.remove(io)
    if mdbc:params.update(Boundary=2,SlipMode=2,NoPenetration=1)
    pars=el(ex,'parameters')
    for k,v in params.items():el(pars,'parameter',key=k,value=v)
    dom=el(pars,'simulationdomain');el(dom,'posmin',x=xmin-buffer-4*dp,y=ymin-4*dp,z=zmin-4*dp);el(dom,'posmax',x=xmax+buffer+4*dp,y=ymax+4*dp,z=zmax)
    ET.indent(root);ET.ElementTree(root).write(out/'Case_Def.xml',encoding='utf-8',xml_declaration=True)
    # Sample beyond the requested run so callers can extend a pilot to a wave test.
    times=np.arange(0,max(duration+1,20),.01);z=np.linspace(zmin,maxwave,101)
    ramp=.5*(1-np.cos(np.pi*np.minimum(times/(3*w['encounter_period_s']),1)))
    eta,u,vz=airy(xmax,z[None,:],times[:,None],w,c)
    u=-U+(u+U)*ramp[:,None];vz*=ramp[:,None]
    with (out/'inlet_velocity.csv').open('w') as f:
        f.write('fmtversion;grid_dpx;grid_dpz;grid_nx;grid_nz;vars\n')
        f.write(f'1;{dp};{z[1]-z[0]};1;{len(z)};velx velz\n\n')
        f.write(';'.join(['time']+[f'vx_z{i}' for i in range(len(z))]+[f'vz_z{i}' for i in range(len(z))])+'\n')
        np.savetxt(f,np.c_[times,u,vz],delimiter=';',fmt='%.9g')
    np.savetxt(out/'inlet_surface.csv',np.c_[times,eta[:,0]*ramp],delimiter=';',header='time;zsurf',fmt='%.9g')
    note=dict(status='generated_not_wave_validated',dp_m=dp,duration_s=duration,hull=hull,rest=rest,boundary='mDBC' if mdbc else 'DBC',solid_interior_fill=hull,particle_mass_policy='default fluid particle mass for interaction; specified physical rigid-body mass',
        free_dofs=['Heave','Pitch'] if hull else [],
        checks_needed=['Empty-flume wave height/phase and reflection before RAO comparison.','Hydrostatic float equilibrium and DBC effective draft.','Lateral extent, bottom treatment and resolution sensitivity.'],
        limitations=['Artificial viscosity alpha=0.01; not a calibrated viscous-resistance model.','One inlet X column: orbital profile is uniform across the 6dp inlet buffer.','Stationary DBC bottom/side walls: local current flume, not exact unbounded translating sea.','Extrapolated outlet without an absorption zone; reflected waves must be assessed.','Startup timing does not establish periodic motion or converged loads.'])
    (out/'setup.json').write_text(json.dumps(note,indent=2));print(out)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--name',default='smoke_hull');a.add_argument('--dp',type=float,default=.06);a.add_argument('--duration',type=float,default=.05);a.add_argument('--empty',action='store_true')
    a.add_argument('--rest',action='store_true');a.add_argument('--mdbc',action='store_true')
    q=a.parse_args();make(q.name,q.dp,q.duration,not q.empty,q.rest,q.mdbc)
