"""Prepare matched-mass cylindrical-tanker GPU SPH comparisons.

Analytic Cartesian particle masks provide auditable geometry. No solver edits.
"""
from pathlib import Path
import json
import os
import subprocess
import xml.etree.ElementTree as ET
import numpy as np
import pyvista as pv

ROOT=Path(__file__).resolve().parent
from tanker_paths import BIN
DP=.04
CASES=['A_unpartitioned','C_separate_compartments','B_perforated_baffles','D_cross_baffles']


def el(p,tag,**attrs):
    return ET.SubElement(p,tag,{k:str(v) for k,v in attrs.items()})


def profile():
    a=.4*9.81;u0=50/3.6;ramp=.4;start=1.5;duration=u0/a+ramp;stop=start+duration
    times=np.unique(np.r_[np.arange(0,12.00001,.01),start,start+ramp,stop-ramp,stop,12.])
    q=times-start
    accel=np.where((q>=0)&(q<ramp),a*.5*(1-np.cos(np.pi*q/ramp)),0.)
    accel=np.where((q>=ramp)&(q<=duration-ramp),a,accel)
    accel=np.where((q>duration-ramp)&(q<duration),a*.5*(1+np.cos(np.pi*(q-duration+ramp)/ramp)),accel)
    return times,accel,{'initial_speed_m_s':u0,'peak_deceleration_m_s2':a,'brake_start_s':start,'brake_stop_s':stop,'ramp_s':ramp}


def lattice(plate_halfwidth=.06, fluid_count=None):
    xx=(np.arange(-78,78)+.5)*DP
    yy=(np.arange(-28,28)+.5)*DP
    zz=yy.copy()
    x,y,z=np.meshgrid(xx,yy,zz,indexing='ij')
    p=np.column_stack((x.ravel(),y.ravel(),z.ravel()))
    x,y,z=p.T;r2=y*y+z*z
    outer=(abs(x)<3.12+1e-7)&(r2<(1.12+1e-7)**2)
    inner=(abs(x)<3)&(r2<1.)
    walls=outer&~inner
    out={}
    for case in CASES:
        plates=np.zeros(len(p),bool)
        if case!='A_unpartitioned':
            trans=(np.minimum(abs(x-1),abs(x+1))<plate_halfwidth-1e-7)&inner
            if case!='C_separate_compartments':
                holes=np.zeros(len(p),bool)
                for cy in (-.42,.42):
                    for cz in (-.4,.4):holes|=(y-cy)**2+(z-cz)**2<.24**2
                trans&=~holes
            plates|=trans
        if case=='D_cross_baffles':plates|=(abs(y)<plate_halfwidth-1e-7)&inner
        out[case]={'walls':p[walls],'plates':p[plates],'available':p[inner&~plates]}
    # Keep exact equal particle masses. Four-fold symmetry balances x/y.
    n=fluid_count if fluid_count is not None else min(np.count_nonzero(v['available'][:,2]<0) for v in out.values())//4*4
    assert n%4==0
    for data in out.values():
        q=data.pop('available')
        positive=q[(q[:,0]>0)&(q[:,1]>0)]
        ix=np.lexsort((positive[:,0],positive[:,1],positive[:,0]**2+positive[:,1]**2,positive[:,2]))
        chosen=positive[ix[:n//4]]
        data['fluid']=np.concatenate([chosen*np.array([sx,sy,1]) for sx,sy in ((1,1),(-1,1),(1,-1),(-1,-1))])
    return out,n


def prepare(case,data,n,suffix='brake_040_v1',plate_thickness=.08):
    f=ROOT/'sph'/f'{case}_{suffix}'
    if (f/'initial_audit.json').exists():
        audit=json.loads((f/'initial_audit.json').read_text())
        if audit['passed']:return json.loads((f/'setup.json').read_text()),audit
    if (f/'output').exists():raise FileExistsError(f'Preserve simulation {f}')
    (f/'generated').mkdir(parents=True,exist_ok=True)
    for label,points in data.items():
        if len(points):
            with (f/f'{label}.vtk').open('w') as stream:
                stream.write('# vtk DataFile Version 3.0\nAnalytic tanker particles\nASCII\nDATASET POLYDATA\n')
                stream.write(f'POINTS {len(points)} float\n')
                np.savetxt(stream,points,fmt='%.8f')
                stream.write(f'VERTICES {len(points)} {2*len(points)}\n')
                np.savetxt(stream,np.column_stack((np.ones(len(points),int),np.arange(len(points)))),fmt='%d')
    root=ET.Element('case');cd=el(root,'casedef');co=el(cd,'constantsdef')
    el(co,'gravity',x=0,y=0,z=-9.81)
    for tag,v,extra in [('rhop0',1000,{}),('rhopgradient',2,{}),('hswl',1,{'auto':'false'}),('gamma',7,{}),
                        ('speedsystem',5,{'auto':'false'}),('coefsound',14,{}),('speedsound',70,{'auto':'false'}),
                        ('coefh',1.2,{}),('cflnumber',.25,{})]:el(co,tag,value=v,**extra)
    el(cd,'mkconfig',boundcount=240,fluidcount=9)
    ge=el(cd,'geometry');de=el(ge,'definition',dp=DP)
    el(de,'pointmin',x=-3.26,y=-1.26,z=-1.26);el(de,'pointmax',x=3.26,y=1.26,z=1.26)
    cmd=el(el(ge,'commands'),'mainlist');el(cmd,'setshapemode').text='dp | bound';el(cmd,'setdrawmode',mode='full')
    for label,tag,mk in [('fluid','setmkfluid',0),('walls','setmkbound',0),('plates','setmkbound',1)]:
        if not len(data[label]):continue
        el(cmd,tag,mk=mk);el(el(cmd,'drawfilevtk',file=str(f/f'{label}.vtk')),'polyselec').text='points'
    ex=el(root,'execution');sp=el(ex,'special');acc=el(el(sp,'accinputs'),'accinput',mkfluid=0)
    el(acc,'acccentre',x=0,y=0,z=0);el(acc,'globalgravity',value=1)
    tt,aa,brake=profile();series=el(acc,'acctimes')
    for t,a in zip(tt,aa):el(series,'timevalue',time=f'{t:.9f}',linx=f'{a:.9f}')
    params=dict(Boundary=1,StepAlgorithm=2,Kernel=2,ViscoTreatment=1,Visco=.01,ViscoBoundFactor=1,
                DensityDT=2,DensityDTvalue=.1,Shifting=0,TimeMax=12,TimeOut=.05,
                CoefDtMin=.01,DtIni=0,DtMin=0,RhopOutMin=800,RhopOutMax=1200,MinFluidStop=.999)
    pars=el(ex,'parameters')
    for k,v in params.items():el(pars,'parameter',key=k,value=v)
    dom=el(pars,'simulationdomain');el(dom,'posmin',x=-3.3,y=-1.3,z=-1.3);el(dom,'posmax',x=3.3,y=1.3,z=1.3)
    ET.indent(root);ET.ElementTree(root).write(f/'Case_Def.xml',encoding='utf-8',xml_declaration=True)
    # Velocity and position are integrated from the same tabulated forcing.
    from scipy.integrate import cumulative_trapezoid
    u=brake['initial_speed_m_s']-cumulative_trapezoid(aa,tt,initial=0)
    u=np.maximum(u,0);distance=cumulative_trapezoid(u,tt,initial=0)
    np.savetxt(f/'braking.csv',np.column_stack((tt,-aa,u,distance)),delimiter=',',header='time_s,vehicle_ax_m_s2,vehicle_speed_m_s,vehicle_x_m',comments='')
    setup=dict(case=f.name,layout=case,dp_m=DP,fluid_particles=n,fluid_mass_kg=n*1000*DP**3,
               wall_particles=len(data['walls']),plate_particles=len(data['plates']),duration_s=12.,output_interval_s=.05,
               tank_length_m=6.,tank_inner_radius_m=1.,nominal_fill_fraction=.5,
               actual_nominal_volume_fraction=n*DP**3/(6*np.pi),initial_fill_top_z_m=float(data['fluid'][:,2].max()+DP/2),
               boundary=f'DBC, three-cell outer wall; {round(plate_thickness/DP)}-cell internal plates',
               numerical_plate_thickness_m=plate_thickness,
               plate_note=f'Voxelized plate thickness {plate_thickness*1000:.0f} mm, not realistic sheet-metal thickness. Openings are resolved. Thin-wall force accuracy is not claimed.',
               brake=brake,frame='Fixed tank, +x inertial forcing equals negative prescribed vehicle acceleration. Gravity retained.',
               fluid_initial_cg_m=data['fluid'].mean(0).tolist(),solver_source='Unmodified existing DualSPHysics 5.4 GPU binary',
               limitations=['Single liquid; gas pressure atmospheric, no compressed air or cavitation.',
                            'Artificial viscosity is numerical stabilization, not a resolved turbulent wall model.',
                            'Pilot spatial resolution, not a mesh-converged safety prediction.',
                            'No vehicle motion feedback yet. All four cases use identical prescribed braking.'])
    (f/'setup.json').write_text(json.dumps(setup,indent=2),encoding='utf-8')
    env=os.environ.copy();env['OMP_NUM_THREADS']='4'
    with (f/'gencase.log').open('w') as log:
        p=subprocess.run([str(BIN/'GenCase_win64.exe'),'Case_Def',str(f/'generated/Case'),'-save:all'],cwd=f,env=env,stdout=log,stderr=subprocess.STDOUT)
    if p.returncode:raise RuntimeError(f'GenCase failed: {f}/gencase.log')
    mesh=pv.read(f/'generated/Case_All.vtk');fluid=mesh['Type']==3
    xyz=mesh.points[fluid]
    audit={'passed':int(fluid.sum())==n and np.linalg.norm(xyz.mean(0)[:2])<1e-4,
           'expected_fluid_particles':n,'generated_fluid_particles':int(fluid.sum()),'total_particles':mesh.n_points,
           'generated_fluid_cg_m':xyz.mean(0).tolist(),'distinct_positions':len(np.unique(np.round(mesh.points/DP,4),axis=0)),
           'declared_dp_m':DP}
    audit['passed']=bool(audit['passed'] and audit['distinct_positions']==mesh.n_points)
    (f/'initial_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    if not audit['passed']:raise RuntimeError(f'Particle audit failed: {f}')
    return setup,audit


if __name__=='__main__':
    states,n=lattice();summaries=[]
    for case in CASES:
        setup,audit=prepare(case,states[case],n)
        summaries.append({'case':setup['case'],'mass_kg':setup['fluid_mass_kg'],'particles':audit['total_particles'],'passed':audit['passed']})
        print(json.dumps(summaries[-1]),flush=True)
    (ROOT/'sph/braking_manifest.json').write_text(json.dumps({'cases':summaries,'status':'generated_and_audited_not_run'},indent=2),encoding='utf-8')
