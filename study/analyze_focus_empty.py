# SPDX-License-Identifier: GPL-3.0-or-later
"""Read actual SPH empty-flume results; do not substitute the target wave."""
from pathlib import Path
import argparse,csv,json
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('case');a=ap.parse_args()
folder=ROOT/'sph'/a.case
setup=json.loads((folder/'setup.json').read_text(encoding='utf-8'))
rows=[r for r in csv.DictReader((folder/'output/RunPARTs.csv').open(),delimiter=';') if r.get('Part','').isdigit()]
times=np.array([float(r['TimeStep [s]']) for r in rows]);parts=np.array([int(r['Part']) for r in rows])
gauges={}
for path in sorted((folder/'output').glob('GaugesSWL_*.csv')):
    lines=path.read_text(encoding='utf-8').splitlines()
    header=next(i for i,s in enumerate(lines) if s.startswith('time'))
    data=np.loadtxt(path,delimiter=';',skiprows=header+1)
    # Native columns: time; elevation is column z (fourth column).
    columns=lines[header].split(';')
    zcol=columns.index('swlz [m]')
    z=data[:,zcol];tt=data[:,0]
    if float(data[0,columns.index('pos0x [m]')])>=setup.get('piston_x_m',12)-.8:continue
    gauges[path.stem]=dict(max_elevation_m=float(z.max()),time_at_max_s=float(tt[np.argmax(z)]),min_elevation_m=float(z.min()))
maxentry=max(gauges.values(),key=lambda q:q['max_elevation_m'])
focus_key='GaugesSWL_swl_'+str(float(setup['focus_x_m'])).replace('-','m').replace('.','p')
center=gauges.get(focus_key,maxentry)['time_at_max_s']
fig,axs=plt.subplots(3,2,figsize=(14,8),sharex=True,sharey=True,facecolor='#080b12')
rho_min=1e9;rho_max=0;speed_max=0
for ax,tt in zip(axs.flat,center+np.linspace(-1,1.5,6)):
    j=np.argmin(abs(times-tt));m=pv.read(folder/f'output/vtks/PartVtk_{parts[j]:04d}.vtk')
    fluid=m['Type']==3;p=m.points[fluid];r=m['Rho'][fluid];v=np.linalg.norm(m['Vel'][fluid],axis=1)
    rho_min=min(rho_min,float(r.min()));rho_max=max(rho_max,float(r.max()));speed_max=max(speed_max,float(v.max()))
    select=(p[:,0]>-2)&(p[:,0]<7)&(p[:,2]>-.6)
    ax.set_facecolor('#080b12');ax.scatter(p[select,0],p[select,2],c=v[select],s=3,cmap='winter',vmin=0,vmax=3)
    ax.axvline(setup['focus_x_m'],color='#efad65',lw=1,ls='--')
    ax.axhline(0,color='#75838c',lw=.5)
    ax.set(xlim=(-2,7),ylim=(-.6,1.15),title=f'Actual SPH particles | t = {times[j]:.2f} s')
    ax.title.set_color('white');ax.tick_params(colors='#aebfc8')
    for sp in ax.spines.values():sp.set_color('#526470')
fig.suptitle('Empty tank: focused-wave commissioning (model scale)',color='white',fontsize=16)
fig.tight_layout();fig.savefig(folder/'wave_profiles.png',dpi=140);plt.close(fig)
report=dict(case=a.case,status='empty_flume_completed' if times[-1]>=setup['duration_s']-.09 else 'partial',
    last_saved_time_s=float(times[-1]),particle_losses=sum(int(float(r['NpOut'])) for r in rows),gauges=gauges,
    sampled_density_range_kg_m3=[rho_min,rho_max],sampled_max_speed_m_s=speed_max,
    note='Gauge elevation is not significant wave height. Inspect particle profiles for actual overturning; no automated breaking classification.')
(folder/'wave_report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
