"""Measure the focused wave actually produced in the empty tank."""
from pathlib import Path
import csv,json,sys
import numpy as np
import pyvista as pv
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
folder=ROOT/'sph'/(sys.argv[1] if len(sys.argv)>1 else 'gunnerus_focus_empty_A085_050_16s')
state=json.loads((folder/'running.json').read_text()); assert state['running'] is False and state['returncode']==0
setup=json.loads((folder/'setup.json').read_text())
rows=[r for r in csv.DictReader((folder/'output/RunPARTs.csv').open(),delimiter=';') if r.get('Part','').isdigit()]
times=np.array([float(r['TimeStep [s]']) for r in rows]); parts=np.array([int(r['Part']) for r in rows])
assert times[-1]>=setup['duration_s']-1e-4
gauges=[]; signals={}
for path in (folder/'output').glob('GaugesSWL_*.csv'):
    q=np.genfromtxt(path,names=True,delimiter=';'); t=q['time_s']; z=q['swlz_m']; x=float(q['swlx_m'][0])
    peak=int(np.argmax(z)); troughs,_=find_peaks(-z,prominence=.01)
    before=troughs[troughs<peak]; after=troughs[troughs>peak]
    local=[int(before[-1])] if len(before) else []
    if len(after): local.append(int(after[0]))
    heights=[float(z[peak]-z[i]) for i in local]
    level=float(np.mean(z[t<.5])); up=np.where((z[:-1]<=level)&(z[1:]>level))[0]+1
    first=up[up<peak]; last=up[up>peak]
    cycle=(int(first[-1]),int(last[0])) if len(first) and len(last) else None
    wave_height=float(np.ptp(z[cycle[0]:cycle[1]])) if cycle else None
    gauges.append(dict(x_m=x,crest_m=float(z[peak]),crest_time_s=float(t[peak]),
        zero_upcrossing_wave_height_m=wave_height,cycle_window_s=None if cycle is None else [float(t[cycle[0]]),float(t[cycle[1]])],
        adjacent_crest_to_trough_heights_m=heights,adjacent_trough_times_s=[float(t[i]) for i in local]))
    signals[x]=(t,z)
near=[q for q in gauges if 6<=q['x_m']<=10]; center=next(q for q in near if q['x_m']==7.5)
peak_time=center['crest_time_s']; plt.style.use('dark_background')
fig,axes=plt.subplots(2,2,figsize=(13,8),facecolor='#080d16')
for x in [12.,10.,7.5,5.]:
    t,z=signals[x]; axes[0,0].plot(t,z,label=f'x={x:g} m')
axes[0,0].set(title='Measured surface elevation',xlabel='Model time [s]',ylabel='m'); axes[0,0].legend(fontsize=8)
sampled=[]
for ax,target in zip([axes[0,1],axes[1,0],axes[1,1]],peak_time+np.array([-.5,0,.5])):
    j=int(np.argmin(abs(times-target))); d=pv.read(folder/f'output/vtks/PartVtk_{parts[j]:04d}.vtk')
    mask=d['Type']==3; xyz=d.points[mask]; speed=np.linalg.norm(d['Vel'][mask],axis=1)
    use=(xyz[:,0]>2)&(xyz[:,0]<13)&(xyz[:,2]>-.7)
    ax.scatter(xyz[use,0],xyz[use,2],c=speed[use],s=3,cmap='winter',vmin=0,vmax=7)
    ax.axvline(7.5,color='#e3b75e',ls='--',lw=1); ax.axhline(0,color='#849cad',lw=.6)
    ax.set(title=f'Actual water particles | t={times[j]:.2f} s',xlim=(2,13),ylim=(-.7,1.5),xlabel='x [m]',ylabel='z [m]')
    sampled.append(dict(time_s=float(times[j]),max_speed_m_s=float(speed.max()),density_range_kg_m3=[float(d['Rho'][mask].min()),float(d['Rho'][mask].max())]))
for ax in axes.flat:ax.grid(alpha=.12)
fig.suptitle('FOCUSED WAVE | EMPTY-TANK COMMISSIONING',fontsize=19)
figure_name=f'{folder.name}_wave_check.png' if len(sys.argv)>1 else 'gunnerus_focused_wave_check.png'
fig.tight_layout(rect=[0,.02,1,.95]); fig.savefig(ROOT/'media'/figure_name,dpi=140); plt.close(fig)
losses=sum(int(r['NpOut'].replace(',','')) for r in rows)
clamps=sum(int(r['DTsMin'].replace(',','')) for r in rows)
report=dict(status='complete',case=folder.name,particle_losses=losses,dt_clamp_events=clamps,
    nominal_focus_station=center,near_encounter_stations=near,sampled_flow=sampled,
    definition='Crest height is above still water. Wave height uses adjacent measured troughs, not global extrema or significant wave height.',
    interpretation='Empty tank only. Actual breaking, deck water and downward vessel force require the ship calculation.',
    commissioning_passed=bool(losses==0 and clamps==0 and center['crest_m']>.5 and max(center['adjacent_crest_to_trough_heights_m'])>.8))
(folder/'wave_report.json').write_text(json.dumps(report,indent=2)); print(json.dumps(report,indent=2))
