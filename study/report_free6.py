# SPDX-License-Identifier: GPL-3.0-or-later
"""Plot all measured rigid DOFs and report short-window wave gauges honestly."""
from pathlib import Path
import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent;folder=ROOT/'sph/free6_history_40'
m=json.loads((folder/'motion_6dof.json').read_text(encoding='utf-8'));s=json.loads((folder/'setup.json').read_text(encoding='utf-8'))
t=np.array(m['times_s']);p=np.array(m['cg_m'])-s['center_of_mass_m'];r=np.array(m['euler_xyz_deg']);vel=np.array(m['velocity_m_s'])
labels=['Surge displacement [m]','Sway [cm]','Heave [cm]','Roll [deg]','Pitch [deg]','Yaw [deg]']
yy=np.c_[p[:,0],p[:,1:]*100,r];plt.style.use('dark_background')
fig,axs=plt.subplots(2,3,figsize=(13,6.7),layout='constrained',facecolor='#080b12')
for i,ax in enumerate(axs.ravel()):
    ax.set_facecolor('#0b1320');ax.plot(t,yy[:,i],color='#83d7ef',lw=1.7);ax.axhline(0,color='#6b8299',lw=.6)
    ax.set_title(labels[i],loc='left',fontsize=12);ax.grid(alpha=.13);ax.set_xlabel('Model time [s]');ax.set_xlim(0,max(t))
    if i==1:lim=max(1,np.max(abs(yy[:,i]))*1.2);ax.set_ylim(-lim,lim)
    if i in [3,5]:lim=max(1,np.max(abs(yy[:,i]))*1.2);ax.set_ylim(-lim,lim)
fig.suptitle('Initial speed only — six free rigid-body DOFs\nIndependent SPH demonstration; 1:10 model, no propulsion',fontsize=16)
fig.savefig(ROOT/'media/turtle_history_six_dof.png',dpi=150);plt.close(fig)
gauges=[];omega=np.sqrt(9.81*(2*np.pi/3)*np.tanh(2*np.pi/3))
for path in sorted((folder/'output').glob('GaugesSWL_*.csv')):
    rows=list(csv.DictReader(path.open(),delimiter=';'));tt=np.array([float(q['time [s]']) for q in rows]);eta=np.array([float(q['swlz [m]']) for q in rows])
    window=tt>=max(0,tt[-1]-4.2);x=tt[window];z=eta[window]
    A=np.c_[np.ones(len(x)),x-x.mean(),np.cos(omega*x),np.sin(omega*x)];q=np.linalg.lstsq(A,z,rcond=None)[0]
    gauges.append(dict(gauge=path.stem,window_s=[float(x[0]),float(x[-1])],fundamental_H_m=float(2*np.hypot(q[2],q[3])),raw_crest_trough_m=float(np.ptp(z)),fit_rms_m=float(np.sqrt(np.mean((A@q-z)**2)))))
summary=dict(case=s.get('hull'),last_saved_time_s=float(t[-1]),initial_velocity_m_s=s['initial_velocity_m_s'],last_estimated_velocity_m_s=vel[-1].tolist(),
    displacement_ranges_m=[{'min':float(p[:,i].min()),'max':float(p[:,i].max())} for i in range(3)],
    angle_ranges_deg=[{'min':float(r[:,i].min()),'max':float(r[:,i].max())} for i in range(3)],
    wave_gauges=gauges,wave_fit_limitation='Short transient window. Total elevation includes body scattering and tank reflection; not incident/reflected-wave separation.',
    status='completed_demonstration_not_validation' if (folder/'summary.json').exists() else 'partial_run_diagnostic',
    comparison_limitation='BEM response is at fixed mean reference speed; SPH coasts. Differences are not a controlled method-accuracy comparison.')
(folder/'response_report.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
