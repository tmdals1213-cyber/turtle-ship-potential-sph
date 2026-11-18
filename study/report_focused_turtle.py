# SPDX-License-Identifier: GPL-3.0-or-later
"""Compact numerical and motion audit for the turtle extreme-wave case."""
from pathlib import Path
import argparse,csv,json,re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('case');a=ap.parse_args()
folder=ROOT/'sph'/a.case;setup=json.loads((folder/'setup.json').read_text(encoding='utf-8'))
motion=json.loads((folder/'focus_motion.json').read_text(encoding='utf-8'));records=motion['records']
rows=[r for r in csv.DictReader((folder/'output/RunPARTs.csv').open(),delimiter=';') if r.get('Part','').isdigit()]
number=lambda r,k:float(r[k].replace(',',''))
t=np.array([r['time_s'] for r in records]);cg=np.array([r['cg_m'] for r in records])
euler=np.array([r['euler_deg'] for r in records]);roof=np.array([r['water_above_roof_particles'] for r in records])
rho=np.array([r['density_range_kg_m3'] for r in records]);cg0=np.array(setup['center_of_mass_m'])
upper=np.array([r['water_above_upper_hull_particles'] for r in records])
log=(folder/'output/Run.out').read_text(encoding='utf-8',errors='replace')
runtime=re.search(r'Simulation Runtime\.+:\s*([0-9.]+)',log,re.IGNORECASE)
steps=re.search(r'Steps of simulation\.+:\s*([0-9,]+)',log,re.IGNORECASE)
stats=dict(case=a.case,status=motion['status'],last_saved_time_s=float(t[-1]),
    initial_particles=int(number(rows[0],'NpSave')),particle_losses=sum(int(number(r,'NpOut')) for r in rows),
    solver_wall_to_last_snapshot_s=number(rows[-1],'SimRuntime [s]'),
    native_steps=int(steps.group(1).replace(',','')) if steps else None,
    sampled_density_range_kg_m3=[float(rho[:,0].min()),float(rho[:,1].max())],
    sampled_max_speed_m_s=max(r['max_speed_m_s'] for r in records),
    max_rigid_fit_rms_m=max(r['fit_rms_m'] for r in records),
    heave_range_m=[float((cg-cg0)[:,2].min()),float((cg-cg0)[:,2].max())],
    pitch_range_deg=[float(euler[:,1].min()),float(euler[:,1].max())],
    surge_range_m=[float((cg-cg0)[:,0].min()),float((cg-cg0)[:,0].max())],
    max_water_above_central_roof_particles=int(roof.max()),
    time_of_max_roof_indicator_s=float(t[roof.argmax()]),
    roof_indicator_definition=motion['roof_indicator'],
    max_water_above_upper_hull_particles=int(upper.max()),
    time_of_max_upper_hull_indicator_s=float(t[upper.argmax()]),
    upper_hull_indicator_definition=motion['upper_hull_indicator'],
    hardware='NVIDIA GeForce RTX 5080',
    limitations=setup['limitations']+['Water above the roof is a geometrical diagnostic; no structural survival or peak-pressure validation.'])
if runtime:stats['native_reported_runtime_s']=float(runtime.group(1))
(folder/'focus_summary.json').write_text(json.dumps(stats,indent=2),encoding='utf-8')
fig,axs=plt.subplots(3,1,figsize=(10,7),sharex=True,facecolor='#080b12')
for ax in axs:
    ax.set_facecolor('#080b12');ax.tick_params(colors='#aebfc8');ax.grid(alpha=.12,color='#aebfc8')
    for sp in ax.spines.values():sp.set_color('#526470')
axs[0].plot(t,cg[:,2]-cg0[2],color='#8dc6d7');axs[0].set_ylabel('Heave [m]',color='white')
axs[1].plot(t,euler[:,1],color='#d4a76c');axs[1].set_ylabel('Pitch [deg]',color='white')
axs[2].plot(t,roof,color='#8dc6d7',label='Above central roof');axs[2].plot(t,upper,color='#d4a76c',label='Near upper surface')
axs[2].set_ylabel('Water above hull\n[particle count]',color='white');axs[2].legend(facecolor='#080b12',labelcolor='white')
axs[2].set_xlabel('Time [s] — model scale',color='white')
fig.suptitle('Measured free-body motion and roof-region water',color='white');fig.tight_layout()
fig.savefig(folder/'focus_motion.png',dpi=150);plt.close(fig)
print(json.dumps(stats,indent=2))
