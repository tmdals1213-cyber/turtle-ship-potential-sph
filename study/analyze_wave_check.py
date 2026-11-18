# SPDX-License-Identifier: GPL-3.0-or-later
"""Fit actual gauges, retaining spatial attenuation and the specified fit window."""
from pathlib import Path
import argparse,json,csv
import numpy as np
ROOT=Path(__file__).resolve().parent
ap=argparse.ArgumentParser();ap.add_argument('case',default='piston_2d',nargs='?');ap.add_argument('--start',type=float,default=9);a=ap.parse_args()
p=ROOT/'sph'/a.case;setup=json.loads((p/'setup.json').read_text(encoding='utf-8'));om=2*np.pi/setup['intrinsic_period_s'];rows=[]
for f in sorted((p/'output').glob('GaugesSWL_*.csv')):
    q=np.genfromtxt(f,names=True,delimiter=';');t=q['time_s'];z=q['swlz_m'];j=t>=a.start
    mat=np.c_[np.cos(om*t[j]),np.sin(om*t[j]),np.ones(sum(j)),t[j]-t[j].mean()]
    c=np.linalg.lstsq(mat,z[j],rcond=None)[0];r=z[j]-mat@c;H=2*np.hypot(*c[:2])
    rows.append(dict(x_m=float(q['swlx_m'][0]),fit_window_s=[a.start,float(t[-1])],H_m=float(H),H_over_target=float(H/setup['target_H_m']),
        mean_surface_m=float(c[2]),mean_drift_m_s=float(c[3]),residual_rms_m=float(np.sqrt(np.mean(r*r))),
        phase_deg=float(np.degrees(np.arctan2(-c[1],c[0])))))
result=dict(case=a.case,fit='Fixed prescribed frequency, constant and linear drift; does not separate reflected waves.',
    target_H_m=setup['target_H_m'],rows=rows,status='wave_maker_check_not_ship_validation')
(p/'wave_fit.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
