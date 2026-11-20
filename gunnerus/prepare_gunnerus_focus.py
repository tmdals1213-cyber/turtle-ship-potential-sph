"""Larger-domain focused crest, using the existing native piston and six-DOF solver."""
import argparse,json,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from prepare_gunnerus_waves import prepare,ROOT,el


def main(empty,amplitude,station,period=2.6,depth=2.):
    period_tag='' if period==2.6 else f'_T{round(period*100):03d}'
    depth_tag='' if depth==2. else f'_D{round(depth*100):03d}'
    name=f'gunnerus_focus_{"empty_" if empty else ""}A{round(amplitude*100):03d}{"" if empty else f"_x{station:g}"}{period_tag}{depth_tag}_050_16s'
    duration=16.; piston=18.; focus=7.5; focus_time=12.
    from place_gunnerus_contact import at_station
    physical=None if empty else at_station(station)
    prepare('empty' if empty else 'waves',name=name,depth=depth,duration=duration,
            domain=[-18.,18.,-5.,5.,-depth,3.],physical_model=physical)
    folder=ROOT/'sph'/name
    freq=np.linspace(.24,.85,101); omega=2*np.pi*freq; fp=1/period
    sigma=np.where(freq<=fp,.07,.09)
    spec=freq**-5*np.exp(-1.25*(fp/freq)**4)*3.3**np.exp(-.5*((freq-fp)/(sigma*fp))**2)
    amp=amplitude*spec/spec.sum()
    k=np.array([brentq(lambda q:9.81*q*np.tanh(q*depth)-w*w,.001,20) for w in omega])
    transfer=4*np.sinh(k*depth)**2/(np.sinh(2*k*depth)+2*k*depth)
    phase=k*(piston-focus)-omega*focus_time
    t=np.arange(0,duration+1.001,.002)
    ramp=.5*(1-np.cos(np.pi*np.minimum(t/2,1)))
    end=.5*(1-np.cos(np.pi*np.minimum(np.maximum(duration+1-t,0)/2,1)))
    disp=-np.sum(amp/transfer*np.sin(t[:,None]*omega+phase),axis=1)*ramp*end
    assert abs(disp).max()<.85, 'Piston stroke must remain above the extended floor'
    np.savetxt(folder/'piston_motion.dat',np.c_[t,disp],fmt='%.10g')
    np.savez(folder/'wave_components.npz',frequency=freq,omega=omega,k=k,amplitude=amp,transfer=transfer,phase=phase,time=t,displacement=disp)
    tree=ET.parse(folder/'Case_Def.xml'); root=tree.getroot()
    if not empty: root.find('.//fillpoint').set('x',str(station))
    root.find('.//constantsdef/speedsound').set('value','80')
    root.find('.//constantsdef/speedsystem').set('value','4')
    root.find('.//parameter[@key="TimeOutExtra"]').set('value','.005')
    gauges=root.find('.//special/gauges')
    for child in list(gauges):
        if child.tag!='default': gauges.remove(child)
    for x in [16.,14.,12.,10.,9.,8.,7.5,7.,6.,5.,4.,2.,0.]:
        gu=el(gauges,'swl',name='x'+str(x).replace('.','p'))
        el(gu,'pointdp',coefdp=.25)
        el(gu,'point0',x=x,y=0 if empty else 3.,z=-1.)
        el(gu,'point2',x=x,y=0 if empty else 3.,z=2.)
    damp=root.find('.//dampingzone')
    damp.find('limitmin').set('x','-11')
    damp.find('limitmax').set('x','-18')
    ET.indent(root); tree.write(folder/'Case_Def.xml',encoding='utf-8',xml_declaration=True)
    s=json.loads((folder/'setup.json').read_text())
    for key in ['target_H_m','wavelength_m','wave_ramp_s']: s.pop(key,None)
    s.update(kind='focused_empty' if empty else 'focused_ship',intrinsic_period_s=period,
        target_linear_crest_amplitude_m=amplitude,focus_x_m=focus,focus_time_s=focus_time,
        piston_x_m=piston,frequency_band_hz=[float(freq[0]),float(freq[-1])],
        wave_input='First-order crest-focused JONSWAP/NewWave piston; measured SPH height and breaking are commissioning outputs.',
        piston_displacement_range_m=[float(disp.min()),float(disp.max())],piston_max_speed_m_s=float(abs(np.gradient(disp,t)).max()),
        initial_station_x_m=0 if empty else station,
        focus_basis='Initial station aligned with the breaking-front region measured in the empty tank; actual ship motion remains free.',
        interpretation='Linear target crest above still water is not crest-to-trough wave height or significant wave height.',
        sound_speed_m_s=80,water_only=True,scalar_output_interval_s=.005)
    s['output']['native_force_motion_extra_interval_s']=.005
    s['limitations']=[x for x in s['limitations'] if not x.startswith('10-second')]
    s['limitations'].append('16-second transient with a focused wave packet; no extra downward body force is prescribed.')
    (folder/'setup.json').write_text(json.dumps(s,indent=2)); print(json.dumps({'case':name,'target_crest_m':amplitude,'duration_s':duration,'max_piston_speed_m_s':s['piston_max_speed_m_s']},indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--empty',action='store_true'); ap.add_argument('--amplitude',type=float,default=.85)
    ap.add_argument('--station-x',type=float,default=2.)
    ap.add_argument('--period',type=float,default=2.6)
    ap.add_argument('--depth',type=float,default=2.)
    a=ap.parse_args(); main(a.empty,a.amplitude,a.station_x,a.period,a.depth)
