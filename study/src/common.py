"""Shared physical conventions; SI units and exp(-i*omega_encounter*t)."""
from pathlib import Path
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]

def config():return json.loads((ROOT/'case.json').read_text())

def physics(c=None):
    c=c or config();h=c['hull'];e=c['environment']
    L,B,T=h['length_m'],h['beam_m'],h['draft_m']
    volume=4*L*B*T/9;mass=e['rho_kg_m3']*volume
    radii=np.array([h['roll_radius_of_gyration_over_beam']*B,h['pitch_radius_of_gyration_over_length']*L,h['yaw_radius_of_gyration_over_length']*L])
    A=2*L*B/3;Iw=B*L**3/30;zb=-3*T/8;zg=h['center_of_mass_m'][2]
    return dict(volume_m3=volume,mass_kg=mass,inertia_kg_m2=(mass*radii**2).tolist(),
        center_of_buoyancy_z_m=zb,waterplane_area_m2=A,
        hydrostatic_stiffness=[[e['rho_kg_m3']*e['g_m_s2']*A,0],[0,e['rho_kg_m3']*e['g_m_s2']*(Iw+volume*(zb-zg))]],
        froude=c['motion']['forward_speed_m_s']/np.sqrt(e['g_m_s2']*L))

def wave(ratio=None,height=None,c=None):
    c=c or config();w=c['waves'];e=c['environment'];m=c['motion']
    ratio=w['baseline_wavelength_over_length'] if ratio is None else ratio
    H=w['baseline_height_m'] if height is None else height
    lam=ratio*c['hull']['length_m'];k=2*np.pi/lam
    omega=np.sqrt(e['g_m_s2']*k*np.tanh(k*e['water_depth_m']))
    we=omega-k*m['forward_speed_m_s']*np.cos(m['wave_direction_rad'])
    return dict(wavelength_m=lam,height_m=H,amplitude_m=H/2,k_rad_m=k,omega_rad_s=omega,
        encounter_omega_rad_s=we,period_s=2*np.pi/omega,encounter_period_s=2*np.pi/we,
        steepness_H_over_lambda=H/lam,wavelength_over_length=ratio)

def half_width(x,z,c=None):
    h=(c or config())['hull'];L,B,T=h['length_m'],h['beam_m'],h['draft_m']
    return .5*B*np.maximum(0,1-(2*np.asarray(x)/L)**2)*np.maximum(0,1-(np.minimum(np.asarray(z),0)/T)**2)

def hull_quads(nx=48,nz=16,above=False,c=None):
    c=c or config();h=c['hull'];L,T,F=h['length_m'],h['draft_m'],h['freeboard_m']
    xs=np.linspace(-L/2,L/2,nx+1)
    zs=np.linspace(-T,0,nz+1)
    if above:zs=np.r_[zs,np.linspace(0,F,max(3,int(nz*F/T))+1)[1:]]
    vertices=[];faces=[];nrow=len(zs)
    for sign in (-1,1):
        offset=len(vertices)
        for x in xs:
            for z in zs:vertices.append([x,sign*float(half_width(x,z,c)),z])
        for i in range(nx):
            for j in range(nrow-1):
                a=offset+i*nrow+j;b=a+nrow;d=a+1;cc=b+1
                faces.append([a,b,cc,d] if sign==-1 else [a,d,cc,b])
    if above:
        off=(nx+1)*nrow
        for i in range(nx):
            a=i*nrow+nrow-1;b=(i+1)*nrow+nrow-1
            faces.append([a,b,off+b,off+a])
    return np.array(vertices,float),np.array(faces,int)

def airy(x,z,t,w=None,c=None):
    """Head waves in ship frame. Uses intrinsic omega for orbital velocity."""
    c=c or config();w=w or wave(c=c);d=c['environment']['water_depth_m']
    k,om,we,a=w['k_rad_m'],w['omega_rad_s'],w['encounter_omega_rad_s'],w['amplitude_m']
    phase=-k*np.asarray(x)-we*np.asarray(t)
    eta=a*np.cos(phase)
    u=-c['motion']['forward_speed_m_s']-a*om*np.cosh(k*(np.asarray(z)+d))/np.sinh(k*d)*np.cos(phase)
    vz=a*om*np.sinh(k*(np.asarray(z)+d))/np.sinh(k*d)*np.sin(phase)
    return eta,u,vz
