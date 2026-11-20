"""Generate and audit the four accepted input layouts without running CFD."""
from pathlib import Path
from prepare_braking import lattice,prepare as prepare_particles
from prepare_coupled_vehicle import prepare as prepare_vehicle
from render_intro_motion import export_geometry

def prepare_all():
    export_geometry()
    initial,n=lattice()
    assert n==137812,'Reference mass changed'
    prepare_particles('A_unpartitioned',initial['A_unpartitioned'],n)
    thick,n=lattice(plate_halfwidth=.081,fluid_count=n)
    order=['A_unpartitioned','C_separate_compartments','B_perforated_baffles','D_cross_baffles']
    for name in order[1:]:prepare_particles(name,thick[name],n,suffix='thick_040_v2',plate_thickness=.16)
    folders=[]
    for name in order:
        folder=prepare_vehicle(name,.69,duration=8,version='v8',
            source_suffix='brake_040_v1' if name==order[0] else 'thick_040_v2',
            lateral_start=1.5,sound_speed=220,cfl=.2)
        folders.append(folder)
    return folders

if __name__=='__main__':
    for folder in prepare_all():print(folder)
