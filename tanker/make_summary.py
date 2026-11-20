"""Rebuild compact metrics from regenerated native fields."""
from pathlib import Path
import csv,hashlib,json,xml.etree.ElementTree as ET
import numpy as np
from render_coupled_results import LABELS
ROOT=Path(__file__).resolve().parent

def summarize():
    verified=json.loads((ROOT/'sph/verified_comparison.json').read_text())
    assert len(verified['cases'])==4 and all(c['qualified'] for c in verified['cases'])
    records=[];fingerprints=[];setups=[]
    for case in verified['cases']:
        folder=ROOT/'sph'/case['case'];s=json.loads((folder/'setup.json').read_text());setups.append(s)
        xml=ET.parse(folder/'Case_Def.xml').getroot()
        floating=xml.find('.//floatings/floating')
        fingerprints.append(dict(case=case['case'],
            input_csv=hashlib.sha256((folder/'input_acceleration.csv').read_bytes()).hexdigest(),
            floating=hashlib.sha256(ET.tostring(floating)).hexdigest(),
            contact_mesh=hashlib.sha256((folder/'truck_contact.vtk').read_bytes()).hexdigest(),
            road_mesh=hashlib.sha256((folder/'road_contact.vtk').read_bytes()).hexdigest()))
        for tag in ['translationDOF','rotationDOF']:
            assert all(floating.find(tag).attrib[k]=='1' for k in 'xyz')
        check=json.loads((folder/'surface_loads/analysis.json').read_text())
        assert check['qualified_for_load_visualization']
        data=np.load(folder/'surface_loads/loads.npz');time=data['time_s']
        use=(time>=1.5)&(time<=3.5001)
        force=data['force_body_filtered_N'].sum(1);moment=data['moment_body_filtered_Nm'].sum(1)
        motion=case['motion']
        records.append(dict(case=case['case'],layout=s['layout'],label=LABELS[s['layout']],
            first_side_lift_s=motion['first_side_lift_time_s'],first_75deg_s=motion['first_75deg_time_s'],
            overturned=motion['overturned'],max_tilt_deg=motion['max_tilt_deg'],
            peak_total_fluid_force_kN=float(np.linalg.norm(force[use],axis=1).max()/1000),
            peak_total_fluid_moment_kNm=float(np.linalg.norm(moment[use],axis=1).max()/1000),
            peak_side_force_kN=float(abs(force[use,1]).max()/1000),
            peak_roll_moment_kNm=float(abs(moment[use,0]).max()/1000),
            local_particle_load_max_N=float(data['per_particle_force_percentiles_N'][:,-1].max()),
            force_check_nrmse=check['native_force_comparison']['filtered_nrmse'],
            moment_check_nrmse=check['native_moment_comparison']['filtered_nrmse'],
            wall_s=case['run']['wall_s'],particle_loss=case['run']['particle_losses']))
    for key in ['input_csv','floating','contact_mesh','road_mesh']:
        assert len({f[key] for f in fingerprints})==1,f'Unequal {key}'
    common_keys=['fluid_mass_kg','dry_mass_kg','dry_inertia_kg_m2','dry_cg_m','contact_mu',
                 'dp_m','sound_speed_m_s','cfl_number','duration_s']
    for key in common_keys:assert all(s[key]==setups[0][key] for s in setups),key
    result=dict(cases=records,input_verification=fingerprints,identical_parameters=common_keys,
        comparison_window_s=[1.5,3.5],load_filter='Cubic Savitzky-Golay across nine 0.05s snapshots (0.4s span).',
        moment_reference='Instantaneous dry-vehicle CG; roll component is body x. Loads include fluid pressure and viscosity; exclude dry-body weight and ground contact.',
        rollover_definition=verified['cases'][0]['motion']['overturn_definition'],
        observation_duration_s=8,
        limitations=['Generic rigid vehicle with fixed wheels, frictional road contact, no suspension, tire handling or structural damage.',
          'Prescribed nonrotating support-frame acceleration; visual road speed is the prescribed input, not solved total sliding speed.',
          'First spatial resolution: DP40mm; four-layer160mm baffle envelopes. Not a calibrated real-vehicle safety threshold.',
          'Same dry mass/inertia across layouts. Added baffle mass not separately modeled; equal water mass produces slightly different initial levels/CGs.',
          'Single liquid with atmospheric gas pressure. Post-roll impacts are illustrative, not validated local stress predictions.'])
    (ROOT/'comparison_summary.json').write_text(json.dumps(result,indent=2))
    with (ROOT/'comparison_summary.csv').open('w',newline='') as out:
        writer=csv.DictWriter(out,fieldnames=records[0].keys());writer.writeheader();writer.writerows(records)
    return result

if __name__=="__main__":summarize()
