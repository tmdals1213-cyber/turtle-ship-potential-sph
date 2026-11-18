# SPDX-License-Identifier: GPL-3.0-or-later
"""Prepend a verified computed SPH hook to the accepted portrait film."""
from pathlib import Path
import argparse,hashlib,json,re,subprocess

ROOT=Path(__file__).resolve().parent;MEDIA=ROOT/'media'
ap=argparse.ArgumentParser();ap.add_argument('--intro',default='turtle_focused_intro');a=ap.parse_args()
intro=MEDIA/f'{a.intro}.mp4';base=MEDIA/'turtle_ship_story_v6_portrait.mp4'
im=json.loads(intro.with_suffix('.json').read_text(encoding='utf-8'))
bm=json.loads(base.with_suffix('.json').read_text(encoding='utf-8'))
assert hashlib.sha256(intro.read_bytes()).hexdigest()==im['sha256']
path=MEDIA/'turtle_ship_story_v7_focused_portrait.mp4'
subprocess.run(['ffmpeg','-v','error','-y','-i',str(intro),'-i',str(base),'-filter_complex',
    '[0:v]setpts=PTS-STARTPTS[v0];[1:v]setpts=PTS-STARTPTS[v1];[v0][v1]concat=n=2:v=1:a=0[v]',
    '-map','[v]','-an','-c:v','libx264','-preset','fast','-crf','19','-pix_fmt','yuv420p','-r','30',
    '-movflags','+faststart',str(path)],check=True)
decoded=subprocess.run(['ffmpeg','-v','error','-nostats','-progress','pipe:1','-i',str(path),'-f','null','-'],capture_output=True,text=True,check=True)
frames=int(re.findall(r'^frame=(\d+)$',decoded.stdout,re.M)[-1])
expected=round(30*(im['video_seconds']+bm['duration_s']));assert frames==expected,(frames,expected)
assert not decoded.stderr.strip(),decoded.stderr
report=dict(duration_s=im['video_seconds']+bm['duration_s'],size=[1080,1920],fps=30,decoded_frames=frames,
    full_decode='passed_without_errors',intro_source=im,base_source=str(base),
    source_simulations=[im['case'],bm['source_simulation']],
    sections=[['computed extreme-wave hook',0,im['video_seconds']]]+
        [[name,start+im['video_seconds'],end+im['video_seconds']] for name,start,end in bm['sections']],
    conditions='Intro: initially resting ship in focused extreme waves. Main SPH result: initial forward speed and regular waves. Same physical hull and mass, separate simulations.',
    sha256=hashlib.sha256(path.read_bytes()).hexdigest())
path.with_suffix('.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
summary=json.loads((ROOT/'sph'/im['case']/'focus_summary.json').read_text(encoding='utf-8'))
description=(ROOT/'templates/turtle_ship_story_v6_description.txt').read_text(encoding='utf-8')
lead=(
    'V7: computed extreme-wave intro, followed by the accepted V6 portrait story.\n\n'
    f"The first {im['video_seconds']:g} video seconds show the SAME turtle hull in a SEPARATE SPH condition: "
    'initially at rest, free six-DOF motion, 24 x 4.8 m tank, 1.2 m water depth, 40 mm particles. '
    'A first-order JONSWAP/NewWave piston signal produces a nonlinear wave group. '
    'Empty-tank commissioning was completed before adding the ship. '
    'The linear target focus time/amplitude are forcing parameters, not measured breaking-wave properties.\n\n'
    f"Intro data interval: {im['physical_interval_s'][0]:.2f} to {im['physical_interval_s'][1]:.2f} model seconds. "
    'Rigid motion comes from the computed body particles. Water follows the same fluid-particle IDs between outputs; '
    'a smoothed three-dimensional surface is reconstructed for rendering. No analytic crash surface or prescribed boat motion is added.\n\n'
    f"Computed case: {summary['last_saved_time_s']:.3f} model seconds, {summary['initial_particles']:,} particles, "
    f"{summary['particle_losses']} excluded particles, approximately {summary['solver_wall_to_last_snapshot_s']/60:.1f} minutes on RTX 5080. "
    'This is an extreme-wave visualization demonstration, not a structural survivability or peak-impact-pressure prediction. '
    'The closed upper hull interacts with water; air entrapment, internal flooding, damage and decorative appendage forces are excluded.\n\n'
    'Wave-generation reference:\nhttps://github.com/DualSPHysics/DualSPHysics/wiki/3.-SPH-formulation\n\n'
    'Main film and model credit:\n\n'
)
path.with_name(path.stem+'_description.txt').write_text(lead+description,encoding='utf-8')
subprocess.run(['ffmpeg','-v','error','-y','-ss',str(im['video_seconds']/2),'-i',str(path),'-frames:v','1','-update','1',
    str(MEDIA/'turtle_v7_encoded_intro_check.png')],check=True)
print(json.dumps(report,indent=2))
