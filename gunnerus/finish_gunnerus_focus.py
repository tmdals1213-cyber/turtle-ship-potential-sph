"""Follow the existing focused-wave job and finish its declared visualization.

No new simulation, retuning, kernel change or automatic debugging is performed.
"""
from pathlib import Path
import json,time,traceback,subprocess,sys,re,argparse
from report_gunnerus_waves import process
from sph_segments import output_dirs

ROOT=Path(__file__).resolve().parent; NAME='gunnerus_focus_A085_x2_050_16s'; FOLDER=ROOT/'sph'/NAME
STATUS=ROOT/'results/gunnerus_focus_status.json'
VIDEO_NAME='gunnerus_focused_waves_story_portrait';PREFIX='gunnerus_focus';DILUTE=False;SURFACE=.5


def save(status,**kwargs):
    value=dict(case=NAME,status=status,updated_epoch_s=time.time(),**kwargs)
    temp=STATUS.with_suffix('.tmp');temp.write_text(json.dumps(value,indent=2));temp.replace(STATUS)


def main():
    while True:
        try:state=json.loads((FOLDER/'running.json').read_text())
        except (json.JSONDecodeError,OSError):time.sleep(10);continue
        if (output_dirs(FOLDER)[-1]/'Run.out').exists() and state.get('running') is False:break
        save('solver_running',solver_state=state);time.sleep(10)
    if state.get('returncode')!=0 or state.get('stop_reason'):raise RuntimeError(str(state))
    save('native_load_postprocessing',solver_state=state)
    summary=process(NAME)
    save('rendering_story',simulation_summary=summary)
    command=[sys.executable,str(ROOT/'render_gunnerus_focus_story.py'),'--case',NAME,'--prefix',PREFIX,'--video-name',VIDEO_NAME,'--surface-level',str(SURFACE)]
    if DILUTE:command.append('--dilute-water')
    subprocess.run(command,check=True)
    video=ROOT/'media'/f'{VIDEO_NAME}.mp4'
    result=subprocess.run(['ffmpeg','-v','error','-nostats','-progress','pipe:1','-i',str(video),'-f','null','-'],capture_output=True,text=True,check=True)
    frames=[int(q) for q in re.findall(r'^frame=(\d+)',result.stdout,re.M)];assert frames[-1]==660
    saved=ROOT/'media'/f'{PREFIX}_encoded_check.png'
    subprocess.run(['ffmpeg','-v','error','-y','-ss','11','-i',str(video),'-frames:v','1','-update','1',str(saved)],check=True)
    save('complete_pending_final_visual_review',simulation_summary=summary,video=str(video),decoded_frames=frames[-1],
        note='Full video decoded; inspect encoded check image before presenting the final new result.')
    print('Focused-wave calculation and story rendering finished.',flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--case',default=NAME);ap.add_argument('--prefix',default=PREFIX)
    ap.add_argument('--video-name',default=VIDEO_NAME);ap.add_argument('--dilute-water',action='store_true');ap.add_argument('--surface-level',type=float,default=SURFACE)
    a=ap.parse_args();NAME=a.case;FOLDER=ROOT/'sph'/NAME;PREFIX=a.prefix;VIDEO_NAME=a.video_name;DILUTE=a.dilute_water;SURFACE=a.surface_level
    if NAME!='gunnerus_focus_A085_x2_050_16s':STATUS=ROOT/'results'/f'{NAME}_status.json'
    try:main()
    except Exception:
        error=traceback.format_exc();save('needs_attention',error=error);print(error,flush=True);raise
