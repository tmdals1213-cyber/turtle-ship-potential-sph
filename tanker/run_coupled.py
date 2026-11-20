"""Preserving runner for native two-way SPH/Chrono tanker cases."""
from pathlib import Path
import argparse, json, os, subprocess, time
from run_braking import rows

ROOT = Path(__file__).resolve().parent
from tanker_paths import SOLVER
BIN = SOLVER.parent


def run(name, smoke=None):
    folder = ROOT / 'sph' / name
    suffix = '_smoke' if smoke is not None else ''
    output = folder / ('output' + suffix)
    result_path = folder / ('run' + suffix + '.json')
    if result_path.exists():
        return json.loads(result_path.read_text())
    if output.exists():
        raise FileExistsError(f'Preserve existing output: {output}')
    assert json.loads((folder / 'initial_audit.json').read_text())['passed']
    setup = json.loads((folder / 'setup.json').read_text())
    duration = float(smoke if smoke is not None else setup['duration_s'])
    exe = SOLVER
    assert exe.exists(), 'Chrono build must finish first'
    env = os.environ.copy()
    env['PATH'] = str(BIN) + os.pathsep + env.get('PATH', '')
    env['OMP_NUM_THREADS'] = '4'
    command = [str(exe), '-gpu:0', str(folder/'generated/Case'), str(output),
               '-dirdataout', 'data', '-sv:binx,info', '-svres', '-toutx:0.01']
    if smoke is not None:
        command += [f'-tmax:{duration}']
    start = time.perf_counter()
    reason = None
    with (folder / ('solver'+suffix+'.log')).open('w') as log:
        proc = subprocess.Popen(command, cwd=folder, env=env, stdout=log, stderr=subprocess.STDOUT)
        while proc.poll() is None:
            rr = rows(output / 'RunPARTs.csv')
            loss = sum(float(r['NpOut']) for r in rr)
            if loss > setup['fluid_particles'] * .001:
                reason = 'More than 0.1 percent particles excluded'
            if len(rr) >= 4 and all(float(r['DTsMin']) > 10 for r in rr[-3:]):
                reason = 'Sustained minimum time-step clamping'
            state = dict(case=name, running=True, pid=proc.pid, output=str(output),
                         wall_s=time.perf_counter()-start,
                         physical_s=float(rr[-1]['TimeStep [s]']) if rr else 0,
                         particle_losses=loss, stop_reason=reason)
            (folder/'running.json').write_text(json.dumps(state, indent=2))
            if reason:
                proc.terminate()
                break
            time.sleep(1)
        code = proc.wait()
    rr = rows(output / 'RunPARTs.csv')
    last = float(rr[-1]['TimeStep [s]']) if rr else 0
    native = (output/'Run.out').read_text(errors='replace') if (output/'Run.out').exists() else ''
    complete = code == 0 and reason is None and last >= duration-1e-6 and 'Finished execution (code=0)' in native
    result = dict(case=name, status='complete' if complete else 'failed_or_incomplete',
                  returncode=code, stop_reason=reason, last_saved_time_s=last,
                  wall_s=time.perf_counter()-start,
                  particle_losses=sum(float(r['NpOut']) for r in rr),
                  dt_min_clamps=sum(float(r['DTsMin']) for r in rr),
                  min_dt_s=min((float(r['DtMin [s]']) for r in rr if float(r['DtMin [s]']) > 0), default=None),
                  solver=str(exe), output=str(output), command=command)
    result_path.write_text(json.dumps(result, indent=2))
    (folder/'running.json').write_text(json.dumps(dict(running=False, **result), indent=2))
    print(json.dumps(result), flush=True)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('case')
    p.add_argument('--smoke', type=float)
    a = p.parse_args()
    run(a.case, a.smoke)
