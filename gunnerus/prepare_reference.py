"""Stage portable reference inputs; never launch a solver or overwrite a run."""
from pathlib import Path
import argparse
import json
import shutil

REFERENCE = Path(__file__).resolve().parent.parent / 'SPH_CFD-incomp-free-surface-6DOF_GPU_Cpp-CUDA-Python_05-gunnerus-focused-wave'
CASES = [('gunnerus_focus_A115_x5.5_D300_050_16s', 'setup'),
         ('gunnerus_focus_empty_A115_D300_050_16s', 'empty_flume/setup')]


def stage(workspace):
    workspace = workspace.resolve()
    targets = [workspace/'geometry'] + [workspace/'sph'/name for name, _ in CASES]
    for target in targets:
        if target.exists():
            raise FileExistsError(f'Preserve existing inputs: {target}')
    workspace.mkdir(parents=True, exist_ok=True)
    for directory in ['sph', 'media', 'results']:
        (workspace/directory).mkdir(exist_ok=True)
    shutil.copytree(REFERENCE/'geometry', workspace/'geometry')
    for path in (workspace/'geometry').rglob('*.json'):
        path.write_text(path.read_text(encoding='utf-8').replace('@GUNNERUS_ROOT@', workspace.as_posix()), encoding='utf-8')
    for name, source in CASES:
        target = workspace/'sph'/name
        shutil.copytree(REFERENCE/source, target)
        for path in target.iterdir():
            if path.suffix in {'.json', '.xml'}:
                path.write_text(path.read_text(encoding='utf-8').replace('@GUNNERUS_ROOT@', workspace.as_posix()), encoding='utf-8')
        settings = json.loads((target/'setup.json').read_text(encoding='utf-8'))
        settings['status'] = 'prepared_not_run'
        (target/'setup.json').write_text(json.dumps(settings, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'workspace': str(workspace), 'cases': [name for name, _ in CASES], 'solver_launched': False}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path, default=Path(__file__).resolve().parent)
    stage(parser.parse_args().workspace)
