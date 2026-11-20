"""Read preserved native output segments as one physical-time history."""
import csv,json
from pathlib import Path
import numpy as np
from read_native_float import read_float


def output_dirs(folder):
    manifest=Path(folder)/'output_segments.json'
    names=json.loads(manifest.read_text())['outputs'] if manifest.exists() else ['output']
    return [Path(folder)/name for name in names]


def part_rows(folder):
    parts={}
    for directory in output_dirs(folder):
        with (directory/'RunPARTs.csv').open() as stream:
            rows=[r for r in csv.DictReader(stream,delimiter=';') if r.get('Part','').isdigit()]
        for row in rows:
            part=int(row['Part'])
            if part in parts:
                assert abs(float(row['TimeStep [s]'])-float(parts[part]['TimeStep [s]']))<1e-6
                continue  # Keep the original checkpoint sample, not the restart placeholder.
            parts[part]={**row,'output_directory':str(directory)}
    result=[parts[p] for p in sorted(parts)]
    assert all(float(b['TimeStep [s]'])>float(a['TimeStep [s]']) for a,b in zip(result,result[1:]))
    return result


def vtk_path(row):
    return Path(row['output_directory'])/f"vtks/PartVtk_{int(row['Part']):04d}.vtk"


def native_history(folder,extra=True):
    chunks=[];metadata=[];last=-np.inf
    filename='PartFloatInfo2.ibi4' if extra else 'PartFloatInfo.ibi4'
    for i,directory in enumerate(output_dirs(folder)):
        data,meta=read_float(directory/'data'/filename,skip_initial=bool(i))
        take=data['time_s']>last+1e-9
        if take.any():
            chunks.append({k:v[take] for k,v in data.items()});last=float(data['time_s'][take][-1])
        metadata.append(meta)
    data={key:np.concatenate([q[key] for q in chunks]) for key in chunks[0]}
    assert np.all(np.diff(data['time_s'])>0)
    assert all(abs(q['mass_kg']-metadata[0]['mass_kg'])<1e-4 for q in metadata)
    meta={**metadata[0], 'sources':metadata,'samples':len(data['time_s']),
        'time_range_s':data['time_s'][[0,-1]].tolist(),
        'max_force_balance_residual_N':max(q['max_force_balance_residual_N'] for q in metadata),
        'restart_policy':'Retain original checkpoint values; discard duplicate resumed initial placeholder; join at actual physical timestamps.'}
    return data,meta


def acceleration_rows(folder,marker):
    result=[];last=-np.inf
    for directory in output_dirs(folder):
        with (directory/f'FloatingAce_mkbound_{marker}.csv').open() as stream:
            rows=[r for r in csv.DictReader(stream,delimiter=';') if r.get('predictor')=='False' and float(r['time [s]'])>max(0,last)+1e-9]
        result.extend(rows)
        if rows:last=float(rows[-1]['time [s]'])
    return result
