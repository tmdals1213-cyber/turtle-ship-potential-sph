"""Read little-endian DualSPHysics 5.4 floating output without altering the solver.

Format reference: upstream JBinaryData.cpp/.h and JPartFloatInfoBi4.cpp at
ef3721a861fda961f0e2f9ec4cd317b19de99086 (DualSPHysics, LGPL-2.1-or-later).
This reader validates lengths/types; it does not scan for arbitrary byte markers.
"""
from pathlib import Path
import argparse, json, struct
import numpy as np

TYPES={2:('i4',1),3:('i1',1),4:('u1',1),5:('i2',1),6:('u2',1),
       7:('i4',1),8:('u4',1),9:('i8',1),10:('u8',1),11:('f4',1),12:('f8',1),
       20:('i4',3),21:('u4',3),22:('f4',3),23:('f8',3)}


class Reader:
    def __init__(self,path):
        self.data=Path(path).read_bytes(); self.pos=64
        assert self.data.startswith(b'#FileJBD JPartFloatInfoBi4')
        assert self.data[60] in (0,10), 'Only little-endian files supported'
        self.count_format='Q' if self.data[61] else 'I'

    def take(self,n):
        assert n>=0 and self.pos+n<=len(self.data), 'Incomplete binary record'
        value=self.data[self.pos:self.pos+n]; self.pos+=n; return value

    def number(self,fmt='I'):
        return struct.unpack('<'+fmt,self.take(struct.calcsize('<'+fmt)))[0]

    def text(self):
        return self.take(self.number()).decode('utf-8')

    def typed(self,typ,count):
        if typ==1: return [self.text() for _ in range(count)]
        dtype,width=TYPES[typ]; dtype=np.dtype('<'+dtype)
        a=np.frombuffer(self.take(dtype.itemsize*width*count),dtype=dtype).copy()
        return a.reshape(count,width) if width>1 else a

    def item(self):
        size=self.number(); end=self.pos+size
        assert self.text()=='\nITEM\n'
        name=self.text(); self.number(); self.number(); self.text(); self.text()
        na,ni,nv=self.number(),self.number(),self.number()
        assert self.pos==end
        values={}; end=self.pos+nv
        if nv:
            assert self.text()=='\nVALUES'
            for _ in range(self.number()):
                key=self.text(); typ=self.number(); value=self.typed(typ,1)[0]
                values[key]=value.item() if isinstance(value,np.generic) else value
        assert self.pos==end
        arrays={}
        for _ in range(na):
            size=self.number(); end=self.pos+size
            assert self.text()=='\nARRAY'
            key=self.text(); self.number(); typ=self.number()
            count=self.number(self.count_format); nbytes=self.number(self.count_format)
            assert self.pos==end
            start=self.pos; arrays[key]=self.typed(typ,count)
            assert self.pos-start==nbytes
        children=[self.item() for _ in range(ni)]
        return dict(name=name,values=values,arrays=arrays,children=children)


def read_float(path,skip_initial=False):
    r=Reader(path); records=[]
    while r.pos<len(r.data): records.append(r.item())
    head=records[0]; assert head['values']['FtCount']==1, 'This project expects one vessel'
    series={}
    for rec in records[1:]:
        a=rec['arrays']; v=rec['values']
        times=np.asarray(a.get('TimeStep',[v.get('TimeStep')]),float)
        if times.ndim!=1 or not np.isfinite(times).all(): raise ValueError('Invalid timestamps')
        n=len(times); series.setdefault('time_s',[]).append(times)
        for key in ['center','fvel','fomega','facelin','faceang','extforcelin','extforceang',
                    'fluforcelin','fluforceang','preacelin','preaceang']:
            assert a[key].shape==(n,3)
            assert np.isfinite(a[key]).all()
            series.setdefault(key,[]).append(a[key])
    data={key:np.concatenate(value) for key,value in series.items()}
    initial_time=float(data['time_s'][0])
    if skip_initial:data={key:value[1:] for key,value in data.items()}
    assert np.all(np.diff(data['time_s'])>0), 'Nonmonotonic native samples'
    mass=float(head['arrays']['Mass'][0]); use=data['time_s']>0
    balance=mass*data['preacelin'][use]-data['fluforcelin'][use]-data['extforcelin'][use]-[0,0,-mass*9.81]
    error=float(np.abs(balance).max())
    assert error<.01, f'Native force balance mismatch: {error} N'
    meta=dict(source=str(path),samples=len(data['time_s']),mass_kg=mass,
        saved_interval_s=head['values']['TimeOut'],time_range_s=data['time_s'][[0,-1]].tolist(),
        max_force_balance_residual_N=error,frame='world',moments_reference='instantaneous CG',
        initial_placeholder_skipped=skip_initial,initial_file_time_s=initial_time,
        note='Native fluid moments are read directly, not inferred using a scalar inertia.')
    return data,meta


if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('case'); ap.add_argument('--extra',action='store_true')
    args=ap.parse_args(); folder=Path(__file__).resolve().parent/'sph'/args.case
    path=folder/'output/data'/('PartFloatInfo2.ibi4' if args.extra else 'PartFloatInfo.ibi4')
    data,meta=read_float(path); suffix='_extra' if args.extra else ''
    np.savez_compressed(folder/f'native_float{suffix}.npz',**data)
    (folder/f'native_float{suffix}.json').write_text(json.dumps(meta,indent=2))
    print(json.dumps(meta,indent=2))
