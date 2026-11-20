"""Boundary-particle load colors, baffle resultants, and vehicle roll loading."""
from pathlib import Path
import argparse,json,subprocess
import numpy as np
import pyvista as pv
from PIL import Image,ImageDraw
from matplotlib import colormaps
from render_coupled_results import CoupledScene,prepare,LABELS,CG0
from render_minimal_intro import W,H,FPS,BG,WHITE,BLUE,ORANGE,text

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'media'
CLIM=5000.


def arrow_mesh(start,vector,scale):
    mag=float(np.linalg.norm(vector));length=min(2.3,mag/scale)
    if length<.025:length=.025
    return pv.Arrow(start=start,direction=vector/max(mag,1e-12) if mag else (0,1,0),scale=length,
                    tip_length=.23,tip_radius=.1,shaft_radius=.035)


def moment_mesh(vector):
    magnitude=float(np.linalg.norm(vector))
    axis=vector/max(magnitude,1e-12) if magnitude else np.array([1.,0,0])
    seed=np.array([0.,1,0]) if abs(axis[1])<.9 else np.array([0.,0,1])
    u=seed-axis*np.dot(seed,axis);u/=np.linalg.norm(u);v=np.cross(axis,u)
    angle=np.linspace(0,max(.05,min(1.4*np.pi,magnitude/80000*1.4*np.pi)),50)
    points=np.array([-.3,0,1.45])+.78*(np.cos(angle)[:,None]*u+np.sin(angle)[:,None]*v)
    line=pv.lines_from_points(points).tube(radius=.035)
    tangent=-np.sin(angle[-1])*u+np.cos(angle[-1])*v
    cone=pv.Cone(center=points[-1]+.07*tangent,direction=tangent,height=.24,radius=.11,resolution=12)
    return line.merge(cone)


def exposed_baffle_faces(xyz,groups,dp=.04):
    """Exposed voxel faces carry their actual boundary particle's load (N).

    This is a boundary-particle load map, not a pressure traction per face.
    The union of all boundary cells removes hidden plate/wall intersections.
    """
    lattice=np.rint((xyz-xyz.min(0))/dp).astype(np.int64)
    stride=lattice.max(0)+4;offset=np.array([1,stride[0],stride[0]*stride[1]])
    keys=lattice@offset;selected=np.flatnonzero(groups>0)
    faces=[];owners=[]
    for axis in range(3):
        tangents=[a for a in range(3) if a!=axis]
        for sign in (-1,1):
            ids=selected[~np.isin(keys[selected]+sign*offset[axis],keys)]
            points=np.repeat(xyz[ids,None,:],4,axis=1)
            points[:,:,axis]+=sign*dp/2
            points[:,:,tangents[0]]+=np.array([-1,1,1,-1])*dp/2
            points[:,:,tangents[1]]+=np.array([-1,-1,1,1])*dp/2
            faces.append(points);owners.append(ids)
    vertices=np.concatenate(faces).reshape(-1,3);owners=np.concatenate(owners)
    cells=np.column_stack((np.full(len(owners),4),np.arange(len(vertices)).reshape(-1,4)))
    return pv.PolyData(vertices,cells.ravel()),owners


class LoadScene(CoupledScene):
    def __init__(self,folder,size):
        self.folder=folder;self.setup=json.loads((folder/'setup.json').read_text())
        layout=self.setup['layout'];super().__init__(size,layout,self.setup.get('numerical_plate_thickness_m',.08))
        self.loads=np.load(folder/'surface_loads/loads.npz')
        report=json.loads((folder/'surface_loads/analysis.json').read_text())
        assert report['qualified_for_load_visualization']
        self.water.GetProperty().SetOpacity(.07)
        # The original plate actors are only a faint geometry guide beneath the
        # real boundary-particle forces. For A, the colored boundary is the wall.
        nplates=0 if layout=='A_unpartitioned' else 3 if layout=='D_cross_baffles' else 2
        for actor in self.actors[6:6+nplates]:actor.GetProperty().SetOpacity(.07)
        groups=self.loads['group'];xyz=self.loads['initial_points'];half=self.setup.get('numerical_plate_thickness_m',.08)/2
        if not nplates:
            rad=np.linalg.norm(xyz[:,1:],axis=1)
            self.mask=((rad<1.041)&(abs(xyz[:,0])<3.041))&((rad>=1.)|(abs(xyz[:,0])>=3.))
        else:
            self.mask=((groups==1)&(abs(xyz[:,0]+1)>half-.04))|((groups==2)&(abs(xyz[:,0]-1)>half-.04))|((groups==3)&(abs(xyz[:,1])>half-.04))
        if nplates:
            self.cloud,self.owners=exposed_baffle_faces(xyz,groups)
            self.cloud.cell_data['load']=np.zeros(self.cloud.n_cells)
        else:
            # Remove the camera-near half to reveal the actual inner-wall field.
            self.mask&=xyz[:,1]>=0
            self.owners=np.flatnonzero(self.mask)
            self.cloud=pv.PolyData(xyz[self.owners]);self.cloud['load']=np.zeros(self.cloud.n_points)
        self.nplates=nplates
        self.load_actor=self.p.add_mesh(self.cloud,scalars='load',clim=(0,np.log10(1+CLIM)),cmap='jet',show_scalar_bar=False,
                                        style='surface' if nplates else 'points',point_size=4,render_points_as_spheres=False,
                                        ambient=1,diffuse=0,specular=0,lighting=False)
        self.actors.append(self.load_actor)
        self.side_arrow=self.add(arrow_mesh((.2,0,1.6),np.array([0,1,0]),100000),BLUE,lighting=False,ambient=1,diffuse=0)
        self.roll_arrow=self.add(moment_mesh(np.zeros(3)),ORANGE,lighting=False,ambient=1,diffuse=0)
        self.group_arrows=[]
        for j,centre in [(1,(-1,0,.25)),(2,(1,0,.25)),(3,(0,0,.5))]:
            if (groups==j).any():self.group_arrows.append((j,centre,self.add(arrow_mesh(centre,np.array([1,0,0]),70000),'#d4e6eb',opacity=.75,lighting=False)))

    def load_frame(self,index,pose,water,display_time=None):
        mesh=pv.read(self.folder/'surface_loads'/f'BoundaryLoad_{index:04d}.vtk');order=np.argsort(mesh['Id'])
        assert np.array_equal(mesh['Id'][order],self.loads['idp'])
        magnitude=np.linalg.norm(mesh['Force'][order],axis=1)
        values=np.log10(1+magnitude[self.owners])
        if self.nplates:self.cloud.cell_data['load']=values
        else:self.cloud.point_data['load']=values
        t=pose['time_s'][index] if display_time is None else display_time
        def interpolated(key):
            return np.array([[np.interp(t,self.loads['time_s'],self.loads[key][:,j,k]) for k in range(3)] for j in range(4)])
        force=interpolated('force_body_filtered_N');moment=interpolated('moment_body_filtered_Nm')
        total_f=force.sum(0);total_m=moment.sum(0)
        self.side_arrow.mapper.SetInputData(arrow_mesh((.2,0,1.6),total_f,100000))
        self.roll_arrow.mapper.SetInputData(moment_mesh(total_m))
        for j,centre,actor in self.group_arrows:actor.mapper.SetInputData(arrow_mesh(centre,force[j],70000))
        # Local load segment isolates loading, so no forward-road scroll is added.
        im=self.actual_frame(pose['rotation'][index],pose['center'][index],water,braking=pose['time_s'][index]>=1.5)
        return im,np.linalg.norm(total_f)/1000,np.linalg.norm(total_m)/1000


def render(names,preview=False):
    folders=[ROOT/'sph'/n for n in names]
    records=[prepare(f) for f in folders];scenes=[LoadScene(f,(510,480)) for f in folders]
    def frame(i,display_time=None):
        im=Image.new('RGB',(W,H),BG);d=ImageDraw.Draw(im)
        text(d,(65,150),'LIQUID LOAD',64,bold=True)
        text(d,(65,235),'ON THE VEHICLE',58,bold=True)
        text(d,(65,345),'TOTAL FORCE',28,fill=BLUE);text(d,(530,345),'TOTAL MOMENT',28,fill=ORANGE)
        text(d,(65,410),f'{records[0][0]["time_s"][i] if display_time is None else display_time:.1f} s',32,fill='#a3bbcc')
        for k,(scene,(p,cache),folder) in enumerate(zip(scenes,records,folders)):
            shot,force,moment=scene.load_frame(i,p,pv.read(cache/f'water_{i:04d}.vtp'),display_time)
            x=30+510*(k%2);y=535+570*(k//2)
            im.paste(shot,(x,y));d=ImageDraw.Draw(im)
            text(d,(x+20,y+10),LABELS[scene.setup['layout']],25,bold=True)
            text(d,(x+20,y+470),f'|F|  {force:.0f} kN',29,fill=BLUE)
            text(d,(x+20,y+510),f'|M|  {moment:.0f} kNm',29,fill=ORANGE)
        cmap=colormaps['jet'];bar=np.uint8(cmap(np.linspace(0,1,820))[None,:,:3]*255)
        im.paste(Image.fromarray(bar).resize((820,18)),(130,1695));d=ImageDraw.Draw(im)
        for value in [0,10,100,1000,CLIM]:
            x=130+820*np.log10(1+value)/np.log10(1+CLIM)
            text(d,(x,1720),f'{value:.0f}',22,fill='#9aafbe',anchor='mt')
        text(d,(130,1765),'LOCAL LOAD [N/PARTICLE]  /  LOG SCALE',23,fill='#9aafbe')
        text(d,(130,1810),'TOTAL FLUID LOADS  /  MOMENT ABOUT VEHICLE CG',20,fill='#748ea3')
        return im
    if preview:
        for i in [40,60,80]:frame(i).save(OUT/f'tanker_surface_loads_t{i/20:.1f}.png')
        frame(60).save(OUT/'tanker_surface_loads_preview.png')
    else:
        path=OUT/'tanker_surface_loads_comparison_portrait.mp4'
        proc=subprocess.Popen(['ffmpeg','-v','error','-y','-f','rawvideo','-pix_fmt','rgb24','-s',f'{W}x{H}',
            '-r','30','-i','-','-an','-c:v','libx264','-preset','fast','-crf','19','-pix_fmt','yuv420p','-movflags','+faststart',str(path)],stdin=subprocess.PIPE)
        # Saved local field / pose at20Hz; filtered resultants interpolate at30Hz.
        for k in range(241):
            t=k/30;i=int(np.argmin(abs(records[0][0]['time_s']-t)))
            proc.stdin.write(frame(i,t).tobytes())
            if k%60==0:print(f'loads {k}/240',flush=True)
        proc.stdin.close()
        if proc.wait():raise RuntimeError('Encoding failed')
        (OUT/'tanker_surface_loads_comparison.json').write_text(json.dumps(dict(video=path.name,cases=names,duration_s=241/30,
          local_field='Actual per-boundary-particle fluid force magnitude from official ComputeForces; N/particle, not pressure.',
          local_colors='Jet, common log10(1+force) scale0-5000N; includes zero and full observed range. Same DP40mm across cases.',
          local_selection='Exposed baffle voxel faces show their boundary particle load; A inner-wall cutaway. Integrate every boundary layer; no face-pressure claim.',
          arrows='White: each baffle net force. Cyan: total fluid-force vector. Orange: right-hand arc about total fluid-moment axis through dry CG; drawn offset above tank for visibility. Numbers are vector magnitudes. Common force scale100kN/m capped2.3m, moment arc80kNm/1.4pi capped1.4pi.',
          temporal='Saved local fields/poses20Hz; filtered integrated load vectors linearly interpolated30Hz.',
          filter='0.4s Savitzky-Golay integrated loads; colors use saved instantaneous local load. Not validated local slam peaks.'),indent=2))
        print(str(path),flush=True)
    for scene in scenes:scene.close()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('cases',nargs='+');p.add_argument('--preview',action='store_true')
    a=p.parse_args();render(a.cases,a.preview)
