"""Common occupancy-surface reconstruction used by accepted films."""
import numpy as np
import pyvista as pv
from scipy.ndimage import gaussian_filter

def surface(points):
    bounds=[(-3.2,3.2),(-1.2,1.2),(-1.2,1.2)]
    counts,edges=np.histogramdd(points,bins=[128,48,48],range=bounds)
    field=gaussian_filter(counts.astype(np.float32),.75)
    grid=pv.ImageData(dimensions=field.shape,spacing=(.05,.05,.05),origin=(-3.175,-1.175,-1.175))
    grid['occupancy']=field.ravel(order='F')
    # A fixed kernel and threshold are used for every case and every time.
    mesh=grid.contour([.55],scalars='occupancy').smooth(n_iter=6,relaxation_factor=.1)
    return mesh
