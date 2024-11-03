import bempp.api
import numpy as np
from scipy.io import savemat
import aux_utils

# construct the plane, in which we want to see the waves
Nx, Ny = 100, 100
extend = 3
slice_to_0 = 'yz'
points = aux_utils.getpointsonslice(extend=extend, n1=Nx, n2=Ny, slice_to_0=slice_to_0)
refinement_of_grid = 2

# define all the parameters for the incoming wave
A = 2.0 # makes the support of the mollifier smaller, but steeper (multiplication of the argument)
d = np.array([1.0, 0.0, 0.0]) # direction of propagation of the plane wave
tlag = 2.0 # the wave lags behind a little bit
T = 5.0 # the final time
M = 100 # number of time steps
tau = T/M
option = 0 # this is for the plane wave
# all things concerning the grid
grid = bempp.api.shapes.regular_sphere(refinement_of_grid) # load the grid. we consider here the unit sphere. argument is the refinement level
dof_pos = grid.vertices # these are the vertices, thats were the dofs are for P
bempp.api.export('sphere.msh', grid=grid) #export the grid for Matlab magic
continuous_space = bempp.api.function_space(grid, "P", 1) # if you change this, then also change dof_pos to grid.centroids.T
ui_on_scat = aux_utils.createIncWaveOnScatterer(dof_pos,T=T, M=M,option=option) # create the incoming wave on the grid

# now finally, define the incoming wave on the plane
x, y, z = points
idx = np.sqrt(x**2 + y**2 + z**2) > 1.0 #exclude all the points inside the scatterer
evaluation_points = points[:, idx]
ui_on_plane = np.zeros((points.shape[1],M+1))
ui_on_plane[:] = np.nan
ui_on_plane[idx,:] = aux_utils.createIncWave(evaluation_points, T=T, M=M,option=option)

# Define the multi step method used in the cq method
delta = lambda z : (1-z) + 1/2*(1-z)**2
#
us_on_plane = np.zeros((points.shape[1],M+1))
us_on_plane[:] = np.nan
us_on_plane[idx,:] = aux_utils.ConvolutionQuadrature(g=-ui_on_scat,
                                   tau=tau,
                                   delta=delta,
                                   points=evaluation_points,
                                   space=continuous_space)

# now save it all
data = {'us_on_plane' : us_on_plane,
        'ui_on_plane' : ui_on_plane,
        'ui_on_scat' : ui_on_scat,
        'n1' : Nx,
        'n2' : Ny,
        'slice_to_0' : slice_to_0,
        'extend' : extend,
        'M' : M,
        'T' : T
        }
savemat('cqdata.mat',data)















