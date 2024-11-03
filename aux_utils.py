import numpy as np
import bempp
from bempp.api.linalg import lu

def planewave(f,d,tlag,x,T,M):
    '''
    Input:
    f: function handle of 1d smooth function
    d: direction of propagation
    tlag: a potential lag (shift to not start at time t=0)
    x: points in space, where you want to have the wave
    normal: the exterior unit normals to the boundary
    T: final time 
    M: number of steps in time
    c: coefficient from the RK method
    Output:
    u: the plane wave at the stages needed for the CQ at the points x
    '''
    t = np.linspace(0,T,M+1);
    d = d[:,np.newaxis]
    d = np.tile(d,(1,x.shape[1]))
    direction = np.sum(x*d,axis=0)
    dir_res = direction[:,np.newaxis]
    Ui = np.zeros([direction.size, t.size])
    Ui = f(dir_res - (t-tlag))
    return Ui

def sphericalwave(f,tlag,x,T,M,RK=0,c=0):
    '''
    Input:
    f: function handle of 1d smooth function
    tlag: a potential lag (shift to not start at time t=0)
    x: points in space, where you want to have the wave
    normal: the exterior unit normals to the boundary
    T: final time 
    M: number of steps in time
    c: coefficient from the RK method
    Output:
    u: the spherical wave at the stages needed for the CQ at the points x
    '''
    t = np.linspace(0,T,M+1);
    Ui = np.zeros([x.shape[1], t.size])
    normx = np.sqrt(np.sum(x**2,0))
    normx = normx[:,np.newaxis]
    Ui = f(t - normx)/(4*np.pi*normx)
    return Ui


def mollifier(x):
    '''
    Input:
    x: some x values
    Output:
    q: this will be the function, in this case a mollifier
    '''
    Index = x**2<0.999999
    xn = x*Index
    q =  np.exp(-1/(1-xn**2))*Index
    return q


def createIncWave(points, A = 2, d = np.array([1.0, 0.0, 0.0]), tlag = 2.0, T=5, M=100,option=0):
    '''
    Create a wave in the plane z=0, which is either an incoming plane wave or a spherical wave. 
    In the latter case, the wave will be the solution of the scattering problem and therefore, well-suited for testing.
    For option=0 the wave is given by
        f(x\cdotd - (t-tlag)) (incoming plane wave with direction of propagation d \in S^2)
    For option=1 the wave is given by
        \phi(t-|x|)/(4\pi|x|) (spherical wave in 3d, see p. 78 Banjai/Sayas Book)
    '''
    fn = lambda x: mollifier(A*x)
    if option == 0:
        vals = planewave(fn, d, tlag, points, T, M)
    elif option == 1:
        vals = sphericalwave(fn, 0, points, T, M)
    return vals

def createIncWaveOnScatterer(points, A = 2.0, d = np.array([1.0, 0.0, 0.0]), tlag = 2.0, T=5, M=100,option = 0):
    '''
    Create a wave on the surface of the scatterer, which is either an incoming plane wave or a spherical wave. 
    In the latter case, the wave will be the solution of the scattering problem and therefore, well-suited for testing.
    For option=0 the wave is given by
        f(x\cdotd - (t-tlag)) (incoming plane wave with direction of propagation d \in S^2)
    For option=1 the wave is given by
        \phi(t-|x|)/(4\pi|x|) (spherical wave in 3d, see p. 78 Banjai/Sayas Book)
    '''
    x, y, z = points
    fn = lambda x: mollifier(A*x)
    vals = np.zeros((points.shape[1],M+1))
    if option == 0:
        vals = planewave(fn, d, tlag, points, T, M)
    elif option == 1:
        vals = sphericalwave(fn, 0, points, T, M)
    return vals


def getpointsonslice(extend=3, n1=100, n2=100, slice_to_0='z'):
    xmin, xmax, ymin, ymax = [-extend, extend, -extend, extend]
    plot_grid = np.mgrid[xmin:xmax:n1 * 1j, ymin:ymax:n2 * 1j]
    if slice_to_0 == 'x':
        points = np.vstack((np.zeros(plot_grid[0].size),
                            plot_grid[0].ravel(),
                            plot_grid[1].ravel()))
    elif slice_to_0 == 'y':
        points = np.vstack((plot_grid[0].ravel(),
                            np.zeros(plot_grid[0].size),
                            plot_grid[1].ravel()))
    elif slice_to_0 =='z':      
        points = np.vstack((plot_grid[0].ravel(),
                            plot_grid[1].ravel(),
                            np.zeros(plot_grid[0].size)))
    elif slice_to_0 == 'yz':
        points1 = np.vstack((plot_grid[0].ravel(),
                                np.zeros(plot_grid[0].size),
                                plot_grid[1].ravel()))
        points2 = np.vstack((plot_grid[0].ravel(),
                            plot_grid[1].ravel(),
                            np.zeros(plot_grid[0].size)))
        points = np.concatenate((points1, points2), axis=1)
    return points


def ConvolutionQuadrature(g,tau,delta,points,space):
    '''
    Computes the convolution quadrature.
    Input:
        g - the input data (usually a incoming wave)
        tau - time step size T/M
        delta - the transition function to determine the wave numbers
        points - we want to have the scattered wave here later on
    Output:
        u_out - the scattered wave at all points determined by 'points' (in the first dimension) 
                and for all times t_n (in the second dimension) 
    '''






