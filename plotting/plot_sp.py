import numpy as np
import load_sp as lsp
import matplotlib.pyplot as plt
import os
import windfarm as wf

def set_equal_tight(ax=plt.gca()):
    ax.set_aspect('equal')
    ax.autoscale(tight=True)

def plot_field_turbines(fieldfile='BL_field.dat', key='u', k=16):
    bl = lsp.load_BLfield_real(fieldfile)
    field = bl[key][:,:,k]
    farm = wf.Windfarm()

    plt.figure()
    plt.imshow(np.flipud(np.transpose(field)), extent=(0, bl['Lx'], 0, bl['Ly']))
    farm.plot_turbines()
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.xlim((0, bl['Lx']))
    plt.ylim((0, bl['Ly']))
    

def plot_bl(bl, kplot=14, iplot=0, jplot=0, key='u', Lz = 1):
    plt.figure()

    plt.subplot(311)
    plt.imshow(np.flipud(np.transpose(bl[key][:,:,kplot])), extent=(0, bl['Lx'], 0, bl['Ly']))
    plt.colorbar()

    plt.subplot(312)
    plt.imshow(np.flipud(np.transpose(bl[key][iplot,:,:])), extent=(0, bl['Ly'], 0, Lz))
    plt.colorbar()

    plt.subplot(313)
    plt.imshow(np.flipud(np.transpose(bl[key][:,jplot,:])), extent=(0, bl['Lx'], 0, Lz))
    plt.colorbar()
    

def plot_turbines_topview(filename):
    turbines = np.loadtxt(filename, skiprows=2)
    for turbine in turbines:
        xcoord = turbine[0]
        ycoord = turbine[1]
        radius = turbine[4]
        ycoords = np.array([ycoord - radius, ycoord + radius])
        plt.plot((xcoord, xcoord), ycoords, 'k', lw=2)

def movie_xy(k, dt, var='u', setuppath='./../', pausetime=0.1, windfarmyaw=True, **kwargs):
    """
    function movie_xy

    Parameters
    ------------
    k: int
        grid index in zmesh where slices should be plotted
    dt: float
        timestep betwee snapshots
    var: str, optional
        variable to be plotted, default is 'u'
    setuppath : str, optional
        path where setupfile is located, default is './../'
    clim: tuple, optional kwarg
        colorbar limits for movie, default is (0,25)
    cmap: str, optional
        colormap used in movie, default is 'jet'
    tstart: float, optional
        initial time for movie snapshots
    tstop: float, optional
        final time for movie snapshots
    """
    setup = lsp.setup(setuppath)
    if 'clim' in kwargs:
        cl = kwargs['clim']
    else:
        cl = (0,25)
    if 'cm' in kwargs:
        cmap = kwargs['cm']
    else:
        cmap = 'jet'

    if windfarmyaw:
        farm = wf.Windfarm(path=setuppath)

    if 'tstop' in kwargs:
        t = np.arange(kwargs['tstart'], kwargs['tstop'], dt)
        print('Making movie for t = ', t)
        for tind, tim in enumerate(t):
            plt.clf()
            print('Plotting t =', tim)
            filename = var+'_zplane_k{:03d}_t_{:4.4f}.dat'.format(k, tim)
            plt.title(tim)
            plot_planexy(filename,show=False,prin=False,clim=cl,cm=cmap)
            if windfarmyaw:
                farm.plot_turbines_yaw(index=tind)
            plt.pause(pausetime)
    else:
        print('Automatic timeloop not yet implemented')


def plot_planexy(filename, Nx=0, Ny=0, show=True, prin=True,**kwargs):
    # First read from setup
    if(os.path.exists('./../NS.setup')):
        if prin:
            print('Reading grid dimensions from setup file')
        setup = lsp.setup('./../')
        Nxd = setup.Nx2
        Nyd = setup.Ny
        if prin:
            print('Nx = ', Nxd)
            print('Ny = ', Nyd)
    else:
        print('Taking grid dimensions from input parameters')
        Nxd = Nx
        Nyd = Ny
#    data = lsp.load_plane(filename, Nxd, Nyd)
    data = lsp.load_plane_single(filename, Nxd, Nyd)
    if 'cm' in kwargs:
        cmap = kwargs['cm']
    else:
        cmap = 'jet'
    plt.imshow(np.flipud(np.transpose(data)), extent=(0, setup.Lx, 0, setup.Ly),cmap=cmap, interpolation='bilinear'); plt.colorbar()
    if 'clim' in kwargs:
        plt.clim(kwargs['clim'])
    if show:
        plt.show()

def plot_planeyz(filename, Ny, Nz):
    return 0

def plot_planexz(filename, Nx, Nz):
    return 0

def make_movie(time_array,N1,N2):
    
    for t in time_array:
        plt.clf()
        tstr = "{:6.4f}".format(t)
        print( 'Loading t = ', tstr)
        basefilename = '_zplane_k013_t_'
        filenameu = 'u'+basefilename+tstr+'.dat'
        filenamev = 'v'+basefilename+tstr+'.dat'
        u = lsp.load_plane(filenameu, N1=N1, N2=N2)
        v = lsp.load_plane(filenamev, N1=N1, N2=N2)
        plt.pcolormesh(np.transpose(np.sqrt(u**2+v**2)))
        plt.clim((0, 30))
        plt.colorbar()
        plt.savefig('u_'+tstr+'.png')

