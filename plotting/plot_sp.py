import matplotlib.pyplot as plt
import numpy as np
import load_sp as lsp
from load_sp import load_plane
import os
from matplotlib import cm 
from matplotlib.widgets import Slider

def plot_turbines_topview(filename):
    turbines = np.loadtxt(filename, skiprows=2)
    for turbine in turbines:
        xcoord = turbine[0]
        ycoord = turbine[1]
        radius = turbine[4]
        ycoords = np.array([ycoord - radius, ycoord + radius])
        plt.plot((xcoord, xcoord), ycoords, 'k', lw=2)

def movie_xy(k, dt, var='u', setuppath='./../', **kwargs):
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
    basename = var+'_zplane_k'+"{0:03d}".format(k)+'_t_'
    setup = lsp.setup(setuppath)
    plt.figure()
    plt.ion()
    if 'clim' in kwargs:
        cl = kwargs['clim']
    else:
        cl = (0,25)
    if 'cm' in kwargs:
        cmap = kwargs['cm']
    else:
        cmap = 'jet'

    if 'tstop' in kwargs:
        t = np.arange(kwargs['tstart'], kwargs['tstop'], dt)
        print('Making movie for t = ', t)
        for tim in t:
            plt.clf()
            print('Plotting t =', tim)
            filename = basename+"{:4.4f}".format(tim)+".dat"
            plt.title(tim)
            plot_planexy(filename,show=False,prin=False,clim=cl,cm=cmap)
            plt.pause(0.1)
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
    data = lsp.load_plane(filename, Nxd, Nyd)
    if 'cm' in kwargs:
        cmap = kwargs['cm']
    else:
        cmap = 'jet'
    plt.imshow(np.flipud(np.transpose(data)), extent=(0, setup.Lx, 0, setup.Ly),cmap=cmap); plt.colorbar()
    if 'clim' in kwargs:
        plt.clim(kwargs['clim'])
    if show:
        plt.show()

def plot_planeyz(filename, Ny, Nz):
    return 0

def plot_planexz(filename, Nx, Nz):
    return 0

def plot3dgui(data,**kwargs):

    fig = plt.figure()
    ax = plt.gca()

    p = ax.imshow(data[:,:,0], interpolation='none')
    
    # Add the slider
    plt.subplots_adjust(bottom=.25)
    slider_ax = plt.axes([.25, .15, .5, .03])
    
    slider = Slider(slider_ax, 'Index', 0, data.shape[2], valinit=0, valfmt='%0.0f')

    def update(val):
        index = val
        p.set_data(data[:,:,index])
        fig.canvas.draw()

    slider.on_changed(update)
    plt.show()
    return p, slider


def cube_show_slider(cube, axis=2, **kwargs):
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider, Button, RadioButtons
    # check dim
    if not cube.ndim == 3:
        raise ValueError("cube should be an ndarray with ndim == 3")
    # generate figure
    fig = plt.figure()
    ax = plt.subplot(111)
    fig.subplots_adjust(left=0.25, bottom=0.25)

    # select first image
    s = [slice(0, 1) if i == axis else slice(None) for i in range(3)]
    im = cube[s].squeeze()

    # display image
    l = ax.imshow(im, **kwargs)
    axcolor = 'lightgoldenrodyellow'
    ax = fig.add_axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
        
    slider = Slider(ax, 'Axis %i index' % axis, 0, cube.shape[axis] - 1,
                                        valinit=5, valfmt='%i')
            
    def update(val):
        ind = int(slider.val)
        s = [slice(ind, ind + 1) if i == axis else slice(None) for i in range(3)]
        im = cube[s].squeeze()
        l.set_data(im, **kwargs)
        fig.canvas.draw()
        
    slider.on_changed(update)
                                                                
    plt.show()

def plot_controls(turbine):
    basefile = 'f_params'
    count = 0
    cont = True
    while cont:
        filename = 'f_params'+str(count)+'.txt'
        if os.path.exists(filename):
            f = np.loadtxt(filename,skiprows=5)[:,1:] 
            plt.plot(f[turbine], label=count)
        else:
            cont = False
        count += 1
    plt.legend(loc=0)
    plt.show()


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

