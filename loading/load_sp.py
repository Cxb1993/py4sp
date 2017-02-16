import numpy as np
import fft_sp as fft
import os
import ops_sp as ops

class setup: 
    def __init__(self, path='./'):
        data = []
        with open(path+'/NS.setup','r') as f:
            for line in f.readlines():
                data.append(line.replace('\n','').split(' '))

        self.case   = data[1][0]
        self.thermo = data[3][0]
        self.ekman  = data[5][0]
        self.Nx2    = int(data[7][0])
        self.Ny     = int(data[7][1])
        self.Lx     = float(data[9][0].replace('d','e'))
        self.Ly     = float(data[9][1].replace('d','e'))
        self.tstart = float(data[15][0].replace('d','e'))
        self.tstop  = float(data[17][0].replace('d','e'))
        self.dt1    = float(data[19][0].replace('d','e'))
        self.dt2    = float(data[21][0].replace('d','e'))
        self.zmeshf = data[25][0]
        if(os.path.exists(path+self.zmeshf)):
            self.zmesh_cc  = np.loadtxt(path+self.zmeshf,skiprows=1)[::2]
            self.zmesh_st  = np.loadtxt(path+self.zmeshf,skiprows=2)[::2]
        else: 
            print('ZMESH file not found.')
            self.zmesh_cc = '0'
            self.zmesh_st = '0'

def load_zmesh(filename):
    dummy = np.loadtxt(filename)
    Nz = dummy[0]
    zmesh_st = dummy[1::2]
    zmesh_cc = dummy[2::2]

    return zmesh_cc, zmesh_st

def load_windfarm(filename):
    dummy = np.loadtxt(filename, skiprows=2)
    windfarm = [] # Empty list
    for turbine_data in dummy:
        windfarm.append(Turbine(turbine_data))
    return windfarm

def load_NSp(filename):
    dummy = np.loadtxt(filename, skiprows=1)
    NSp = {}
    NSp['t'] = dummy[:,0]
    NSp['U'] = dummy[:,1]
    NSp['V'] = dummy[:,2]
    NSp['W'] = dummy[:,3]
    NSp['Etot'] = dummy[:,4]
    NSp['Eturb'] = dummy[:,5]
    NSp['sf_mean_x'] = dummy[:,6] # based on integral of right-hand side stress terms
    NSp['sf_mean_y'] = dummy[:,7] # based on integral of right-hand side stress terms
    NSp['sf_Um'] = dummy[:,8] # based on log law an <U>[0]
#    NSp['sf_sgs'] = dummy[:,6]
    return NSp

def load_BLcc(filename):
    dummy = np.loadtxt(filename, skiprows=2)
    BLinst = {}
    BLinst['z']  = dummy[:,0]
    BLinst['U']  = dummy[:,1]
    BLinst['V']  = dummy[:,2]
    BLinst['W']  = dummy[:,3]
    BLinst['uu'] = dummy[:,4]
    BLinst['vv'] = dummy[:,5]
    BLinst['ww'] = dummy[:,6]
    BLinst['uv'] = dummy[:,7]
    BLinst['uw'] = dummy[:,8]
    BLinst['vw'] = dummy[:,9]
    return BLinst

def load_plane_ascii(filename, N1, N2):
    # Load the data into memory
    u = np.loadtxt(filename)
    # Cast into vector
    u = np.reshape(u, [np.size(u), 1])
    # Rearrange, taking into account Fortran memory layout
    u = np.reshape(u, [N1, N2], order='F')
    return u

def load_plane(filename, N1, N2):
    with open(filename,'rb') as binfile:
        u = np.fromfile(binfile, dtype=np.float64).reshape((N1, N2), order='F')
    return u

def load_spectral_field(filename, N1, N2, N3, N4):
    dummy = np.loadtxt(filename, skiprows=7)
    spectral = {}
    d = dummy[:,0] + dummy[:,1]*1j
    spectral['u'] = np.reshape(d[0:N1*N2*N3], [N1, N2, N3], order='F')
    spectral['v'] = np.reshape(d[N1*N2*N3:2*N1*N2*N3], [N1, N2, N3], order='F')
    spectral['w'] = np.reshape(d[2*N1*N2*N3:], [N1, N2, N3-1], order='F')
    return spectral

def load_field(filename, N1, N2, N3, N4):
    # Load the data into memory
    u = np.loadtxt(filename)
    # Cast into vector
    u = np.reshape(u, [np.size(u), 1])
    # Rearrange, taking into account Fortran memory layout
    u = np.reshape(u, [N1, N2, N3, N4], order='F')
    return u

def load_windpower(filename,Nrows,Ncols):
    dummy = np.loadtxt(filename,comments='%')
    time = dummy[:,0]
    powerdum = -dummy[:,2::2]
    # Reshape data into windfarm format, this assumes the turbines are numbered row per row
    power = np.zeros((Nrows, Ncols, time.size))
    for num in range(Nrows*Ncols):
        row = int(num)//int(Ncols)
        col = num - row*Ncols
        power[row, col, :] = powerdum[:, num]
    return time, power

def load_BLfieldstat(filename, N1, N2, N3):
    N4 = 11
    # Load the data into memory
    u = np.loadtxt(filename, skiprows=1)
    # Cast into vector
    u = np.reshape(u, [np.size(u), 1])
    # Create dictionary
    stat = {}
    u = np.reshape(u, [N1, N2, N3, N4], order='F')
    stat['u'] = u[:,:,:,0]
    stat['v'] = u[:,:,:,1]
    stat['w'] = u[:,:,:,2]
    stat['uu'] = u[:,:,:,3]
    stat['vv'] = u[:,:,:,4]
    stat['ww'] = u[:,:,:,5]
    stat['uv'] = u[:,:,:,6]
    stat['uw'] = u[:,:,:,7]
    stat['vw'] = u[:,:,:,8]
    stat['p'] = u[:,:,:,9]
    stat['wst'] = u[:,:,:,10]
    return stat

def load_1D_spectrum(filename, L, N, Nz):
    dummy = np.loadtxt(filename, comments='%')
    spec = {}
    spec['k']  = np.array([(i)/L*(2*np.pi) for i in range(N//2)])
    spec['z']  = dummy[0:Nz-1,0]
    spec['uu'] = dummy[0:Nz-1,1:]
    spec['vv'] = dummy[Nz:2*Nz-1,1:]
    spec['zst']= dummy[2*Nz-1:,0]
    spec['ww'] = dummy[2*Nz-1:,1:]
#    if(np.size(spec['uu'][1,:]) > np.size(spec['k'])):
#        spec['uu'] = np.delete(spec['uu'],np.size(spec['uu'][1,:])-1, axis=1)
#        spec['vv'] = np.delete(spec['vv'],np.size(spec['uu'][1,:])-1, axis=1)
#        spec['ww'] = np.delete(spec['ww'],np.size(spec['uu'][1,:])-1, axis=1)
    return spec

def load_field_bin(filename, N1, N2, N3):
    stat = {}
    with open(filename, 'rb') as binfile:
        dumm = np.fromfile(binfile, dtype=np.float64)
        shape = (N1,N2,N3,3)
        dumm = dumm.reshape(shape, order='F')
        stat['u']   = dumm[:,:,:,0]
        stat['v']   = dumm[:,:,:,1]
        stat['w']   = dumm[:,:,:,2]
    return stat

def load_3Dvelocity(filename, N1, N2, N3, N4=3, precision='double', verbose=True):
    shape = (N1, N2, N3, N4)
    if precision=='double':
        dtype = np.float64
    elif precision=='single':
        if verbose:
            print('Single precision')
        dtype = np.float32
    with open(filename, 'rb') as binfile:
        velocity = np.fromfile(binfile, dtype=dtype, count=N1*N2*N3*N4).reshape(shape, order='F')
    return velocity

def load_BLfieldstat_bin(filename, N1, N2, N3, N4=11, Nload=11):
    stat = {}
    with open(filename, 'rb') as binfile:
        stat['nsamp'] = np.fromfile(binfile, dtype=np.int32, count=1)
        stat['time_interv'] = np.fromfile(binfile, dtype=np.float32, count=1)
        stat['time_incurr'] = np.fromfile(binfile, dtype=np.float32, count=1)
        dumm = np.fromfile(binfile, dtype=np.float64, count=N1*N2*N3*Nload)
        shape = (N1,N2,N3,Nload)
        dumm = dumm.reshape(shape, order='F')
        names = ['u', 'v', 'w', 'uu', 'vv', 'ww', 'uv', 'uw', 'vw', 'p', 'wst']
        for i in range(Nload):
            stat[names[i]] = dumm[:,:,:,i]
    return stat

def load_BLfield_real(filename, setuppath='./', **kwargs):
    if('verbose' in kwargs):
        verbose = kwargs['verbose']
    else:
        verbose = True
    if('Nx' in kwargs):
        BL = load_BLfield(filename,Nx=kwargs['Nx'], Ny=kwargs['Ny'], Nz=kwargs['Nz'], verbose=verbose)
    else:
        BL = load_BLfield(filename, verbose=verbose)

            
    
    if verbose:
        print('Performing c2r ffts')
    if('k' in kwargs):
        k = kwargs['k']
        BL['u']  = fft.c2r(BL['uu'][:,:,k], BL['Nx2'], BL['Ny'])
        BL['v']  = fft.c2r(BL['vv'][:,:,k], BL['Nx2'], BL['Ny'])
        BL['w']  = fft.c2r(BL['ww'][:,:,k], BL['Nx2'], BL['Ny'])
        del BL['uu'], BL['vv'], BL['ww'], BL['kx'], BL['ky']
    else:
        BL['u']  = fft.c2r(BL['uu'], BL['Nx2'], BL['Ny'])
        BL['v']  = fft.c2r(BL['vv'], BL['Nx2'], BL['Ny'])
        BL['w']  = fft.c2r(BL['ww'], BL['Nx2'], BL['Ny'])
        del BL['uu'], BL['vv'], BL['ww'], BL['kx'], BL['ky']

    # Read the setup file
    if not setuppath==None:
        stp = setup(path=setuppath)
    
        if not (stp.Nx2==BL['Nx2'] and 
                stp.Ny==BL['Ny'] ):
                print('NS.setup and fieldfile dimensions inconsistent.')
        else:
            # Interpolate w to cc
            BL['wcc'] = ops.interpolate_st2cc(BL['w'], stp.zmesh_cc, stp.zmesh_st)

    return BL

def load_BLfield_real_ascii(filename, N1, N2, N3):
    N1h = N1/2+1
    BL = {}
    dum = np.loadtxt(filename, comments='%')
    dum = dum[:,0] + 1j*dum[:,1]

    amount = N1h*N2*N3
    shape  = (N1h, N2, N3)
    shape2 = (N1h, N2, N3-1)
    uu = dum[:amount].reshape(shape, order='F')
    vv = dum[amount:2*amount].reshape(shape, order='F')
    ww = dum[2*amount:].reshape(shape2, order='F')
    
    BL['u']  = fft.c2r(uu, N1, N2)
    BL['v']  = fft.c2r(vv, N1, N2)
    BL['w']  = fft.c2r(ww, N1, N2)

    return BL

def cube_show_slider(cube, axis=2, **args):
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
    s = [slice(0, 1) if i == axis else slice(None) for i in xrange(3)]
    im = cube[s].squeeze()

    # display image
    l = ax.imshow(im, **args)
    axcolor = 'lightgoldenrodyellow'
    ax = fig.add_axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
        
    slider = Slider(ax, 'Axis %i index' % axis, 0, cube.shape[axis] - 1,
                                        valinit=0, valfmt='%i')
            
    def update(val):
        ind = int(slider.val)
        s = [slice(ind, ind + 1) if i == axis else slice(None) for i in xrange(3)]
        im = cube[s].squeeze()
        l.set_data(im, **args)
        fig.canvas.draw()
        
    slider.on_changed(update)
                                                                
    plt.show()

def load_BLfield(filename, real=False, log=False, cut=False, verbose=True, **kwargs):
    BL = {}
    with open(filename, 'rb') as binfile:
        BL['time'] = np.fromfile(binfile, dtype=np.float64, count=1)
        BL['Lx'] = np.fromfile(binfile, dtype=np.float64, count=1)
        BL['Ly'] = np.fromfile(binfile, dtype=np.float64, count=1)
        BL['Nx2'] = np.fromfile(binfile, dtype=np.int32, count=1)
        BL['Ny'] = np.fromfile(binfile, dtype=np.int32, count=1)
        BL['Nz'] = np.fromfile(binfile, dtype=np.int32, count=1)
        BL['thetaground'] = np.fromfile(binfile, dtype=np.float64, count=1)
        dum = np.fromfile(binfile,dtype=np.complex128)
    if('Nx' in kwargs):
        BL['Nx2'] = kwargs['Nx']
        BL['Ny'] = kwargs['Ny']
        BL['Nz'] = kwargs['Nz']
    if('Lx' in kwargs):
        BL['Lx'] = kwargs['Lx']
    if('Ly' in kwargs):
        BL['Ly'] = kwargs['Ly']
    if verbose:
        print( '######################' )
        print( 'BL_field.dat data:' )
        print( 'time         = ', BL['time'] )
        print( 'Lx           = ', BL['Lx'] )
        print( 'Ly           = ', BL['Ly'] )
        print( 'Nx2          = ', BL['Nx2'] )
        print( 'Ny           = ', BL['Ny'] )
        print( 'Nz           = ', BL['Nz'] )
        print( 'thetaground  = ', BL['thetaground'] )
        print( 'restsize     = ', dum.size )
        print( '######################' )

    N1 = int(BL['Nx2']/2+1)
    N2 = int(BL['Ny'])
    N3 = int(BL['Nz'])

    amount = N1*N2*N3
    amount2 = N1*N2*(N3-1)
    shape  = (N1, N2, N3)
    shape2 = (N1, N2, N3-1)
    BL['uu'] = dum[:amount].reshape(shape, order='F')
    BL['vv'] = dum[amount:2*amount].reshape(shape, order='F')
    BL['ww'] = dum[2*amount:2*amount+amount2].reshape(shape2, order='F')
    BL['kx'] = np.array([(i)/BL['Lx']*(2*np.pi) for i in range(N1)])
    BL['ky'] = np.array([(i)/BL['Ly']*(2*np.pi) for i in range(int(-N2/2+1), int(N2/2+1))])
    
    BLpostkeys = ['uu','vv','ww']
    if cut:
        BL['kx'] = BL['kx'][:-1]
        BL['ky'] = BL['ky'][:-1]

        for key in BLpostkeys:
        # Shift spectrum to match correct locations
            BL[key] = np.fft.fftshift(BL[key],axes=(1,))
        # Remove the defunct wavenumbers
            BL[key] = BL[key][:N1-1, 1:] 
    if real:
        for key in BLpostkeys:
            BL[key] = np.real(BL[key])

    if log:
        for key in BLpostkeys:
        # Take the logs
            BL[key] = np.log(np.abs(BL[key]))   



    return BL

def head_BLfield(filename):
    head = {}
    with open(filename, 'rb') as binfile:
        head['time'] = np.fromfile(binfile, dtype=np.float64, count=1)
        head['Lx'] = np.fromfile(binfile, dtype=np.float64, count=1)
        head['Ly'] = np.fromfile(binfile, dtype=np.float64, count=1)
        head['Nx2'] = np.fromfile(binfile, dtype=np.int32, count=1)
        head['Ny'] = np.fromfile(binfile, dtype=np.int32, count=1)
        head['Nz'] = np.fromfile(binfile, dtype=np.int32, count=1)
        head['thetaground'] = np.fromfile(binfile, dtype=np.float64, count=1)
    return head
