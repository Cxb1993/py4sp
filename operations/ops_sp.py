import numpy as np
from scipy.interpolate import interp1d

def interpolate_st2cc(field, z_cc, z_st):
    field_interp = np.zeros((field.shape[0], field.shape[1], z_st.size), dtype=np.float64)

    # !! w is stored on z_cc locations, i.e. in between z_st..

    # add a slab on top and bottom of zeros..
    field_dum = np.zeros((field.shape[0], field.shape[1], z_cc.size))
    field_dum[:,:,0]    = np.zeros((field.shape[0], field.shape[1]))
    field_dum[:,:,-1]   = np.zeros((field.shape[0], field.shape[1]))
    field_dum[:,:,1:-1] = field

    # get the interpolating function
    f = interp1d(z_cc, field_dum, axis=2)

    # interpolate
    field_interp = f(z_st)

    return field_interp

def interpolate_cc2st(field, z_st, z_cc):
    field_interp = np.zeros((field.shape[0], field.shape[1], z_st.size), dtype=np.float64)

    print('Returning zeros in interpolate_cc2st')
    return field_interp

def stack_velocity(bl):
    """
    function stack_velocity

    Stacks velocity fields from bl dictionary to a 4D numpy array.

    """
    return np.stack((bl['u'], bl['v'], bl['wcc']), axis=3)

def calc_vorticity(bl, dx=1, dy=1, dz=1):
    """
    function calc_vorticity

    Calculates vorticity field using 2nd order FD approximations for gradients

    """
    # Vorticity = [ dw/dy - dv/dz, du/dz - dw/dx, dv/dx - du/dy ]

    # First, calculate the gradient tensor
    if dx==1: # no kwargs given.. automatically determine grid spacings
        dx = bl['Lx']/bl['Nx2']; dy = bl['Ly']/bl['Ny']; dz = 1/bl['Nz']
    velocity = stack_velocity(bl)
    grad = calc_gradient(velocity, dx, dy, dz)

    # Next, compose the vorticity field
    vorticity = np.zeros(velocity.shape, dtype=np.float64)
    vorticity[:,:,:,0] = grad[:,:,:,2,1] - grad[:,:,:,1,2] # dw/dy - dv/dz
    vorticity[:,:,:,1] = grad[:,:,:,0,2] - grad[:,:,:,2,0] # du/dz - dw/dx
    vorticity[:,:,:,2] = grad[:,:,:,1,0] - grad[:,:,:,0,1] # dv/dx - du/dy

    return vorticity

def calc_lambda2(bl, dx, dy, dz):
    """
    function calc_lambda2

    Calculates second-largest eigenvalue of S^2 + Omg^2, with S: rate of strain tensor, Omega: rate of rotation tensor
    Used as a vortex identification criterion.

    """
    S, Omg = calc_S_Omg(stack_velocity(bl), dx, dy, dz)
    A = S**2 + Omg**2

    lambda2 = np.zeros(bl['u'].shape)

    # This might be nasty... MANY calls to eig...
    print('Calculating lambda2')
    for i in range(lambda2.shape[0]):
        for j in range(lambda2.shape[1]):
            for k in range(lambda2.shape[2]):
                lambda2[i,j,k] = np.linalg.eigvals(A[i,j,k])[1]

    return lambda2

def calc_Q(bl, dx, dy, dz):
    """
    function calc_Q

    Calculates second invariant of S^2 + Omg^2, with S: rate of strain tensor, Omega: rate of rotation tensor
    Used as a vortex identification criterion.

    Q = 1/2 * ( | Omega |**2  - | S |**2 )

    """

    S, Omg = calc_S_Omg(stack_velocity(bl), dx, dy, dz)

    print(np.mean(S))

    print('----------------')

    print(np.mean(Omg))

    Q = np.zeros(bl['u'].shape)

    # Not sure whether there's a NumPy function that can calculate norm of a tensor... so manually..
    Omg_sq = np.zeros(bl['u'].shape)
    S_sq = np.zeros(bl['u'].shape)

    print('Calculating Q invariant')
    for i in range(3):
        for j in range(3):
            Omg_sq = Omg_sq + Omg[:,:,:,i,j]**2
            S_sq   = S_sq + S[:,:,:,i,j]**2

    Q = 0.5*(Omg_sq - S_sq)

    return Q

def calc_S_Omg(field, dx, dy, dz):
    """
    function calc_S_Omg

    Calculates (a)symmetric parts of grad(u)
    """
    gradu = calc_gradient(field, dx, dy, dz)
    gradut = np.zeros(gradu.shape)
    
    # perform a manual tranpose over final 2 dimensions
    for i in range(gradut.shape[3]):
        for j in range(gradut.shape[4]):
            gradut[:,:,:,i,j] = gradu[:,:,:,j,i]

    S = 1/2*(gradu + gradut)
    Omg = 1/2*(gradu - gradut)

    return S, Omg

def calc_gradient(field, dx, dy, dz):
    """
    function calc_gradient

    Calculates gradient tensor using 2nd order FD approximations

    Parameters
    ------------------
    field: 4D np.array, type = np.float64
        field of which gradient will be calculated
        
    dx: np.float
        grid spacing in x direction

    dy: np.float
        grid spacing in y direction

    dz: np.float
        grid spacing in z direction

    Returns
    ------------------
    gradient: 4D np.array, type = np.float64
        gradient tensor
    """
    gradient = np.zeros((field.shape[0], field.shape[1], field.shape[2], field.shape[3], field.shape[3]))
    for component_i in range(field.shape[3]): # This will virtually always be three, or two for 2D spaces
        g = np.gradient(field[:,:,:,component_i])
        for component_j in range(field.shape[3]):
            gradient[:, :, :, component_i, component_j] = g[component_j]
    return gradient



def dspectral_1D(field, dx, axis=0):
    """
    function dspectral_1D WM - UNTESTED, TEST THIS ON A SIMPLE SIGNAL e.g. SINE WAVE FIRST

    Calculates 1D derivative of periodic fields using spectral numerics

    Parameters
    --------------
    field: np.array
        real space array for which derivative will be calculated

    axis: int, optional
        axis over which derivative will be calculated (default = first axis)

    Returns
    --------------
    dfield: np.array
        real space array, derivative of input field
    """
    field_hat = np.fft.rfft(field, axis=axis)
    k         = np.fft.rfftfreq(np.shape(field)[axis], dx)*(2*np.pi)

    field_hatprime = np.zeros(field_hat.shape, dtype=np.complex128)
    if axis==0:
        for j, kj in enumerate(k):
            field_hatprime[j,:,:] = 1j*kj*field_hat[j,:,:]
    elif axis==1:
        for j, kj in enumerate(k):
            field_hatprime[:,j,:] = 1j*kj*field_hat[:,j,:]
    else:
        print('Thou shalt not calculate spectral derivatives in the vertical direction.')
        return 0

    dfield = np.fft.irfft(field_hatprime)

