import numpy as np
import fft_sp as fft
import os

def create_vapor_file(filename, bl, varstring='u'):

    # bl is a container dictionary, can originate from read_BLfield, or handmade, the following should be contained:
    Nx = bl['Nx2'][0]; Ny = bl['Ny'][0]; Nz = bl['Nz'][0]
    Lx = bl['Lx'][0];  Ly = bl['Ly'][0]; Lz = 1.0 # hardcoded height!
    field = bl[varstring]

    # First create the vdf container file
    print(Nx, Ny, Nz, Lx, Ly, Lz)
    vdf_createstring = 'vdfcreate -dimension {:d}x{:d}x{:d} -extents 0:0:0:{:4.4f}:{:4.4f}:{:4.4f} -vars3d '.format(Nx, Ny, Nz, Lx, Ly, Lz)+varstring+' '+filename

    print('VDF Creation: ')
    print(vdf_createstring)

    # Execute
    os.system(vdf_createstring)

    # Next, dump the field contents to file
    tmp_file = 'temp.tmp'
    with open(tmp_file, 'wb') as binfile_write:
        field.astype(np.float32).tofile(binfile_write)

    # Finally, add the fields to the vapor vdf container
    vdf_addstring = 'raw2vdf -varname '+varstring+' '+filename+' '+tmp_file
    print('VDF Add fields: ')
    print(vdf_addstring)
    # Execute
    os.system(vdf_addstring)
    os.system('rm '+tmp_file)

def create_constant_field(filename, Nx, Ny, Nz, Lx, Ly, Lz, value):
    """
    function create_constant_field

    Parameters
    -------------
    filename: str
        filename to which field will be written
    Nx, Ny, Nz: int
        grid size
    Lx, Ly, Lz:
        field size
    value: 
        constant value throughout field
    """
    print('##########################################')
    print('# Writing field to file ', filename)
    print('# ------------------------------------------')
    print('# Field information')
    print('# ------------------------------------------')
    print('# Lx    = ', Lx)
    print('# Ly    = ', Ly)
    print('# Nx2   = ', Nx)
    print('# Ny    = ', Ny)
    print('# Nz    = ', Nz)
    print('##########################################')

    with open(filename, 'wb') as binfile:
        np.array([0]).astype(dtype=np.float64).tofile(binfile)
        np.array([Lx]).astype(dtype=np.float64).tofile(binfile)
        np.array([Ly]).astype(dtype=np.float64).tofile(binfile)
        np.array([Nx]).astype(dtype=np.int32).tofile(binfile)
        np.array([Ny]).astype(dtype=np.int32).tofile(binfile)
        np.array([Nz]).astype(dtype=np.int32).tofile(binfile)
        np.array([0]).astype(dtype=np.float64).tofile(binfile)

        shapeu = (Nx, Ny, Nz)
        shapew = (Nx, Ny, Nz-1)
        u = fft.r2c(np.ones(shapeu, dtype=np.float64)*value, Nx, Ny).astype(np.complex128)
        v = fft.r2c(np.zeros(shapeu, dtype=np.float64), Nx, Ny).astype(np.complex128)
        w = fft.r2c(np.zeros(shapew, dtype=np.float64), Nx, Ny).astype(np.complex128)

        u.reshape((u.size,1), order='F').tofile(binfile)
        v.reshape((v.size,1), order='F').tofile(binfile)
        w.reshape((w.size,1), order='F').tofile(binfile)

def write_BLfield(field, filename, spectral=False):
    """
    function write_BLfield

    Parameters
    ------------
    field: dict
        structure produced by load_BLfield function from load_sp.py
    filename: str
        filename to which the current field data should be written
    spectral: bool, optional
        flag whether field data is in Fourier space (default is False)
    """
    print('##########################################')
    print('# Writing field to file ', filename)
    print('# Input field is spectral? ', spectral)
    print('# ------------------------------------------')
    print('# Field information')
    print('# ------------------------------------------')
    print('# time  = ', field['time'])
    print('# Lx    = ', field['Lx'])
    print('# Ly    = ', field['Ly'])
    print('# Nx2   = ', field['Nx2'])
    print('# Ny    = ', field['Ny'])
    print('# Nz    = ', field['Nz'])
    print('# theta = ', field['thetaground'])
    print('##########################################')

    with open(filename,'wb') as binfile:
        field['time'].astype(dtype=np.float64).tofile(binfile)
        field['Lx'].astype(dtype=np.float64).tofile(binfile)
        field['Ly'].astype(dtype=np.float64).tofile(binfile)
        field['Nx2'].astype(dtype=np.int32).tofile(binfile)
        field['Ny'].astype(dtype=np.int32).tofile(binfile)
        field['Nz'].astype(dtype=np.int32).tofile(binfile)
        field['thetaground'].astype(dtype=np.float64).tofile(binfile)

        if spectral:
            field['uu'].reshape((field['uu'].size,1), order='F').tofile(binfile)
            field['vv'].reshape((field['vv'].size,1), order='F').tofile(binfile)
            field['ww'].reshape((field['ww'].size,1), order='F').tofile(binfile)

        else:
            # Convert to wavenumber space
            uu =fft.r2c(field['u'], field['Nx2'], field['Ny']).astype(np.complex128)
            vv =fft.r2c(field['v'], field['Nx2'], field['Ny']).astype(np.complex128)
            ww =fft.r2c(field['w'], field['Nx2'], field['Ny']).astype(np.complex128)
            # Reshape and write to file
            uu.reshape((uu.size,1), order='F').tofile(binfile)
            vv.reshape((vv.size,1), order='F').tofile(binfile)
            ww.reshape((ww.size,1), order='F').tofile(binfile)



