#! /usr/bin/env python
import numpy as np
import load_sp as lsp
import os
import sys

time_vect = np.linspace(0, 1.8, 601)

Nt = time_vect.size
Nx = 384; Lx = 10
Ny = 256; Ly = 3.6
Nz = 144; Lz = 1
Nel = Nx*Ny*Nz

# Flag to check whether initial files are double or single precision
input_double = False

################################
# Create the vdf file
################################
vdf_name = 'windfarm_opt.vdf'
vdf_createstring = 'vdfcreate -numts '+str(len(time_vect))+' -dimension '+str(Nx)+'x'+str(Ny)+'x'+str(Nz)+' -extents 0:0:0:'+str(Lx)+':'+str(Ly)+':'+str(Lz)+' -vars3d u:vort '+vdf_name
print(vdf_createstring)
os.system(vdf_createstring)


for ind, t in enumerate(time_vect):

    filename_read_velocity  = 'velocity_field_t_'+"{:4.4f}".format(t)+'.dat'
    filename_read_vorticity = 'vorticity_mag_t_'+"{:4.4f}".format(t)+'.dat'
    filename_writeu = 'vapor_u_t_'+"{:4.4f}".format(t)+'.raw'
    filename_writev = 'vapor_v_t_'+"{:4.4f}".format(t)+'.raw'
    filename_writew = 'vapor_w_t_'+"{:4.4f}".format(t)+'.raw'
    filename_writevort = 'vapor_vort_t_'+"{:4.4f}".format(t)+'.raw'

    print('Reading file: ', filename_read_velocity)

    if input_double:
        ################################
        # Read from double precision files
        ################################
        # Velocity
        with open(filename_read_velocity, 'rb') as binfile_read:
            u = np.fromfile(binfile_read, dtype=np.float32, count=Nel)

        # Vorticity
        with open(filename_read_vorticity, 'rb') as binfile_read:
            vort = np.fromfile(binfile_read)

        ################################
        # Write to new file using single precision
        ################################
        # Velocity
        with open(filename_writeu, 'wb') as binfile_write:
            u.astype(np.float32).tofile(binfile_write)
        with open(filename_writev, 'wb') as binfile_write:
            v.astype(np.float32).tofile(binfile_write)
        with open(filename_writew, 'wb') as binfile_write:
            w.astype(np.float32).tofile(binfile_write)

        # Vorticity
        with open(filename_writevort, 'wb') as binfile_write:
            vort.astype(np.float32).tofile(binfile_write)
    else:
        filename_writeu = filename_read_velocity
        filename_writevort = filename_read_vorticity

    ################################
    # Populate the vdf file        #
    ################################
    os.system('raw2vdf -ts '+str(ind)+' -varname u '+vdf_name+' '+filename_writeu)
    os.system('raw2vdf -ts '+str(ind)+' -varname v '+vdf_name+' '+filename_writev)
    os.system('raw2vdf -ts '+str(ind)+' -varname w '+vdf_name+' '+filename_writew)
    os.system('raw2vdf -ts '+str(ind)+' -varname vort '+vdf_name+' '+filename_writevort)

    # Remove the vapor raw file
    if input_double:
        os.system('rm '+filename_writeu)
        os.system('rm '+filename_writevort)
    print('----------------------------------------------------------------------------')
