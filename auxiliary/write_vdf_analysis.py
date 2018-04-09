#! /usr/bin/env python
import numpy as np
import load_sp as lsp
import ops_sp as ops
import os
import sys

time_vect = np.arange(start=0, stop=.6, step=.0015)

Nt = time_vect.size
Nx = 384; Lx = 10
Ny = 256; Ly = 3.6
Nz = 144; Lz = 1
Nel = Nx*Ny*Nz

# Directories
scratch_dir = os.environ['VSC_SCRATCH']+'/'
input_dir  = scratch_dir + '/OPTIMIZATION_THRUST_ANALYSE/run_aligned_Ct3_tau_0/instantaneous/'
output_dir = scratch_dir + '/OPTIMIZATION_THRUST_ANALYSE/run_aligned_Ct3_tau_0/vapor/'

if not os.path.exists(output_dir):
    os.mkdir(output_dir)

varnames = ['u', 'v', 'w', 'vort', 'l2', 'q']
varnames_string = ''
for varname in varnames:
    varnames_string += varname+':'
varnames_string = varnames_string[:-1]

################################
# Create the vdf file
################################
vdf_name = output_dir+'vis.vdf'
vdf_createstring = 'vdfcreate -numts '+str(len(time_vect))+' -dimension '+str(Nx)+'x'+str(Ny)+'x'+str(Nz)+' -extents 0:0:0:'+str(Lx)+':'+str(Ly)+':'+str(Lz)+' -vars3d '+varnames_string+' '+vdf_name
print(vdf_createstring)
os.system(vdf_createstring)

# Loop over timesteps
for ind, t in enumerate(time_vect):

    # velocity
    filename_read_velocity  = input_dir+'velocity_field_t_'+"{:4.4f}".format(t)+'.dat'
    with open(filename_read_velocity, 'rb') as binfile_read:
        vel = np.fromfile(binfile_read, dtype=np.float32)
        u = vel[:Nx*Ny*Nz]
        v = vel[Nx*Ny*Nz:2*Nx*Ny*Nz]
        w = vel[2*Nx*Ny*Nz:]
    with open(scratch_dir+'u.tmp','wb') as binfile_write:
        u.tofile(binfile_write)
    with open(scratch_dir+'v.tmp','wb') as binfile_write:
        v.tofile(binfile_write)
    with open(scratch_dir+'w.tmp','wb') as binfile_write:
        w.tofile(binfile_write)

    # vort
    filename_read_vorticity = input_dir+'vorticity_mag_field_t_'+"{:4.4f}".format(t)+'.dat'
    with open(filename_read_vorticity, 'rb') as binfile_read:
        vort = np.fromfile(binfile_read, dtype=np.float32)
    with open(scratch_dir+'vort.tmp','wb') as binfile_write:
        vort.tofile(binfile_write)

    # l2
    filename_read_l2 = input_dir+'lambda2_field_t_'+"{:4.4f}".format(t)+'.dat'
    with open(filename_read_l2, 'rb') as binfile_read:
        l2 = np.fromfile(binfile_read, dtype=np.float32)
    with open(scratch_dir+'l2.tmp','wb') as binfile_write:
        l2.tofile(binfile_write)

    #q
    filename_read_q = input_dir+'q_field_t_'+"{:4.4f}".format(t)+'.dat'
    with open(filename_read_q, 'rb') as binfile_read:
        q = np.fromfile(binfile_read, dtype=np.float32)
    with open(scratch_dir+'q.tmp','wb') as binfile_write:
        q.tofile(binfile_write)


    ################################
    # Populate the vdf file        #
    ################################
    os.system('raw2vdf -ts '+str(ind)+' -varname u '+vdf_name+' '+scratch_dir+'u.tmp')
    os.system('raw2vdf -ts '+str(ind)+' -varname v '+vdf_name+' '+scratch_dir+'v.tmp')
    os.system('raw2vdf -ts '+str(ind)+' -varname w '+vdf_name+' '+scratch_dir+'w.tmp')
    os.system('raw2vdf -ts '+str(ind)+' -varname vort '+vdf_name+' '+scratch_dir+'vort.tmp')
    os.system('raw2vdf -ts '+str(ind)+' -varname l2 '+vdf_name+' '+scratch_dir+'l2.tmp')
    os.system('raw2vdf -ts '+str(ind)+' -varname q '+vdf_name+' '+scratch_dir+'q.tmp')

    # Remove the tmp files
    os.system('rm '+scratch_dir+'u.tmp')
    os.system('rm '+scratch_dir+'v.tmp')
    os.system('rm '+scratch_dir+'w.tmp')
    os.system('rm '+scratch_dir+'vort.tmp')
    os.system('rm '+scratch_dir+'l2.tmp')
    os.system('rm '+scratch_dir+'q.tmp')
