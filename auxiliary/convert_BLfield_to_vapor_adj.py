import numpy as np
import load_sp as lsp
import os
import ops_sp as osp

# INPUT PARAMETERS DEFINING DOMAIN ETC.
#-----------------------------------------------------
scratchdir='/scratch/leuven/306/vsc30627/'
paths = [
'/scratch/leuven/306/vsc30627/OPTIMIZATION_YAW/turbulent_T_0.3/arch_opt/window1/run_aligned_induction_ct2/'
    ]
vdfpaths = [
    scratchdir + 'vdf_adjoint_turb/',
    ]

for path, vdf_path in zip(paths, vdfpaths):

    if not os.path.exists(vdf_path):
        os.mkdir(vdf_path)
    #path = '/scratch/leuven/306/vsc30627/OPTIMIZATION_YAW/laminar_T_0.3/dirs_fwd/run_aligned_yaw_0_3/'
    #vdf_path = '/scratch/leuven/306/vsc30627/vdf_files_lam_0_3/'
    setup = lsp.setup(path=path)
    filename_bl = path+'BL_field.dat' 
    filename_vdf= vdf_path+'vdf_bl.vdf'
    
    Lx, Ly, Lz = setup.Lx, setup.Ly, 1.
    Nx, Ny, Nz = setup.Nx2, setup.Ny, setup.zmesh_st.size
    
    varnames = 'u:v:w:vort:vortmag:xi1:xi2:xi3'
    
    
    
    # CREATE THE VDF FILE
    #-----------------------------------------------------
    vdf_createstring = 'vdfcreate -numts 1 -dimension '+str(Nx)+'x'+str(Ny)+'x'+str(Nz) \
                    +' -extents 0:0:0:'+str(Lx)+':'+str(Ly)+':'+str(Lz)+' -vars3d '+varnames+' '+filename_vdf
    print(vdf_createstring)
    os.system(vdf_createstring)
    
    # LOAD DATA & CONVERT DATA TO RAW FILES
    #-----------------------------------------------------
    bl = lsp.load_BLfield_real(filename_bl, setuppath=path)
    adj = lsp.load_BLfield_real(path+'adjoint_field.dat', setuppath=path)
    #vort = np.linalg.norm(osp.calc_vorticity(bl),axis=3)
    vort = osp.calc_vorticity(bl)
    vortmag = np.linalg.norm(vort, axis=3)
    vort = vort[:,:,:,2]
    datanames = ['u','v','w','vort', 'vortmag','xi1','xi2','xi3']
    dataset = [bl['u'],bl['v'],bl['w'],vort,vortmag, adj['u'], adj['v'], adj['w']]
    filenames_raw = ['u_raw', 'v_raw', 'w_raw', 'vort_raw', 'vortmag_raw', 'xi1r', 'xi2r', 'xi3r']
    
    for data, file in zip(dataset, filenames_raw):
        print('Writing to '+file)
        with open(scratchdir+file, 'wb') as binfile_write:
            write_data = data.flatten(order='F')
            write_data.astype(np.float32).tofile(binfile_write)
    
    # POPULATE THE VDF FILE
    #-----------------------------------------------------
    for varname, file in zip(datanames, filenames_raw):
        print('raw2vdf -ts 0 -varname '+varname+' '+filename_vdf+' '+scratchdir+file)
        os.system('raw2vdf -ts 0 -varname '+varname+' '+filename_vdf+' '+scratchdir+file)
        os.system('rm '+scratchdir+file)


