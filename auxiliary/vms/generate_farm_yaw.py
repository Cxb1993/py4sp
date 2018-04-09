"""
    Script to generate vms file with non-rotating disks at turbine locations 
"""
import numpy as np
import windfarm as wf

dirs = [
        '/scratch/leuven/306/vsc30627/OPTIMIZATION_YAW/laminar_T_0.3/dirs_fwd/run_aligned_yaw/',
        '/scratch/leuven/306/vsc30627/OPTIMIZATION_YAW/laminar_T_0.3/dirs_fwd/run_aligned_yaw_0_3/',
        '/scratch/leuven/306/vsc30627/OPTIMIZATION_YAW/laminar_T_0.3/dirs_fwd/run_aligned_induction_ct2/',
        '/scratch/leuven/306/vsc30627/OPTIMIZATION_YAW/laminar_T_0.3/dirs_fwd/run_aligned_induction_ct2_yaw_0_3/',
        '/scratch/leuven/306/vsc30627/OPTIMIZATION_YAW/laminar_T_0.3/dirs_fwd/ref/',
        '/scratch/leuven/306/vsc30627/OPTIMIZATION_YAW/turbulent_T_0.3/dirs_fwd/run_aligned_yaw/',
        '/scratch/leuven/306/vsc30627/OPTIMIZATION_YAW/turbulent_T_0.3/dirs_fwd/run_aligned_yaw_0_3/',
        '/scratch/leuven/306/vsc30627/OPTIMIZATION_YAW/turbulent_T_0.3/dirs_fwd/run_aligned_induction_ct2/',
        '/scratch/leuven/306/vsc30627/OPTIMIZATION_YAW/turbulent_T_0.3/dirs_fwd/run_aligned_induction_ct2_yaw_0_3/',
        '/scratch/leuven/306/vsc30627/OPTIMIZATION_YAW/turbulent_T_0.3/dirs_fwd/ref/',
        ]
names = [
        'lam_yaw',
        'lam_yaw_0_3',
        'lam_ind',
        'lam_yawind',
        'lam_ref',
        'turb_yaw',
        'turb_yaw_0_3',
        'turb_ind',
        'turb_yawind',
        'turb_ref',
        ]

for path, outputfile in zip(dirs, names):
    modelfile  = "Disco_antiwarping-10mmx0.6mm.stl"
    
    # Expand this to read from windfarm file...
    farm = wf.Windfarm(path=path)
    if 'lam' in path:
        z = 0.5
    else:
        z = 0.1
    sc = 0.010
    
    yaw_angles = np.loadtxt(path+'/Turbine_yaw.dat', skiprows=7)[-1,1:]
    #yaw_angles = yaw_angles*0
    
    with open(outputfile+'.vms','w') as file:
        
        # Write the header
        file.write('<?xml version="1.0" encoding="ISO-8859-1" standalone="yes"?>\n')
        file.write('\n')
        file.write('<ModelScene>\n')
        file.write('\n')
        file.write('<File>'+modelfile+'</File>\n')
        file.write('\n')
        file.write('<!-- ts 0 -->')
    
        file.write('<Timestep>\n')
        file.write('\n')
        cnt = 0
        for yaw, turbine in zip(yaw_angles, farm.turbines):
    
            print('Writing turbine ', cnt)
            cnt+=1
    
            file.write('\n')
            file.write('<Transform>\n')
            file.write('<translate>'+"{:10.4f}".format(turbine.x)+"{:10.4f}".format(turbine.y)+"{:10.4f}".format(z)+'</translate>\n')
            file.write('<rotate> 0 1 0 90 </rotate>\n')
            file.write('<rotate> 1 0 0 '+str(-yaw)+' </rotate>\n')
            file.write('<scale>'+"{:10.4f}".format(sc)+"{:10.4f}".format(sc)+"{:10.4f}".format(sc)+'</scale>\n')
            file.write('<scale>0.5 0.5 0.5</scale>\n')
            file.write('</Transform>\n')
            file.write('\n')
    
        file.write('\n')
        file.write('</Timestep>\n')
        file.write('\n')
        file.write('</ModelScene>\n')
