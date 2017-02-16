import load_sp as lsp
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

class Windfarm:
    def __init__(self, Nrows=0, Ncols=0, path='./', windpowerfile='Windpower.dat',windfarmfile='windfarm.setup', controlfile='control.dat', timestep=.00075, ct_filt_init=2, load_measurements=False, tau=-1, verbose=True):
        path                             = path+'/'
        self.path                        = path
        if Nrows==0 or Ncols==0:
            Nrows, Ncols                 = countrows(path+windfarmfile, verbose)
        self.Nrows                       = Nrows
        self.Ncols                       = Ncols
        self.turbines, self.timeconstant = load_windfarm(path+windfarmfile, self.Nrows, self.Ncols)
        if not tau == -1:
            self.timeconstant=tau/1000
        self.wptime, self.windpower, self.thrustforce = load_windpower(path+windpowerfile, self.Nrows, self.Ncols)
        self.meanwindpower               = np.mean(self.windpower, 2)
        self.totalpower                  = np.sum(np.sum(self.windpower,0),0)
        self.power_per_row               = np.mean(self.meanwindpower,1) 
        self.power_per_col               = np.mean(self.meanwindpower,0) 
        self.rows                        = np.arange(Nrows)
        self.cols                        = np.arange(Ncols)
        self.time_integrated_power       = calc_time_integrated_power(self.totalpower, self.wptime)

        # Some stuff on controls
        self.controls                    = load_controls(path+controlfile, self.Nrows, self.Ncols)
        self.controls_filt               = time_filter_controls(self.controls, self.timeconstant, timestep, self.Nrows, self.Ncols, path=path, ct_filt_init=2)
        #self.diskvel                     = self.windpower/self.thrustforce

        # Load the disk velocity file if it exists
        self.diskvelocity               = load_upstream_velocity_file(path+'Velocity_upstream_0D.dat', self.Nrows, self.Ncols)

        # Upstream measurements (if existing)
        if load_measurements:
            self.upstream_vel                = load_upstream_vel(path, self.Nrows, self.Ncols)
            self.upstream_shear_0D           = load_upstream_velocity_file(path+'Shear_upstream_0D.dat', self.Nrows, self.Ncols)
            self.upstream_vel_top_0D         = load_upstream_velocity_file(path+'Velocity_upstream_top_0D.dat', self.Nrows, self.Ncols)
            self.upstream_vel_bot_0D         = load_upstream_velocity_file(path+'Velocity_upstream_bot_0D.dat', self.Nrows, self.Ncols)
            self.upstream_vel_left_0D        = load_upstream_velocity_file(path+'Velocity_upstream_left_0D.dat', self.Nrows, self.Ncols)
            self.upstream_vel_right_0D       = load_upstream_velocity_file(path+'Velocity_upstream_right_0D.dat', self.Nrows, self.Ncols)
            self.upstream_shear_1D           = load_upstream_velocity_file(path+'Shear_upstream_1D.dat', self.Nrows, self.Ncols)
            self.upstream_vel_top_1D         = load_upstream_velocity_file(path+'Velocity_upstream_top_1D.dat', self.Nrows, self.Ncols)
            self.upstream_vel_bot_1D         = load_upstream_velocity_file(path+'Velocity_upstream_bot_1D.dat', self.Nrows, self.Ncols)
            self.upstream_vel_left_1D        = load_upstream_velocity_file(path+'Velocity_upstream_left_1D.dat', self.Nrows, self.Ncols)
            self.upstream_vel_right_1D       = load_upstream_velocity_file(path+'Velocity_upstream_right_1D.dat', self.Nrows, self.Ncols)

            self.reynolds_1D                 = load_upstream_Reynolds_file(path+'Reynolds_upstream_1D.dat',       self.Nrows, self.Ncols)
            if not type(self.reynolds_1D) == 'int':
                self.tke_1D                  = 1/2*(self.reynolds_1D[:,:,:,0,0] + self.reynolds_1D[:,:,:,1,1] + self.reynolds_1D[:,:,:,2,2])
            self.reynolds_left_1D            = load_upstream_Reynolds_file(path+'Reynolds_upstream_left_1D.dat',  self.Nrows, self.Ncols)
            self.reynolds_right_1D           = load_upstream_Reynolds_file(path+'Reynolds_upstream_right_1D.dat', self.Nrows, self.Ncols)
            self.reynolds_top_1D             = load_upstream_Reynolds_file(path+'Reynolds_upstream_top_1D.dat',   self.Nrows, self.Ncols)
            self.reynolds_bot_1D             = load_upstream_Reynolds_file(path+'Reynolds_upstream_bot_1D.dat',   self.Nrows, self.Ncols)

            self.reynolds_0D                 = load_upstream_Reynolds_file(path+'Reynolds_upstream_0D.dat',       self.Nrows, self.Ncols)
            if not type(self.reynolds_0D) == 'int':
                self.tke_0D                  = 1/2*(self.reynolds_0D[:,:,:,0,0] + self.reynolds_0D[:,:,:,1,1] + self.reynolds_0D[:,:,:,2,2])
            self.reynolds_left_0D            = load_upstream_Reynolds_file(path+'Reynolds_upstream_left_0D.dat',  self.Nrows, self.Ncols)
            self.reynolds_right_0D           = load_upstream_Reynolds_file(path+'Reynolds_upstream_right_0D.dat', self.Nrows, self.Ncols)
            self.reynolds_top_0D             = load_upstream_Reynolds_file(path+'Reynolds_upstream_top_0D.dat',   self.Nrows, self.Ncols)
            self.reynolds_bot_0D             = load_upstream_Reynolds_file(path+'Reynolds_upstream_bot_0D.dat',   self.Nrows, self.Ncols)

        self.Sx = 0.
        self.Sy = 0.


        

    def plot_power_per_row(self, show=True, lw=1, ls='-', mfc='k', mec='k', mew=1, xlabel='row', **kwargs):
        if 'normalization' in kwargs:
            normalization = kwargs['normalization']
        else:
            normalization = self.power_per_row[0]

        if 'label' in kwargs:
            plt.plot(self.rows+1, self.power_per_row/normalization, ls, lw=lw, label=kwargs['label'],mec=mec,mfc=mfc,mew=mew)
        else:
            plt.plot(self.rows+1, self.power_per_row/normalization, ls, lw=lw, mec=mec,mfc=mfc,mew=mew)
        plt.xlabel(xlabel)
        plt.ylabel(r'$P_i/P_{1}$')
        if show:
            plt.show()

        plt.xlim(0.5,self.Nrows+0.5)

    def plot_total_power(self, show=True, label='', normalization=1, xlabel='time', **kwargs):
        plot_power = self.windpower
        if 'timefactor' in kwargs:
            timefactor=kwargs['timefactor']
        else:
            timefactor=1

        if 'row' in kwargs and 'col' in kwargs:
            ylab = 'P row = '+str(kwargs['row']) + ' col = '+str(kwargs['col'])
            plot_power = plot_power[kwargs['row']-1, kwargs['col']-1, :]
        elif 'row' in kwargs:
            ylab = 'P row = '+str(kwargs['row'])
            plot_power = np.sum(plot_power[kwargs['row']-1, :, :],0)
        elif 'col' in kwargs:
            plot_power = np.sum(plot_power[:,kwargs['col']-1, :],0)
            ylab = 'P col = '+str(kwargs['col'])
        else:
            plot_power = np.sum( np.sum( plot_power, 0), 0)
            ylab = r'$P_{tot}$'

        if 'lw' in kwargs:
            linewidth = kwargs['lw']
        else:
            linewidth = 1

        if 'ls' in kwargs:
            linestyle = kwargs['ls']
        else:
            linestyle = '-'
        
        print(normalization)

        if 'color' in kwargs:
            plt.plot(self.wptime*timefactor, plot_power/normalization,linestyle+kwargs['color'], label=label, lw=linewidth)
        else:
            plt.plot(self.wptime*timefactor, plot_power/normalization,linestyle, label=label, lw=linewidth)
        plt.xlabel(xlabel)
        plt.ylabel(ylab)
        if show:
            plt.show()

    def plot_turbines(self, show=True,color='k',transpose=False, cols=np.arange(1000), how='topview'):
        for turbine in self.turbines:
            if how=='topview':
                if turbine.col in cols:
                    if transpose:
                        plt.plot((turbine.y-turbine.r, turbine.y+turbine.r),(turbine.x, turbine.x), color,lw=2)
                    else:
                        plt.plot((turbine.x, turbine.x), (turbine.y-turbine.r, turbine.y+turbine.r),color,lw=2)
            elif how=='sideview':
                if turbine.col == 0:
                    plt.plot((turbine.x, turbine.x), (turbine.H-turbine.r, turbine.H+turbine.r),color,lw=2)

            elif how=='frontview':
                print('Frontview not implemented yet in windfarm module.')

        if show:
            plt.show()

    def plot_turbines_yaw(self, index, color='k'):
        if not os.path.exists(self.path+'/Turbine_yaw.dat'):
            print('Yawing file not found...')

        else:
            yaw = np.loadtxt(self.path+'/Turbine_yaw.dat', skiprows=7)[:,1:]
            for it, turbine in enumerate(self.turbines):
                if it==0:
                    c = 'r'
                elif it==1:
                    c = 'g'
                elif it==2:
                    c = 'b'
                else:
                    c = 'k'
                angle = yaw[index, it]*np.pi/180
                plt.plot((turbine.x - turbine.r*np.sin(angle), turbine.x + turbine.r*np.sin(angle)),
                         (turbine.y + turbine.r*np.cos(angle), turbine.y - turbine.r*np.cos(angle)), c, lw=2)



class Turbine:
    def __init__(self, turb, row, col):
        self.x = turb[0]
        self.y = turb[1]
        self.H = turb[2]
        self.r = turb[4]
        self.row = row
        self.col = col


def countrows(filename, verbose):
    if os.path.exists(filename):
        Nturb = np.int(np.genfromtxt(filename, max_rows=1))
        dummyfarm = np.loadtxt(filename, skiprows=2)
        Ncols = np.size(np.extract(dummyfarm[:,0]==dummyfarm[0,0], dummyfarm[:,0]))
        Nrows = np.size(np.extract(dummyfarm[:,1]==dummyfarm[0,1], dummyfarm[:,1]))
        
        if not Nturb == Nrows*Ncols:
            print('Automatic row & column identification failed.')
            return 0,0
        else:
            if verbose:
                print('%%%%%%%%%%%%%%%%%%%%')
                print(' Nturb = '+str(Nturb))
                print(' Nrows = '+str(Nrows))
                print(' Ncols = '+str(Ncols))
                print('%%%%%%%%%%%%%%%%%%%%')
            return Nrows, Ncols
    else: 
        print('Windfarm file not found.', filename)
        return 0, 0
        

def load_windfarm(filename, Nrows, Ncols):
    if os.path.exists(filename):
        dummy = np.loadtxt(filename, skiprows=2)
        timeconstant = np.genfromtxt(filename, skip_header=1, usecols=2, max_rows=1)
    #    rows, cols = get_rows_and_cols(dummy[0,:], dummy[1,:], Nrows, Ncols)
    #    Nrows = rows.size
    #    Ncols = cols.size
    #    print 'Amount of rows    = ', Nrows
    #    print 'Amount of columns = ', Ncols
        windfarm = [] # Empty list
        if dummy.size==11:
            windfarm.append(Turbine(dummy, 1, 1))
        else:
            for index, turbine_data in enumerate(dummy):
                row, col = get_row_col(index, Ncols)
                windfarm.append(Turbine(turbine_data, row, col))
        return windfarm, timeconstant
    else:
        print( 'windfarm.setup not found', filename )
        return 0, 0

def load_controls(controlfile, Nrows, Ncols):
    if os.path.exists(controlfile):
        control_dummy = np.loadtxt(controlfile, skiprows=5)[:,1:]
        Nt = control_dummy.shape[0]
        controls = np.zeros((Nrows, Ncols, Nt)) 
        for num in range(Nrows*Ncols):
            row, col = get_row_col(num, Ncols)
            controls[row, col, :] = control_dummy[:, num]
        return controls
    else:
        print( 'Controls file not found.', controlfile)
        return np.zeros((Nrows, Ncols, 1))

def load_windpower(filename, Nrows, Ncols):
    if os.path.exists(filename):
        dummy = np.loadtxt(filename,comments='%')
        time = dummy[:,0]
        powerdum = -dummy[:,2::2]
        forcedum = -dummy[:,1::2]
        # Reshape data into windfarm format, this assumes the turbines are numbered row per row
        power = np.zeros((Nrows, Ncols, time.size))
        force = np.zeros((Nrows, Ncols, time.size))
        for num in range(Nrows*Ncols):
            row, col = get_row_col(num, Ncols)
            power[row, col, :] = powerdum[:, num]
            force[row, col, :] = forcedum[:, num]

        # Now we remove the zeros that are written at the start of every simulation..
        indices = []
        for i in range(time.size):
            if np.all( power[:,:,i] == 0):
                indices.append(i)

        time = np.delete(time, indices)
        power = np.delete(power, indices, axis=2)
        force = np.delete(force, indices, axis=2)
        

        return time, power, force
    else:
        print( 'Windpower file not found.' , filename)
        return np.zeros((1)), np.zeros((Nrows, Ncols, 1)), np.zeros((Nrows, Ncols, 1))

def get_row_col(index, Ncols):
    row = int(index//Ncols)
    col = index - row*Ncols
    return row, col

def calc_time_integrated_power(power, time):
    J = 0
    for i in range(1, time.size):
        J += power[i]*(time[i] - time[i-1])
    return J

def load_upstream_Reynolds_file(filename, Nrows, Ncols):
    if os.path.exists(filename):
        dum = np.loadtxt(filename, skiprows=1)[:,1:]
        uu = dum[:,0::6]
        vv = dum[:,1::6]
        ww = dum[:,2::6]
        uv = dum[:,3::6]
        uw = dum[:,4::6]
        vw = dum[:,5::6]
    else:
        return 0
    Ntime = uu.shape[0]
    Nturb = uu.shape[1]

    reynolds = np.zeros((Nrows, Ncols, Ntime, 3, 3))
    for row in range(Nrows):
        for col in range(Ncols):
            turb = row*Ncols+col
            reynolds[row, col, :, 0, 0] = uu[:,turb]
            reynolds[row, col, :, 1, 1] = vv[:,turb] 
            reynolds[row, col, :, 2, 2] = ww[:,turb]
            reynolds[row, col, :, 0, 1] = uv[:,turb]
            reynolds[row, col, :, 0, 2] = uw[:,turb]
            reynolds[row, col, :, 1, 2] = vw[:,turb]
    return reynolds


def load_upstream_velocity_file(filename, Nrows, Ncols):
    if os.path.exists(filename):
        dum = np.loadtxt(filename, skiprows=1)[:,1:]
        u = dum[:,0::3]
        v = dum[:,1::3]
        w = dum[:,2::3]
    else:
        return 0
    Ntime = u.shape[0]
    Nturb = u.shape[1]
                
    vel = np.zeros((Nrows, Ncols, Ntime, 3))
    for row in range(Nrows):
        for col in range(Ncols):
            turb = (row)*Ncols+col
            vel[row, col, :, 0]= u[:,turb]
            vel[row, col, :, 1]= v[:,turb]
            vel[row, col, :, 2]= w[:,turb]
    return vel

def load_upstream_vel(path, Nrows, Ncols):
    filename = path+'Velocity_upstream_0D.dat' 
    if os.path.exists(filename):
        dum_0 = np.loadtxt(filename, skiprows=1)[:,1:]
        u_0 = dum_0[:,0::3]
        v_0 = dum_0[:,1::3]
        w_0 = dum_0[:,2::3]
    else:
        return 0
    filename = path+'Velocity_upstream_1D.dat' 
    if os.path.exists(filename):
        dum_1 = np.loadtxt(filename, skiprows=1)[:,1:]
        u_1 = dum_1[:,0::3]
        v_1 = dum_1[:,1::3]
        w_1 = dum_1[:,2::3]
    else:
        return 0
    filename = path+'Velocity_upstream_2D.dat' 
    if os.path.exists(filename):
        dum_2 = np.loadtxt(filename, skiprows=1)[:,1:]
        u_2 = dum_2[:,0::3]
        v_2 = dum_2[:,1::3]
        w_2 = dum_2[:,2::3]
    else:
        return 0
    filename = path+'Velocity_upstream_3D.dat' 
    if os.path.exists(filename):
        dum_3 = np.loadtxt(filename, skiprows=1)[:,1:]
        u_3 = dum_3[:,0::3]
        v_3 = dum_3[:,1::3]
        w_3 = dum_3[:,2::3]
    else:
        return 0
    filename = path+'Velocity_upstream_4D.dat' 
    if os.path.exists(filename):
        dum_4 = np.loadtxt(filename, skiprows=1)[:,1:]
        u_4 = dum_4[:,0::3]
        v_4 = dum_4[:,1::3]
        w_4 = dum_4[:,2::3]
    else:
        return 0
    filename = path+'Velocity_upstream_5D.dat' 
    if os.path.exists(filename):
        dum_5 = np.loadtxt(filename, skiprows=1)[:,1:]
        u_5 = dum_5[:,0::3]
        v_5 = dum_5[:,1::3]
        w_5 = dum_5[:,2::3]
    else:
        return 0
    filename = path+'Velocity_upstream_6D.dat' 
    if os.path.exists(filename):
        dum_6 = np.loadtxt(filename, skiprows=1)[:,1:]
        u_6 = dum_6[:,0::3]
        v_6 = dum_6[:,1::3]
        w_6 = dum_6[:,2::3]
    else:
        return 0

    u = np.stack((u_0, u_1, u_2, u_3, u_4, u_5, u_6), axis=-1) 
    v = np.stack((v_0, v_1, v_2, v_3, v_4, v_5, v_6), axis=-1) 
    w = np.stack((w_0, w_1, w_2, w_3, w_4, w_5, w_6), axis=-1) 

    Ntime = u.shape[0]
    Nturb = u.shape[1]
    Nmeas = u.shape[2]
                
    upstream_vel = np.zeros((Nrows, Ncols, Ntime, Nmeas, 3))
    for row in range(Nrows):
        for col in range(Ncols):
            turb = (row)*Ncols+col
            upstream_vel[row, col, :, :, 0]= u[:,turb,:]
            upstream_vel[row, col, :, :, 1]= v[:,turb,:]
            upstream_vel[row, col, :, :, 2]= w[:,turb,:]

    return upstream_vel



def time_filter_controls(controls, timeconstant, timestep, Nrows, Ncols, path='./', ct_filt_init=2):
    if ct_filt_init == 0:
        # First read initial conditions
        if os.path.exists(path+'ct_filt_initial_cond.dat'):
            initial_condition = np.loadtxt(path+'ct_filt_initial_cond.dat', skiprows=1)[1:]
            initial_condition = np.reshape(initial_condition, (Nrows, Ncols), order='C')
        else:
            print('Initial conditions not present!')
            return 0
    else:
        # Every turbines initial condition is the one provided in ct_filt_init
        initial_condition = ct_filt_init*np.ones((Nrows, Ncols))

    # Perform the time integration
    Ntime = controls.shape[2]
    alfa = timestep/(timestep+timeconstant)
    controls_filt = np.zeros( controls.shape )

    # The first one...
    controls_filt[:,:,0] = (1-alfa)*initial_condition + alfa*controls[:,:,0]
    # And the rest...
    for t in range(1,Ntime):
        controls_filt[:,:,t] = (1-alfa)*controls_filt[:,:,t-1] + alfa*controls[:,:,t]
    
    return controls_filt
