import load_sp as lsp
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

class Windfarm:
    def __init__(self, Nrows=1, Ncols=1, path='./', windpowerfile='Windpower.dat',windfarmfile='windfarm.setup', controlfile='control.dat', timestep=.00075, ct_filt_init=0):
        self.path                        = path+'/'
        self.Nrows                       = Nrows
        self.Ncols                       = Ncols
        self.turbines, self.timeconstant = load_windfarm(path+windfarmfile, self.Nrows, self.Ncols)
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
        self.controls_filt               = time_filter_controls(self.controls, self.timeconstant, timestep, self.Nrows, self.Ncols, path=path)
        #self.diskvel                     = self.windpower/self.thrustforce

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
            plt.plot(self.wptime, plot_power/normalization,linestyle+kwargs['color'], label=label, lw=linewidth)
        else:
            plt.plot(self.wptime, plot_power/normalization,linestyle, label=label, lw=linewidth)
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


class Turbine:
    def __init__(self, turb, row, col):
        self.x = turb[0]
        self.y = turb[1]
        self.H = turb[2]
        self.r = turb[4]
        self.row = row
        self.col = col


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
        print( 'windfarm.setup not found' )
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
        print( 'Controls file not found.')
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
        print( 'Windpower file not found.' )
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

def time_filter_controls(controls, timeconstant, timestep, Nrows, Ncols, path='./', ct_filt_init=0):
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

#def get_rows_and_cols(x, y):
#    # Works only for regular windfarms!
#    xt1 = x[0]
#    yt1 = y[0]
#    Nrows = y.count(y1) # aka how many times does the same spanwise location occur?
#    Ncols = x.count(x1)
#    if not Nrows*Ncols == x.size: 
#        print 'Windfarm seems to be slightly irregular, row and column classification will be based on spanwise locations'
#        Ncols = x.size/Nrows
#    rows = np.zeros(Nrows)
#    cols = np.zeros(Ncols)

