import matplotlib.pyplot as plt
import numpy as np
import windfarm as wf
from sklearn import datasets
from sklearn import linear_model
import auxiliary as aux

plt.close('all')

dir = '../rerun/correlations_2/run_aligned_Ct3_tau_0/'
farm = wf.Windfarm(path=dir)

# Set up some parameters
incl_upstream_vel = True
incl_upstream_acc = False
incl_windpower    = False
Ncols   = 6 
Nrows   = 12
Nt = 2400
dt = farm.wptime[1] - farm.wptime[0]


for row in range(Nrows):
    print(' Training algorithm for row '+str(row))
    Nfeatures = 0
    if incl_upstream_vel:
        print(' Including upstream velocity in observations ')
        Nfeatures += 21
    if incl_upstream_acc:
        print(' Including upstream acceleration in observations ')
        Nfeatures += 21
    if incl_windpower:
        print(' Including windpower in observations ')
        Nfeatures += 72
    observations = np.zeros([Nt*Ncols, Nfeatures])
    
    # Fill the observations matrix
    # ----------------------------------
    for col in range(Ncols):
        acceleration = aux.calc_acceleration(farm.upstream_vel[row, col, :, :, :],dt)
        for t in range(Nt):
            index = t+Nt*col
            # Observations: 
            # 1. Upstream velocities: u,v,w at 0D, 1D, 2D, 3D, 4D, 5D, 6D in previous timestep (21 entries)
            if incl_upstream_vel:
                observations[index,0:21]   = farm.upstream_vel[row, col, t, :, :].flatten() # Note farm.upstream_vel starts from t=0!
            # 2. Upstream accelerations: du/dt, dv/dt, dw/dt at same locations (21 entries)
            if incl_upstream_acc:
                observations[index,21:42]  = acceleration[t,:,:].flatten()
            # 3. Power signal of entire windfarm in previous timestep (72 entries)
            if incl_windpower:
                observations[index,42:114] = farm.windpower[:,:,t-1].flatten() # Check if this is actually previous timestep..
    
    # Fill the target matrix
    # ----------------------------------
    target = np.zeros([Nt*Ncols])
    for t in range(Nt):
        for turb in range(Ncols):
            target[t+Nt*turb] = farm.controls[row, turb, t]
    
    # Define train and test subsets
    # ----------------------------------
    Nn = 200
    observations_train = observations[:-Nn,:]
    observations_test  = observations[-Nn:,:]
    target_train = target[:-Nn]
    target_test  = target[-Nn:]
    
    # Do some machine learning now...
    # ----------------------------------
    regr = linear_model.LinearRegression()
    
    # Train the model
    # ----------------------------------
    print('#-------------------------#')
    print('#     Training model      #')
    print('#-------------------------#')
    regr.fit(observations_train, target_train)
    
    # Output the results
    # ----------------------------------
    print('#-------------------------#')
    print('# Regression coefficients #')
    print('#-------------------------#')
    print(regr.coef_)
    
    print('#--------------------------#')
    print('# Regression error & score #')
    print('#--------------------------#')
    print(' MSE ')
    print(np.mean((regr.predict(observations_test)-target_test)**2))
    print(' Variance score ')
    print(regr.score(observations_test, target_test))
    
    # Plot the results 
    # ----------------------------------
    plt.figure()
    plt.title('Results of machine learning prediction: row '+str(row))
    plt.plot(regr.predict(observations_test), '-or', label='Prediction')
    plt.plot(target_test, '-ok', label='Actual')
    plt.legend(numpoints=1, loc=0).draw_frame(False)

