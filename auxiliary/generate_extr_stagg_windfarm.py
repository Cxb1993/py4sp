import numpy as np
import matplotlib.pyplot as plt

def generate(dx_turb=.004, Ct_tower=1.2, 
        Ct=4./3., timeconstant=.005, powerconstant=0., turbinescaling=1000., 
        plot=False, staggered=False, **kwargs):


    # Domain parameters
    Lfringe = 1.
    Lx = 10.
    Ly = 3.6
    H = 1.
    
    # Turbine parameters
    D = .1
    zhub = .1
    
    # Farm parameters
    Nrows = 12
    Ncols = 10
    
    # Spacing
    if 'xturb' in kwargs:
        print('Taking turbine location from input values')
        Ncols = kwargs['Ncols']
        Nrows = kwargs['Nrows']
        xlocs = kwargs['xturb']
        ylocs = kwargs['yturb']
        Nturb = xlocs.size
        if not Nturb==Ncols*Nrows:
            print('Something wrong in input turbine locations!')
    else:
        print('Using default turbine locations')
        Nturb = Nrows*Ncols
        Sx = 3.6*D
        
        Sy = Ly/Ncols
        
        print('Sx/D = ', Sx/D)
        print('Sy/D = ', Sy/D)

        print(Sx/Sy)
        
        ylocs = np.zeros(Ncols)
        for col in range(Ncols):
            offset = Sy/2
            ylocs[col] = offset+col*Sy

        print('Y locations: ')
        print(ylocs)
        
        xlocs = np.zeros(Nrows)
#        offset = 1.5*Sx
        offset = .9
        for row in range(Nrows):
            xlocs[row] = offset+row*Sx
        
        print('X locations: ')
        print(xlocs)
    
    # Now we generate the file
    if staggered: 
        print('###################################')
        print('#  Generating extremely staggered windfarm  #')
        print('###################################')
    else: 
        print('###################################')
        print('#   Generating aligned windfarm   #')
        print('###################################')
    
    if plot:
        plt.figure()

    eps = 1e-6
    # Windfarm file
    if 'filename_farm' in kwargs:
        filenamef = kwargs['filename_farm']
    else:
        filenamef = 'windfarm.setup'
    with open(filenamef,'w') as windfarmfile:
        windfarmfile.write(str(Nturb)+'\n')
        windfarmfile.write("{:10.7f}".format(Ct) + "{:15.7f}".format(powerconstant) + "{:15.7f}".format(timeconstant) + "{:15.4f}".format(turbinescaling)+'\n')
        for turbine in range(Nturb):
            row = int(turbine//Ncols)
            col = turbine - row*Ncols
            xt = xlocs[row]
            yt = ylocs[col]
            if staggered:
                yt_1 = yt+Sy/3*row
                if col==4:
                    print('pre', yt_1, Ly)
                yt = yt_1 - int(10*yt_1)//int(10*Ly)*Ly
                if col==4:
                    print('post', yt, Ly)
#                if np.-1jjuabs(yt_1-Ly) <= eps:
#                    print(yt_1)
#                    print(Ly)
#                    print('--')
#                    yt = yt_1 - Ly
#                else:
#                    yt = yt_1
            windfarmfile.write("{:10.7f}".format(xt) + "{:15.7f}".format(yt) + "{:15.7f}".format(zhub) + "{:15.7f}".format(.015)+"{:15.7f}".format(D/2)+"{:15.7f}".format(20.)+"{:15.7f}".format(0.)+"{:15.7f}".format(1.5)+'\n')
            if plot:
                plt.plot((xt, xt), (yt-D/2, yt+D/2), 'k', lw=2)
    
    # Towerfarm file
    if 'filename_tower' in kwargs:
        filenamet = kwargs['filename_tower']
    else:
        filenamet = 'tower.setup'
    with open(filenamet, 'w') as towerfile:
        towerfile.write(str(Nturb)+'\n')
        towerfile.write("{:10.7f}".format(Ct_tower) + "{:15.7f}".format(timeconstant)+'\n')
        for tower in range(Nturb):
            row = int(tower//Ncols)
            col = tower - row*Ncols
            xt = xlocs[row]+dx_turb
            yt = ylocs[col]
            if staggered and not np.mod(row,2) == 0:
                yt = yt+Sy/2
                yt = yt - int(yt/Ly)*Ly
            towerfile.write("{:10.7f}".format(xt) + "{:15.7f}".format(yt) + "{:15.7f}".format(zhub) + "{:15.7f}".format(.0045000)+"{:15.7f}".format(zhub)+"{:15.7f}".format(1.5)+'\n')
            if plot:
                pass
#                plt.plot((xt, xt),(yt,yt),'or',mec='r')

    if plot:
        plt.ion()
        plt.axvline(x=9, linestyle='--', color='r')
        plt.gca().set_aspect('equal')
        plt.gca().autoscale(tight=True)
        plt.ylim((0,Ly))
        plt.xlim((0,Lx))
        plt.xlabel(r'$x$ [km]')
        plt.ylabel(r'$y$ [km]')
        plt.yticks((0,1,2,3))
        plt.xticks((0,1,2,3,4,5,6,7,8,9,10))
        plt.savefig('aligned.png',dpi=300)
        plt.savefig('aligned.eps')
        plt.show()



# Default generator...
# generate()
