import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erf, erfinv
plt.close('all')

z0 = 10**(-4)
zhub = 0.1
D = 0.1
H = 1
Hh = 1

Nz = 144

# Generate the base profile
uhub = 8.5
kappa = 0.4
utau = uhub*kappa/(np.log(10**3))
z = np.linspace(0, 1, num=Nz*2+1)
zcc = z[1::2]
ubase = utau/kappa*np.log(zcc/z0)

# Generate the gust profile
eps = 2*erfinv(.99)
alpha = zhub - D/2
beta  = zhub + D/2
deltag = 2*D
Ug = 5

# Model spatial behavior
fz = (zcc/zhub*np.exp(1-zcc/zhub))**2
#fz =  Ug*erf(eps*zcc/2/alpha)

# Model temporal behavior with a cosine
Tend = 1
time = np.linspace(0, Tend,num=100)
Tstart = 0.1
Tstop  = 0.5
gt = np.zeros(time.size)

for i, ti in enumerate(time):
    if ti>=Tstart and ti<=Tstop:
        gt[i] = 1/2*(1 - np.cos(2*np.pi*(ti-Tstart)/(Tstop-Tstart)))

# Put space and time together
ugust = np.zeros((zcc.size, time.size))
for i, gi in enumerate(gt):
    ugust[:,i] = Ug*fz*gt[i]

ucorr = np.zeros(ugust.shape)
hz = np.exp(-zcc/Hh)-1
for i, ti in enumerate(time):
    sumh = np.sum(hz)
    sumgust = np.sum(ugust[:,i])
    Uc = sumgust/sumh
    ucorr[:,i] = -Uc*hz

plt.figure()
for i in range(time.size):
    plt.clf()
    plt.subplot(131)
    u = ubase + ugust[:,i] + ucorr[:,i]
    plt.plot(u, zcc,'b')
    plt.plot(ubase, zcc,'k')
    plt.xlim((4, 16))
    plt.xlabel('U [m/s]')
    plt.ylabel('z/H')
    plt.axhline(y=0.05, linestyle=':', color='k')
    plt.axhline(y=0.15, linestyle=':', color='k')

    plt.subplot(132)
    plt.plot(ugust[:,i],zcc,'r')
    plt.xlim((-Ug,Ug))
    plt.xlabel('U_gust [m/s]')
    plt.ylabel('z/H')
    plt.axhline(y=0.05, linestyle=':', color='k')
    plt.axhline(y=0.15, linestyle=':', color='k')

    plt.subplot(133)
    plt.plot(ucorr[:,i],zcc,'g')
    plt.xlim((-Ug,Ug))
    plt.xlabel('U_corr [m/s]')
    plt.ylabel('z/H')
    plt.suptitle('T = '+ str(time[i]))
    plt.axhline(y=0.05, linestyle=':', color='k')
    plt.axhline(y=0.15, linestyle=':', color='k')
    plt.pause(0.05)
    print(np.sum(u))

