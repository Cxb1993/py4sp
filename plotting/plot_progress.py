import numpy as np
import matplotlib.pyplot as plt
import load_sp as lsp

tend = 2.5
dt = .5
t    = np.arange(0,tend,dt)

kappa = .4
z01 = 1.e-4
z02 = z01/2

NSp = lsp.load_NSp('NS_post1.dat')

plt.figure()
for tme in t:
    tstr = "{:7.7f}".format(tme)
    print 'Loading time = '+ tstr

    name = 'BL_instcc_t='+tstr+'.dat'
    BL = lsp.load_BLcc(name)

    plt.subplot(221)
    uhub = BL['U'][13]
    print 'Uhub = ', uhub
    print 'uu_m = ', np.mean(BL['uu'])
    eturb = np.sqrt(1./3.*(BL['uu']+BL['vv']+BL['ww']))/uhub
    plt.plot(BL['z'], eturb, label='t = '+tstr)
#    plt.plot(BL['z'], BL['uu'], label='t = '+tstr)
    plt.xlabel('z')
    plt.ylabel('uu')
    
    plt.subplot(222)
    plt.semilogx(BL['z'], BL['U'])
    plt.xlabel('z')
    plt.ylabel('U')

plt.subplot(223)
plt.plot(NSp['t'], NSp['U'],'k')
plt.xlabel('t')
plt.ylabel('Ub')
ax2 = plt.gca().twinx()
ax2.plot(NSp['t'], NSp['Eturb'],'r')
ax2.set_ylabel(r'$E_t$')
for tl in ax2.get_yticklabels():
    tl.set_color('r')

plt.subplot(224)
#plt.plot(NSp['t'], NSp['Eturb'],label=r'$E_t$')
plt.plot(NSp['t'], NSp['sf_mean_x'], label=r'$\tau_w$ (RHS)')
#plt.plot(NSp['t'], NSp['sf_Um'], label=r'$\tau_w$ (U)')
plt.xlabel('t')
plt.ylabel(r'$\tau_w$ (RHS)')
plt.ylim((-0.8, -1.1))

r = 0.08/2
plt.subplot(221)
# plt.plot(( 0.1+r, 0.1+r) , (0, .10) , ':k')
# plt.plot(( 0.1-r, 0.1-r) , (0, .10) , ':k')
plt.legend()

plt.subplot(222)
u1 = 1/kappa*np.log(BL['z']/z01)
u2 = 1/kappa*np.log(BL['z']/z02)
plt.semilogx(BL['z'], u1,':k')
plt.semilogx(BL['z'], u2,':k')
plt.tight_layout()


plt.show()
