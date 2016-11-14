import numpy as np
import matplotlib.pyplot as plt

z0 = 10**(-4)
zhub = 0.1

Nz = 144

uhub = 8.5
kappa = 0.4

utau = uhub*kappa/(np.log(10**3))

z = np.linspace(0, 1, num=Nz*2+1)
zcc = z[1::2]
u = utau/kappa*np.log(zcc/z0)

plt.close('all')
plt.figure()
plt.plot(u, zcc, 'r')
plt.xlabel('U [m/s]')
plt.ylabel('z/H')
plt.axhline(y=0.1)

filename = 'BL_INFLOW_144'
with open(filename, 'w') as file:
    file.write("# Inflow file with boundary layer inflow with uhub = 8.5m/s")
    file.write("# U            V            W")
    for u_k in u:
        ustr = "{:7.7f}".format(u_k)
        file.write(ustr+"   0.000000   0.000000\n")
