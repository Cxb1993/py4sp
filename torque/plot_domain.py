import windfarm as wf
import numpy as np
import matplotlib.pyplot as plt

fs = 12

WF = wf.Windfarm(Nrows=12, Ncols=6)

# Plotting
plt.close('all')
plt.figure()
# Plot turbines
WF.plot_turbines()
# Plot fringe region
plt.axvline(x=9, linestyle='--', color='r')


plt.gca().set_aspect('equal')
plt.gca().autoscale(tight=True)
plt.xlim((0, 10))
plt.ylim((0, 3.6))
plt.xlabel(r'$x$ [km]', size=fs)
plt.ylabel(r'$y$ [km]', size=fs)

# Add some text
for i in range(12):
    turb = WF.turbines[6*i]
    plt.text(x=turb.x, y=3.7, s='R'+str(i+1), horizontalalignment='center', size=fs-2)

turb_above  = WF.turbines[67]
turb_before = WF.turbines[66-6]
turb_last   = WF.turbines[66]

dx = 0.5
dy = .05
plt.annotate('', xy=(turb_last.x+dx-.35, turb_last.y-dy), xycoords='data',
        xytext=(turb_above.x+dx-.35, turb_above.y+dy), textcoords='data',
        arrowprops={'arrowstyle': '<->'})
plt.annotate('', xy=(turb_last.x+dy, turb_last.y+dx+.3), xycoords='data',
        xytext=(turb_before.x-dy, turb_before.y+dx+.3), textcoords='data',
        arrowprops={'arrowstyle': '<->'})
plt.text(x=(turb_last.x+turb_before.x)/2, y=turb_last.y+dx+.35, s='6$D$', horizontalalignment='center', verticalalignment='bottom')
plt.text(x=turb_last.x+dx+.12, y=.5*(turb_last.y+turb_above.y), s='6$D$', horizontalalignment='right', verticalalignment='center')



plt.savefig('simulation_domain.eps')
