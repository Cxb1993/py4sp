#!/usr/bin/env python

"""
plot_adjointfield.py: Script to plot a topview of the adjoint field for a windfarm simulation

parameters are hard-coded for optimization cases run for APS conference
"""
import matplotlib.pyplot as plt
import numpy as np
import load_sp as lsp
import windfarm as wf

# Some initialization and loading stuff
setup = lsp.setup(path='./')
adjoint_field = lsp.load_BLfield_real('adjoint_field.dat')

plt.figure()
plt.imshow( np.transpose(np.fliplr(adjoint_field['u'][:,:,20]), extent=(0, setup.Lx, 0, setup.Ly), cmap='coolwarm' ))
