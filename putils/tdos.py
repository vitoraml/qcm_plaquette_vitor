"""
Project: Plots the TDOS associated to a CDMFT calculation

Description:

Reads the spectral data from a .tsv file and generates the plot.

Authors:
- Vitor Assunção Moreira Lima (vitor.aml@hotmail.com)
-

Created: 2025-05-05
Last Modified: 2025-06-01
Version: 1.0.0

License: MIT

Dependencies:
- Python3
- numpy
- matplotlib
"""
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#================================================================================
# Reorganize the size of the labels 
def set_size(w,h, ax=None):
    #""" w, h: width, height in inches """
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)

# Reads the density of states from a converged calculation
f1 = open(r'dos_converged.tsv','r')

x1, y1 = [], []

next(f1)
for line in f1:
    v1=[float(n) for n in line.split()]
    x1.append(v1[0])
    y1.append(v1[1])

f1.close()

# Plot

fig, ax=plt.subplots()
plt.xlim([-1,1])
plt.ylim([0.0,0.25])
plt.xlabel('$\omega$')
plt.ylabel('N($\omega$)')
plt.axvline(0,c='k',ls='dotted')
ax.plot(x1,y1,c='k',ls='-')
#legend = ax.legend()

set_size(5,2)
plt.tight_layout()
plt.savefig('dos.png')
plt.show()

