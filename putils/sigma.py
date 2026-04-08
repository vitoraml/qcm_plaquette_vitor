"""
Project: Plots the cluster self-energy function in the edges of the Brillouin zone (exact)

Description:

Reads the spectral data from the e*.tsv files and generates the plot.

Authors:
- Vitor Assunção Moreira Lima (vitor.aml@hotmail.com)
-

Created: 2025-06-01
Last Modified: 2025-06-01
Version: 1.0.0

License: MIT

Dependencies:
- Python3
- numpy
- matplotlib
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
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

# Reads the components from the output files

w=np.arange(-1,1,0.01)

f1=open(r'e11.tsv','r')
f2=open(r'e12.tsv','r')
f3=open(r'e13.tsv','r')
f4=open(r'e14.tsv','r')

y1, y2, y3, y4 = [], [], [], []
g00, gp0, g0p, gpp = [], [], [], []

for line in f1:
    v=[float(n) for n in line.split()]
    y1.append(v[2])
f1.close()

for line in f2:
    v=[float(n) for n in line.split()]
    y2.append(v[2])
f2.close()

for line in f3:
    v=[float(n) for n in line.split()]
    y3.append(v[2])
f3.close()

for line in f4:
    v=[float(n) for n in line.split()]
    y4.append(v[2])
f4.close()

# Calculates the momentum quantities
for j in range(0,len(y1)):
    g00.append(y1[j]+y2[j]+y3[j]+y4[j])
    gp0.append(y1[j]-y2[j]+y3[j]-y4[j])
    g0p.append(y1[j]+y2[j]-y3[j]-y4[j])
    gpp.append(y1[j]-y2[j]-y3[j]+y4[j])
    

# Plots the figures

fig, ax=plt.subplots()

plt.xlim([-1,1])
plt.ylim([-6.0,0.0])
plt.xlabel('$\omega$')
plt.ylabel('$\Sigma$[$k_x$,$k_y$]($\omega$)')
plt.axvline(0,c='k',ls='dotted')
cmap = mpl.colormaps['inferno']

ax.plot(w,g00,c='k',ls='-', lw=2, label=r'[0,0]')
ax.plot(w,gp0,c='r',ls='-', lw=2, label=r'[$\pi$,0]')
ax.plot(w,g0p,c='b',ls='--', lw=2, label=r'[0,$\pi$]')
ax.plot(w,gpp,c='g',ls='dotted', lw=2, label=r'[$\pi$,$\pi$]')

plt.legend(loc='upper right',fontsize=7)
plt.savefig('sigma_cluster_bz.png', dpi=300)
plt.show()
