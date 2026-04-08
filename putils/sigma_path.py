"""
Project: Plots the periodized self-energy in a Brillouin zone k-path

Description:

Reads the spectral data from the sigma.tsv file and generates the plot.

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
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as patches
from matplotlib.colors import hsv_to_rgb
#================================================================================
# Reorganize the size of the labels 
def set_size(w,h, ax=None):
    """ w, h: width, height in inches """
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)

# Opens the full spectral file and defines the grids (same size as in the .py file)

f=open(r'sigma.tsv','r')

spec = [[],[],[],[]]

w=np.arange(-1,1,0.01)

nk = 201

kx = ky = np.linspace(0,np.pi,int(nk/2+1))

k = np.linspace(0,np.pi,int(nk/2+1))

offset=0.08

# Data analysis

# Read spectral function

for line in f:
    v=[float(n) for n in line.split()]
    spec[0].append(v[0])     # w
    spec[1].append(v[1])     # kx
    spec[2].append(v[2])     # ky
    spec[3].append(-v[3])     # spec
f.close()

# X-M direction

xm = []

# Filters the k-points in the X-M path and writes to a list
for i in range(0,len(w)*len(kx)*len(ky)):
    if spec[1][i] == 0.5:
        xm.append(spec[3][i])

# Rewrite the list as a numpy multidimensional array
xm = np.reshape(xm,(len(w),int(nk/2+1)))

# Plots the spectral function frequency dependence at every k-point in the path

fig, ax=plt.subplots()
plt.xlim([-1,1])
#plt.ylim([0.0,1.0])
plt.xlabel('$\omega$')
plt.ylabel('$\Sigma$($\omega$)')
plt.axvline(0,c='k',ls='dotted')

n=0

for j in range(0,int(nk/2+1)):
    xm_p=[]
    for i in range(0,len(w)):
        xm_p.append(xm[i][j]+offset*n)
    ax.plot(w,xm_p,c='k',ls='-', lw=1.5)
    n = n + 1


set_size(6,4)
plt.tight_layout()
#plt.legend(loc='upper right',fontsize=7)
plt.savefig('xm_sigma.png', dpi=300)
plt.show()

# For the other directions the process is the same

# Y-M direction

ym = []

for i in range(0,len(w)*len(kx)*len(ky)):
    if spec[2][i] == 0.5:
        ym.append(spec[3][i])

ym = np.reshape(ym,(len(w),int(nk/2+1)))

# Plot

fig, ax=plt.subplots()
plt.xlim([-1,1])
#plt.ylim([0.0,1.0])
plt.xlabel('$\omega$')
plt.ylabel('$\Sigma$($\omega$)')
plt.axvline(0,c='k',ls='dotted')

n=0

for j in range(0,int(nk/2+1)):
    ym_p=[]
    for i in range(0,len(w)):
        ym_p.append(ym[i][j]+offset*n)
    ax.plot(w,ym_p,c='k',ls='-', lw=1.5)
    n = n + 1


set_size(6,4)
plt.tight_layout()
#plt.legend(loc='upper right',fontsize=7)
plt.savefig('ym_sigma.png', dpi=300)
plt.show()

# G-M

gm = []

for i in range(0,len(w)*len(kx)*len(ky)):
    if spec[1][i] == spec[2][i]:
        gm.append(spec[3][i])

gm = np.reshape(gm,(len(w),int(nk/2+1)))

# Plot

fig, ax=plt.subplots()
plt.xlim([-1,1])
#plt.ylim([0.0,1.0])
plt.xlabel('$\omega$')
plt.ylabel('$\Sigma$($\omega$)')
plt.axvline(0,c='k',ls='dotted')

n=0

for j in range(0,int(nk/2+1)):
    gm_p=[]
    for i in range(0,len(w)):
        gm_p.append(gm[i][j]+offset*n)
    ax.plot(w,gm_p,c='k',ls='-', lw=1.5)
    n = n + 1


set_size(6,4)
plt.tight_layout()
#plt.legend(loc='upper right',fontsize=7)
plt.savefig('gm_sigma.png', dpi=300)
plt.show()

# X-Y

xy = []

for i in range(0,len(w)*len(kx)*len(ky)):
    if (spec[1][i]+spec[2][i]) == 0.5:
        xy.append(spec[3][i])

xy = np.reshape(xy,(len(w),int(nk/2+1)))

# Plot

fig, ax=plt.subplots()
plt.xlim([-1,1])
#plt.ylim([0.0,1.0])
plt.xlabel('$\omega$')
plt.ylabel('$\Sigma$($\omega$)')
plt.axvline(0,c='k',ls='dotted')

n=0

for j in range(0,int(nk/2+1)):
    xy_p=[]
    for i in range(0,len(w)):
        xy_p.append(xy[i][j]+offset*n)
    ax.plot(w,xy_p,c='k',ls='-', lw=1.5)
    n = n + 1


set_size(6,4)
plt.tight_layout()
#plt.legend(loc='upper right',fontsize=7)
plt.savefig('xy_sigma.png', dpi=300)
plt.show()
