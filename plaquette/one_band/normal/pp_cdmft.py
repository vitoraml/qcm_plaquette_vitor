### This program read the converged bath from cdmft.py output and generates a more complet set of output information

import numpy as np
import pandas as pd
import os
import pyqcm
from model import model
from pyqcm.cdmft import CDMFT
#================================================================================
pyqcm.set_global_parameter('max_iter_lanczos', 4000)
pyqcm.set_global_parameter('Ground_state_method','P')

var=[]
for j in range(1,2):
    for i in range(1,9):
        var.append('eb{:d}_{:d}'.format(i,j))
        var.append('tb{:d}_{:d}'.format(i,j))

# Here we read the bath parameters from cdmft.tsv like in bath.py

parf = '../cdmft_converged.tsv'
bp = pd.read_csv(parf, sep='\t')

eb = bp.filter(like='eb').iloc[:, 0:16]
tb = bp.filter(like='tb').iloc[:, 0:32]

eb = eb[eb.columns[::2]]
tb = tb[tb.columns[::2]]

params = pd.concat([eb,tb],axis = 1)
names = np.array(params.columns, dtype=object)
values = np.array(params)[-1]
st = names[0] + ' = ' + str(values[0])+'\n'
index=0
for n in names[1:]:
    index += 1
    st += n + ' = ' + str(values[index])+'\n'
varia_values = st

band_params = """
U = 8
t = 1.0
t1 = 0.3
mu = 1.5
"""

model.set_target_sectors('R0:N12:S0')			# If we already know the converged sector we can inform it here directly
model.set_parameters(band_params + varia_values)

# Defines a real frequency grid
w=np.arange(-1,1,0.01)

# Defines a matsubara frequency grid
fwn=open(r'wn.tsv','w')
N_mats=100			    # Number of frequencies in the grid
b=100				    # Inverse temperature beta
wn=[]

# Generate the matsubara frequencies wn and write them to a file
for n in range(0,N_mats):       
    wn.append((2*n+1)*np.pi/b)
    fwn.write(str(wn[n])+'\n')	    

fwn.close()

# Defines a k-point grid in the BZ and an auxiliary grid in the first quadrant

nk = 201
kx = ky = np.linspace(-0.5,0.5,nk)  # In units of 2*pi
kx_red = ky_red = np.linspace(0.0,0.5,int(nk/2+1))

# One iteration CDMFT cycle just to generate again the instance with the previously converged bath
CDMFT(model, varia=var, accur=1*1e-6, convergence='self-energy', accur_bath=1*1e-10, depth=2, accur_dist=1e-10, maxiter = 137000, max_value=10000, max_function_eval=5000000000, alpha=0.0, initial_step=0.0, beta=100, wc=4, grid_type='sharp')
#================================================================================
### Here the output generation will be described in more detail

# Starts the instance and generate the cluster or lattice averages
I = pyqcm.model_instance(model)
I.Green_function_solve()
I.averages(pr=True)
I.GS_consistency(check_ground_state=True,threshold=0.0001)

## Spectral plots (basic outputs)

# Calculates the first BZ momentum distribution curve at a given frequency 'freq' and in a plane 'plane' and plots in a file 'file'
I.mdc(nk=200, eta=0.1, orb=None, spin_down=False, quadrant=False, opt='GF', k_perp=0, freq=0.0, max=None, plane='xy', size=1.0, band_basis=False, sym=None, file='mdc.pdf')

# Calculates the DOS in a -w to w window and plots in a file 'file'
I.plot_DoS(w=4, eta=0.1, sum=False, progress=True, colors=None, file='dos.pdf')

# Calculates the spectral feature 'opt' between -wmax and wmax in a given BZ path 'path' and plots in a file 'file' and writes the data points to 'file'.tsv. The spectral options are the spectral function (A) or the self-energy (self).
I.spectral_function(wmax=0.5, eta=0.01, path=([0.0,0.0,0.0],[1.0,1.0,0.0]), nk=60, orb=None, offset=8, opt='A', Nambu_redress=False, inverse_path=False, title='G-M', file='spectral_GM.pdf')
os.system('mv spectral_data.tsv spectral_GM.tsv')

# Generates the non-interacting FS of the lattice model at a given BZ plane 'plane' and plots in a file 'file' and a spectral_data file
I.Fermi_surface(nk=200, orb=None, zone=((0,0),1), plane='xy', k_perp=0.0, file='fs_lattice.pdf')

# Plots the eigenvalues of the inverse Green function in w=0 for a given plane 'plane' in a file 'file'
I.G_dispersion(nk=400, orb=0, period='G', contour=False, inv=False, zone=((0, 0), 1), k_perp=0.0, plane='xy', file='g_disp.pdf')

# Calculates the Green function average
G_avg=I.Green_function_average(clus=0,spin_down=False)
print()
print('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'+'\n')
print('Green function frequency average'+'\n')
print(G_avg)
print('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'+'\n')

# Calculates the Green function density
G_dens=I.Green_function_density(clus=0)
print()
print('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'+'\n')
print('Impurity density (Obtained through G) ='+str(G_dens)+'\n')
print('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'+'\n')

# Plots the Luttinger surface in a BZ plane 'plane' at w=0 and writes to a file 'file'
I.Luttinger_surface(nk=400, orb=1, zone=((0, 0), 1), k_perp=0, plane='xy', file='luttinger_surf.pdf')

# Calculates the Grand potential free energy of the system using the Potthoff_functional
Grand_pot = I.Potthoff_functional(hartree=None, file='potthoff.tsv', symmetrized_operator=None, consistency_check=True)
print()
print('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'+'\n')
print('Free energy (Potthoff) = '+str(Grand_pot)+'\n')
print('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'+'\n')

# Generates the hopping matrix
t_mat=I.cluster_hopping_matrix(clus=0,spin_down=False,full=0)
print()
print('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'+'\n')
print('t matrix'+'\n')
print(t_mat)
print('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'+'\n')

# Calculates the cluster spectral feature 'opt' between -wmax and wmax in a given BZ path 'path' and plots in a file 'file'. The spectral options are the Green funtion (None), the self-energy (self) or hybridization (hyb).
I.cluster_spectral_function(wmax=0.5, eta=0.01, clus=0, spin_down=False, imaginary=False, offset=8, opt=None, blocks=False, color='k', file='spectral_cluster.pdf')

# Generates the density matrix in a subsystem A associated to sites 'sites'
n_mat=I.density_matrix(clus=0,sites=[1])
print()
print('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'+'\n')
print('n matrix'+'\n')
print(n_mat[0])
print('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'+'\n')

# Calculates the dispersion relation (free electron) in a BZ plane 'plane' and saves in a file 'file'.
I.plot_dispersion(nk=200, orb=None, spin_down=False, contour=False, zone=((0, 0), 1), k_perp=0, plane='xy', file='disp.pdf')

# Calculates the potential energy of the system using Tr(E*G)
Pot_ene = I.potential_energy()
print()
print('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'+'\n')
print('Potential energy = '+str(Pot_ene)+'\n')
print('\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'+'\n')


# To get some quantities such as the cluster quantities, we need to define a parameter that will store the output from the instance that we requested and write it to a file 'by hand'. It's important to remember that when dealing with the cluster Green function and the cluster Self-energy, they come in a matricial form, so on needs to be careful when extracting it's components

# Generates the k-dependent quasiparticle weight from the lattice self-energy and writes to a file

f=open(r'qp_k.tsv','w')

for i in range(0,len(kx_red)):
    for j in range(0,len(ky_red)):
        qp = I.QP_weight(np.array([kx_red[i],ky_red[j],0.0]), eta=0.1, orb=1, spin_down=False)
        f.write(str(kx[i])+'        '+str(ky[j])+'        '+str(qp[0])+'\n')

f.close()

# Generates the k-dependent energy dispersion relation and writes to a file

f=open(r'disp.tsv','w')

for i in range(0,len(kx_red)):
    for j in range(0,len(ky_red)):
        disp = I.dispersion(np.array([kx_red[i],ky_red[j],0.0]), label=0, spin_down=False)[0]
        f.write(str(kx[i])+'        '+str(ky[j])+'        '+str(disp[0])+'\n')

f.close()

# Generates the k-dependent spectral gap and writes to a file

f=open(r'gap.tsv','w')

for i in range(0,len(kx_red)):
    for j in range(0,len(ky_red)):
        gap = I.gap(np.array([kx_red[i],ky_red[j],0.0]), orb=1, threshold=0.001)
        f.write(str(kx[i])+'        '+str(ky[j])+'        '+str(gap)+'\n')

f.close()

# Generates the Hybridization matrix in a frequency and momentum grid and writes to a file

# Real w 
v11=open(r'v11.tsv','w')
v12=open(r'v12.tsv','w')
v13=open(r'v13.tsv','w')
v14=open(r'v14.tsv','w')

for n in range(0,len(w)):
    for i in range(0,len(kx_red)):
        for j in range(0,len(ky_red)):
            V=I.V_matrix(w[n], np.array([kx_red[i],ky_red[j],0]), spin_down=False)
            V11 = V[0,0]
            V12 = V[0,1]
            V13 = V[0,2]
            V14 = V[0,3]
            v11.write(str(w[n])+'        '+str(kx_red[i])+'        '+str(ky_red[j])+'        '+str(np.real(np.array([V11]))[0])+'        '+str(np.imag(np.array([V11]))[0])+'\n')
            v12.write(str(w[n])+'        '+str(kx_red[i])+'        '+str(ky_red[j])+'        '+str(np.real(np.array([V12]))[0])+'        '+str(np.imag(np.array([V12]))[0])+'\n')
            v13.write(str(w[n])+'        '+str(kx_red[i])+'        '+str(ky_red[j])+'        '+str(np.real(np.array([V13]))[0])+'        '+str(np.imag(np.array([V13]))[0])+'\n')
            v14.write(str(w[n])+'        '+str(kx_red[i])+'        '+str(ky_red[j])+'        '+str(np.real(np.array([V14]))[0])+'        '+str(np.imag(np.array([V14]))[0])+'\n')

v11.close()
v12.close()
v13.close()
v14.close()

#Im iwn

v11=open(r'v11_mats.tsv','w')
v12=open(r'v12_mats.tsv','w')
v13=open(r'v13_mats.tsv','w')
v14=open(r'v14_mats.tsv','w')

for n in range(0,len(wn)):
    for i in range(0,len(kx_red)):
        for j in range(0,len(ky_red)):
            V=I.V_matrix(wn[n], np.array([kx_red[i],ky_red[j],0]), spin_down=False)
            V11 = V[0,0]
            V12 = V[0,1]
            V13 = V[0,2]
            V14 = V[0,3]
            v11.write(str(wn[n])+'        '+str(kx_red[i])+'        '+str(ky_red[j])+'        '+str(np.real(np.array([V11]))[0])+'        '+str(np.imag(np.array([V11]))[0])+'\n')
            v12.write(str(wn[n])+'        '+str(kx_red[i])+'        '+str(ky_red[j])+'        '+str(np.real(np.array([V12]))[0])+'        '+str(np.imag(np.array([V12]))[0])+'\n')
            v13.write(str(wn[n])+'        '+str(kx_red[i])+'        '+str(ky_red[j])+'        '+str(np.real(np.array([V13]))[0])+'        '+str(np.imag(np.array([V13]))[0])+'\n')
            v14.write(str(wn[n])+'        '+str(kx_red[i])+'        '+str(ky_red[j])+'        '+str(np.real(np.array([V14]))[0])+'        '+str(np.imag(np.array([V14]))[0])+'\n')

v11.close()
v12.close()
v13.close()
v14.close()

# Obtain the real frequency cluster Green's functions, self-energies and hybridizations

f11=open(r'g11.tsv','w')
f12=open(r'g12.tsv','w')
f13=open(r'g13.tsv','w')
f14=open(r'g14.tsv','w')
e11=open(r'e11.tsv','w')
e12=open(r'e12.tsv','w')
e13=open(r'e13.tsv','w')
e14=open(r'e14.tsv','w')
h11=open(r'h11.tsv','w')
h12=open(r'h12.tsv','w')
h13=open(r'h13.tsv','w')
h14=open(r'h14.tsv','w')

for n in range(0,len(w)):
    Gw = I.cluster_Green_function(w[n]+0.1j, clus=0)
    G11 = Gw[0,0]
    G12 = Gw[0,1]
    G13 = Gw[0,2]
    G14 = Gw[0,3]
    f11.write(str(w[n])+'        '+str(np.real(np.array([G11]))[0])+'        '+str(np.imag(np.array([G11]))[0])+'\n')
    f12.write(str(w[n])+'        '+str(np.real(np.array([G12]))[0])+'        '+str(np.imag(np.array([G12]))[0])+'\n')
    f13.write(str(w[n])+'        '+str(np.real(np.array([G13]))[0])+'        '+str(np.imag(np.array([G13]))[0])+'\n')
    f14.write(str(w[n])+'        '+str(np.real(np.array([G14]))[0])+'        '+str(np.imag(np.array([G14]))[0])+'\n')
    Ew = I.cluster_self_energy(w[n]+0.1j, clus=0)
    E11 = Ew[0,0]
    E12 = Ew[0,1]
    E13 = Ew[0,2]
    E14 = Ew[0,3]
    e11.write(str(w[n])+'        '+str(np.real(np.array([E11]))[0])+'        '+str(np.imag(np.array([E11]))[0])+'\n')
    e12.write(str(w[n])+'        '+str(np.real(np.array([E12]))[0])+'        '+str(np.imag(np.array([E12]))[0])+'\n')
    e13.write(str(w[n])+'        '+str(np.real(np.array([E13]))[0])+'        '+str(np.imag(np.array([E13]))[0])+'\n')
    e14.write(str(w[n])+'        '+str(np.real(np.array([E14]))[0])+'        '+str(np.imag(np.array([E14]))[0])+'\n')
    Hw = I.hybridization_function(w[n]+0.1j, clus=0)
    H11 = Hw[0,0]
    H12 = Hw[0,1]
    H13 = Hw[0,2]
    H14 = Hw[0,3]
    h11.write(str(w[n])+'        '+str(np.real(np.array([H11]))[0])+'        '+str(np.imag(np.array([H11]))[0])+'\n')
    h12.write(str(w[n])+'        '+str(np.real(np.array([H12]))[0])+'        '+str(np.imag(np.array([H12]))[0])+'\n')
    h13.write(str(w[n])+'        '+str(np.real(np.array([H13]))[0])+'        '+str(np.imag(np.array([H13]))[0])+'\n')
    h14.write(str(w[n])+'        '+str(np.real(np.array([H14]))[0])+'        '+str(np.imag(np.array([H14]))[0])+'\n')

f11.close()
f12.close()
f13.close()
f14.close()
e11.close()
e12.close()
e13.close()
e14.close()
h11.close()
h12.close()
h13.close()
h14.close()

# Obtain the matsubara frequency cluster Green's functions, self-energies and hybridizations

f11=open(r'g11_mats.tsv','w')
f12=open(r'g12_mats.tsv','w')
f13=open(r'g13_mats.tsv','w')
f14=open(r'g14_mats.tsv','w')
e11=open(r'e11_mats.tsv','w')
e12=open(r'e12_mats.tsv','w')
e13=open(r'e13_mats.tsv','w')
e14=open(r'e14_mats.tsv','w')
h11=open(r'h11_mats.tsv','w')
h12=open(r'h12_mats.tsv','w')
h13=open(r'h13_mats.tsv','w')
h14=open(r'h14_mats.tsv','w')

for n in range(0,len(wn)):
    Gw = I.cluster_Green_function(wn[n]*1j, clus=0)
    G11 = Gw[0,0]
    G12 = Gw[0,1]
    G13 = Gw[0,2]
    G14 = Gw[0,3]
    f11.write(str(wn[n])+'        '+str(np.real(np.array([G11]))[0])+'        '+str(np.imag(np.array([G11]))[0])+'\n')
    f12.write(str(wn[n])+'        '+str(np.real(np.array([G12]))[0])+'        '+str(np.imag(np.array([G12]))[0])+'\n')
    f13.write(str(wn[n])+'        '+str(np.real(np.array([G13]))[0])+'        '+str(np.imag(np.array([G13]))[0])+'\n')
    f14.write(str(wn[n])+'        '+str(np.real(np.array([G14]))[0])+'        '+str(np.imag(np.array([G14]))[0])+'\n')
    Ew = I.cluster_self_energy(wn[n]*1j, clus=0)
    E11 = Ew[0,0]
    E12 = Ew[0,1]
    E13 = Ew[0,2]
    E14 = Ew[0,3]
    e11.write(str(wn[n])+'        '+str(np.real(np.array([E11]))[0])+'        '+str(np.imag(np.array([E11]))[0])+'\n')
    e12.write(str(wn[n])+'        '+str(np.real(np.array([E12]))[0])+'        '+str(np.imag(np.array([E12]))[0])+'\n')
    e13.write(str(wn[n])+'        '+str(np.real(np.array([E13]))[0])+'        '+str(np.imag(np.array([E13]))[0])+'\n')
    e14.write(str(wn[n])+'        '+str(np.real(np.array([E14]))[0])+'        '+str(np.imag(np.array([E14]))[0])+'\n')
    Hw = I.hybridization_function(wn[n]*1j, clus=0)
    H11 = Hw[0,0]
    H12 = Hw[0,1]
    H13 = Hw[0,2]
    H14 = Hw[0,3]
    h11.write(str(wn[n])+'        '+str(np.real(np.array([H11]))[0])+'        '+str(np.imag(np.array([H11]))[0])+'\n')
    h12.write(str(wn[n])+'        '+str(np.real(np.array([H12]))[0])+'        '+str(np.imag(np.array([H12]))[0])+'\n')
    h13.write(str(wn[n])+'        '+str(np.real(np.array([H13]))[0])+'        '+str(np.imag(np.array([H13]))[0])+'\n')
    h14.write(str(wn[n])+'        '+str(np.real(np.array([H14]))[0])+'        '+str(np.imag(np.array([H14]))[0])+'\n')

f11.close()
f12.close()
f13.close()
f14.close()
e11.close()
e12.close()
e13.close()
e14.close()
h11.close()
h12.close()
h13.close()
h14.close()

# Generates the real frequency momentum-dependent periodized spectral function A(k,w) and imaginary part of the self-energy E(k,w) in the first BZ quadrant

f1=open(r'spectral.tsv','w')
f2=open(r'sigma.tsv','w')

for n in range(0,len(w)):
    for i in range(0,len(kx_red)):
        for j in range(0,len(ky_red)):
            G_per = I.periodized_Green_function(w[n]+0.1*1j, np.array([kx_red[i],ky_red[j],0.0]))
            E_per = I.self_energy(w[n]+0.1*1j, np.array([kx_red[i],ky_red[j],0.0]))
            spec =-1*np.imag(G_per)/np.pi
            sigma = np.imag(E_per)
            f1.write(str(round(w[n],4))+'        '+str(round(kx_red[i],6))+'        '+str(round(ky_red[j],6))+'        '+str(float(spec))+'\n')
            f2.write(str(round(w[n],4))+'        '+str(round(kx_red[i],6))+'        '+str(round(ky_red[j],6))+'        '+str(float(sigma))+'\n')

f1.close()
f2.close()

print('End of run')
