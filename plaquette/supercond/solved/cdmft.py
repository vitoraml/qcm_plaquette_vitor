import numpy as np
import pandas as pd
import os
import pyqcm
from model import model			# Imports the lattice_model defined in model.py
from pyqcm.cdmft import CDMFT		# Imports the ED-CDMFT routine from pyqcm
from scipy.optimize import brentq
#================================================================================
### Parameters and variables ###

# Here we provide the setup and input parameters for the CDMFT run

# Pyqcm general flags
pyqcm.set_global_parameter('max_iter_lanczos', 4000)		# Max Lanczos iterations per CDMFT iteration (def=600)
pyqcm.set_global_parameter('Ground_state_method','P')       # Set ED algorithm to PRIMME (def='L' -> Lanczos)

# Cluster model parameters
var=[]								# Defines an auxiliar list called var
for j in range(1,2):						# Appends to var the names of the bath parameters
    for i in range(1,9):
        var.append('eb{:d}_{:d}'.format(i,j))
        var.append('tbl{:d}_{:d}'.format(i,j))
        var.append('tbr{:d}_{:d}'.format(i,j))
        var.append('sbl{:d}_{:d}'.format(i,j))
        var.append('sbr{:d}_{:d}'.format(i,j))

# Sets the bath parametrization values
varia_values = """
eb1_1  =  -0.36052554119999997
eb2_1  =  -0.36714693719999997
eb3_1  =  -0.28882774780000003
eb4_1  =  -1.104689469
eb5_1  =  0.28895992239999996
eb6_1  =  1.057548416
eb7_1  =  0.3338090609
eb8_1  =  1.104773726
tbl1_1  =  0.008822598792
tbl2_1  =  0.22523938920000003
tbl3_1  =  -0.20672429899999997
tbl4_1  =  -0.08260495222
tbl5_1  =  0.16155771
tbl6_1  =  -0.5595635917
tbl7_1  =  -9.545390293000001e-05
tbl8_1  =  0.4135437228
tbr1_1  =  0.008933873349000001
tbr2_1  =  0.2253640178
tbr3_1  =  0.2066009053
tbr4_1  =  0.08259035436
tbr5_1  =  0.1616020358
tbr6_1  =  0.5586754699
tbr7_1  =  -9.346360571e-06
tbr8_1  =  0.4147647583
sbl1_1  =  0.2446412154
sbl2_1  =  -0.009564807891000001
sbl3_1  =  0.1614036478
sbl4_1  =  -0.4141815457
sbl5_1  =  -0.2067227257
sbl6_1  =  -7.005617851999999e-05
sbl7_1  =  0.2490237356
sbl8_1  =  0.08265706661
sbr1_1  =  0.24449300300000001
sbr2_1  =  -0.009637863231999999
sbr3_1  =  -0.16154167460000002
sbr4_1  =  0.4142173339
sbr5_1  =  -0.20682568910000002
sbr6_1  =  -9.04441114e-05
sbr7_1  =  -0.24894709510000002
sbr8_1  =  0.08260868796
"""

# Lattice model parameters (chemical potential 'mu' is mandatory). Since now we're dealing with a possible superconducting (sc) solution, we defined the superconducting order parameter (D) in the model and here we atribute him a very small value to not affect the hamiltonian. This way, if a sc state is the ground state of the system, the expected value of D will be our order parameter. 
#In some sense we're just 'allowing' the system to be superconducting, but not forcing it to. If the provided D in the lattice model is too big (in some of my tests more than 1e-3), this can lead to erroneous results and an overestimation of the superconducting order if compared to the literature.
band_params = """
U = 8
t = 1.0
t1 = 0.3
mu = 1.5
D = 1e-6
"""

# Lattice model configuration
# Since we're allowing creation/destruction of pairs, there's no 'constant particle sector' anymore. Our proposed allowed symmetry for the sc order is singlet or S0 triplet, so we, in principle, don't need to look in other spin sectors.
model.set_target_sectors('R0:S0')	# Sector to search for 
model.set_parameters(band_params + varia_values)		# Insert the bath and lattice parameters defined before in the model
#================================================================================
### CDMFT run ###

target = 0.8631967215
bracket = [1.4, 1.6]

# Main CDMFT
# When including sc we need also to converge the superconducting order parameter in the self-consistent cycle
def task(mu):
    model.set_parameter('mu', mu)
    CDMFT(model, varia=var, accur=(1*1e-6,1*1e-8), convergence=('self-energy','D'), accur_bath=1*1e-10, depth=2, accur_dist=1e-10, maxiter = 137000, max_value=10000, max_function_eval=5000000000, alpha=0.0, initial_step=0.05, beta=100, wc=2, grid_type='sharp')
    I = pyqcm.model_instance(model)
    n = I.averages()['mu']
    print('mu = {:g}, n = {:g}'.format(I.parameters()['mu'], n))
    return n-target

x0, r = brentq(task, bracket[0], bracket[1], xtol=1e-4, maxiter=100, full_output=True, disp=True)

print(r.flag)
if not r.converged:
    raise RuntimeError('the root finding routine could not find a solution!')

# CDMFT format:

# model 			  -> lattice model
# varia 			  -> list of bath free parameters (def=explicit list)
# accur				  -> convergence tolerance (def=1*1e-4)
# convergence         -> parameters that will be converged in the cycle
# accur_bath		  -> tolerance distance function minimization(def=1*1e-3)
# depth				  -> number of previous iterations that convergence will be checked (def=2)
# accur_dist		  -> relative tolerance distance function (def=1*1e-5)
# maxiter			  -> CDMFT max iteration (def=32)
# max_value			  -> max allowed value for the free parameters before crash (def=100)
# max_function_eval	  -> max number of distance function evaluations (def=500000)
# alpha               -> precentual 'weight' of the previous iteration in the new one (def=0.0)
# initial_step		  -> percentual 'size' in parameter space of the initial iteration step (def=0.1)
# beta                -> inverse fictitious temperature for matsubara (def=50)
# wc                  -> matsubara frequency cutoff for the weight function (def=2)
# grid_type           -> weight choice for each matsubara frequency inside the wc window (def='sharp')

# The converged solution is written in cdmft.tsv file and every iteration solution is written in cdmft_iter.tsv

#The options above are the main ones, for more information go to the documentation
#================================================================================
### Basic outputs ###

# The instance is our best friend, but she's lazy (really lazy). Since most quantities are calculated 'on demmand' to save cpu time, RAM and disk space, with the model information stored and updated in our instance, we need to ask her exactly what we need as output.

# Some outputs are more straightfoward, while others need a more 'detailed' request, here we will focus on the first type ones, but in pp_cdmft there are examples on how to extract information in the second type ones. 

I = pyqcm.model_instance(model)									# Generates an instance from our lattice model
I.Green_function_solve()									    # Calculate and store the Green's function explicitly (non-lazy)
I.averages(pr=True)										        # Prints the averages of the lattice model operators
I.plot_DoS(w=4, eta=0.1, sum=False, progress=True, colors=None, file='dos_converged.pdf')			    # Plot the DOS in a real frequency grid, also writes a dos.tsv data file
I.mdc(nk=200, eta=0.1, opt='GF', k_perp=0, freq=0.0, plane='xy', size=1.0, file='mdc_converged.pdf')	# Plot a FS BZ planar cut in a given real frequency, other functions can be plotted
I.mdc_anomalous(nk=200, w=0.1j, opt='GF', orbitals=(1,1), k_perp=0, plane='xy', file='mdc_anomalous_converged.pdf')    # Plot the anomalous part of the spectral function in a BZ plane at a given complex frequency 'w'

# Rename the main output files for a more practical data treatment
os.system('mv cdmft.tsv cdmft_converged.tsv')
os.system('mv dos.tsv dos_converged.tsv')

#The options above are just examples, for more information go to the documentation
#================================================================================
