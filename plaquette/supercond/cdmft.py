import numpy as np
import pandas as pd
import os
import pyqcm
from model import model			# Imports the lattice_model defined in model.py
from pyqcm.cdmft import CDMFT		# Imports the ED-CDMFT routine from pyqcm
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
eb1_1  =  -0.2593375538
eb2_1  =  -0.2541980757
eb3_1  =  -0.2541980657
eb4_1  =  -0.2113528555
eb5_1  =  0.07336270796
eb6_1  =  0.7736925191
eb7_1  =  0.7736925097
eb8_1  =  1.01528381
tb1_1  =  0.2067308646
tb2_1  =  0.2198793405
tb3_1  =  0.2198793378
tb4_1  =  -0.2283733783
tb5_1  =  0.1680405159
tb6_1  =  -0.4067758811
tb7_1  =  -0.4067758835
tb8_1  =  0.5325887048
sbl1_1  =  0.3160640239
sbl2_1  =  4.300798014e-06
sbl3_1  =  2.86634553e-06
sbl4_1  =  -0.379465184
sbl5_1  =  -0.002426391765
sbl6_1  =  -3.390293888e-06
sbl7_1  =  0.3892796949
sbl8_1  =  -0.001956724715
sbr1_1  =  -0.3160641859
sbr2_1  =  -1.033980344e-07
sbr3_1  =  2.876551019e-06
sbr4_1  =  0.3794651519
sbr5_1  =  -0.002426248721
sbr6_1  =  3.599971297e-06
sbr7_1  =  -0.3892776593
sbr8_1  =  -0.001956563353
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

# Main CDMFT
# When including sc we need also to converge the superconducting order parameter in the self-consistent cycle
CDMFT(model, varia=var, accur=(1*1e-6,1*1e-6), convergence=('self-energy','D'), accur_bath=1*1e-10, depth=2, accur_dist=1e-10, maxiter = 137000, max_value=10000, max_function_eval=5000000000, alpha=0.0, initial_step=0.01, beta=100, wc=4, grid_type='sharp')

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
I.mdc(nk=200, w=0.1j, opt='GF', orbitals=(1,1), k_perp=0, plane='xy', file='mdc_anomalous_converged.pdf')    # Plot the anomalous part of the spectral function in a BZ plane at a given complex frequency 'w'

# Rename the main output files for a more practical data treatment
os.system('mv cdmft.tsv cdmft_converged.tsv')
os.system('mv dos.tsv dos_converged.tsv')

#The options above are just examples, for more information go to the documentation
#================================================================================
