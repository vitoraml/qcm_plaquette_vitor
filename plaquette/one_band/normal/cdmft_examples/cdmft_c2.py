import numpy as np
import pandas as pd
import os
import pyqcm
from model import model			# Imports the lattice_model defined in model.py
from pyqcm.cdmft import CDMFT		# Imports the ED-CDMFT routine from pyqcm
#================================================================================
### Parameters and variables ###

# Here we provide the setup and input parameters for the CDMFT run, it should converge in around 17 iterations with the provided seed

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

# Sets the bath parametrization values
varia_values = """
eb1_1  =  -0.2579227295
eb2_1  =  0.7743075062
eb3_1  =  0.07243100443
eb4_1  =  -0.2533705966
eb5_1  =  0.7743051498
eb6_1  =  -0.2111158235
eb7_1  =  -0.2533716499
eb8_1  =  1.015886176
tbl1_1  =  0.2065410313
tbl2_1  =  -0.4068388202
tbl3_1  =  0.1674822293
tbl4_1  =  -0.2197630016
tbl5_1  =  0.4068387291
tbl6_1  =  0.2283752831
tbl7_1  =  0.2197810452
tbl8_1  =  0.5326832017
tbr1_1  =  0.2065025962
tbr2_1  =  0.4068387735
tbr3_1  =  0.1674822779
tbr4_1  =  0.2197987499
tbr5_1  =  0.4068385021
tbr6_1  =  -0.2283752531
tbr7_1  =  0.2197811329
tbr8_1  =  -0.532683386
"""

# Lattice model parameters (chemical potential 'mu' is mandatory)
band_params = """
U = 8
t = 1.0
t1 = 0.3
mu = 1.5
"""

# Lattice model configuration
model.set_target_sectors('R0:N10:S0/R0:N12:S0/R0:N14:S0')	# List of sectors to search for the GS
model.set_parameters(band_params + varia_values)		# Insert the bath and lattice parameters defined before in the model
#================================================================================
### CDMFT run ###

# Main CDMFT
CDMFT(model, varia=var, accur=1*1e-6, convergence='self-energy', accur_bath=1*1e-10, depth=2, accur_dist=1e-10, maxiter = 137000, max_value=10000, max_function_eval=5000000000, alpha=0.01, initial_step=0.01, beta=100, wc=4, grid_type='sharp')

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

# Rename the main output files for a more practical data treatment
os.system('mv cdmft.tsv cdmft_converged.tsv')
os.system('mv dos.tsv dos_converged.tsv')

#The options above are just examples, for more information go to the documentation
#================================================================================
