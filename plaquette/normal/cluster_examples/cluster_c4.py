import pyqcm
#================================================================================
### Defining the clusters ###

# The cluster definition don't have, 'a priori' any geometrical interpretation, but if one wants to implicity include some geometrical constrains, the definitions of the connections needs to be carefully done, taking into account the geometry of the system !!!

ns = 4							# Number of impurity sites
nb = 8							# Number of bath sites
no = ns + nb
cm = pyqcm.cluster_model(ns, nb, 'cluster')		# Class that defines the model for the impurity+bath 

# The cluster model indexes the impurity sites and the bath sites in a list with format [imp up, bath up, imp dw, bath dw], in our case is a 24 element list (12 up and 12 dw)

cm.varia = []						# Creates an empty list of bath parameters (varia <=> variables of the bath), we will add elements to it below
#================================================================================
### Bath related operators ###

# Here they are defined in a compact way, but you can list one-by-one deffining each new_operator explicitly

# Generates epsilon matrix (On-site energies)

for i in range(1,nb+1):			
    name = 'eb'+str(i)
    lab = i+ns
    cm.new_operator(name, 'one-body', [	    		# Defines a new operator connecting each site with itself
        (lab, lab, 1.0),                    		# Up
        (lab + no, lab + no, 1.0)           		# Dw
    ])
    cm.varia += [name]			    		# Add the names of the operators to the bath parameters

# Generates theta matrix (Impurity-bath hoppings)

def new_tb(x, seq):					# Function that for each 'x' (bath site) defines its connection with the impurity sites, the amplitudes come from a list called 'seq' 
    elem = []
    for i in [1,2,3,4]:                 # Im
        elem.append((i, x+ns, seq[i-1]))        # Up
        elem.append((i+no, x+ns+no, seq[i-1]))      # Dw

    cm.new_operator('tb'+str(x), 'one-body', elem)  # Defines a new operator connecting every cluster site to 'x'
    cm.varia += ['tb'+str(x)]       # Add the names of the operators to the bath parameters

new_tb(1, [1, 1, 1, 1])
new_tb(2, [1, 1,-1,-1])
new_tb(3, [1,-1,-1, 1])
new_tb(4, [1, 1, 1, 1])
new_tb(5, [1, 1,-1,-1])
new_tb(6, [1,-1,-1, 1])
new_tb(7, [1,-1, 1,-1])
new_tb(8, [1,-1, 1,-1])
#================================================================================
