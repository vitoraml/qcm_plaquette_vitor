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

def new_tb(x, seq):                 # Function that for each 'x' (bath site) defines its connection with the impurity sites, the amplitudes come from a list called 'seq'
    eleml = []
    for i in [1,3]:
        eleml.append((i, x+ns, seq[i-1]))         # Up
        eleml.append((i+no, x+ns+no, seq[i-1]))   # Dw
    elemr = []
    for i in [2,4]:
        elemr.append((i, x+ns, seq[i-1]))         # Up
        elemr.append((i+no, x+ns+no, seq[i-1]))   # DW

    cm.new_operator('tbl'+str(x), 'one-body', eleml) # Defines a new operator connecting cluster sites 1 and 3 (left sites) to 'x'
    cm.new_operator('tbr'+str(x), 'one-body', elemr) # Defines a new operator connecting cluster sites 2 and 4 (right sites) to 'x'
    cm.varia += ['tbl'+str(x), 'tbr'+str(x)] # Add the names of the operators to the bath parameters

for i in range(1,5):
    new_tb(i, [1, 1, 1, 1])
for i in range(5,9):
    new_tb(i, [1, 1, -1, -1])

# Generates delta matrix (Anomalous impurity-bath hoppings)

def new_pairing(x, seq):                 # Function that for each 'x' (bath site) with spin y, defines its connection with the spin -y impurity sites, the amplitudes come from a list called 'seq'
    eleml = []
    for i in [1,3]:
        eleml.append((i, x+ns+no, seq[i-1]))      # Up
        eleml.append((x+ns, i+no, seq[i-1]))      # Dw
    elemr = []
    for i in [2,4]:
        elemr.append((i, x+ns+no, seq[i-1]))      # Up
        elemr.append((x+ns, i+no, seq[i-1]))      # Dw
    
    cm.new_operator('sbl'+str(x), 'anomalous', eleml) # Defines a new operator connecting cluster sites 1 and 3 (left sites) spin y to spin -y 'x'
    cm.new_operator('sbr'+str(x), 'anomalous', elemr) # Defines a new operator connecting cluster sites 2 and 4 (right sites) spin y to spin -y 'x'
    cm.varia += ['sbl'+str(x), 'sbr'+str(x)] # Add the names of the operators to the bath parameters

for i in range(1,5):
    new_pairing(i, [1, 1, 1, 1])
for i in range(5,9):
    new_pairing(i, [1, 1, -1, -1])
#================================================================================
