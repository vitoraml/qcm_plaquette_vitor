import pyqcm
from cluster import cm 									# Imports our cluster model class defined in cluster.py 
#================================================================================
### Defining the lattice model ###

# Here we give a geometry to the impurity

pos_cu = [[0,0,0], [1,0,0], [0,1,0], [1,1,0]] 						# Impurity sites positions
clus_cu = pyqcm.cluster(cm, pos_cu)							# Attach the geometry to our previously defined cluster model (cluster = cluster model + geometry)
model = pyqcm.lattice_model('model', clus_cu, [[2,0,0],[0,2,0]],[[1,0,0],[0,1,0]])	# Defines the lattice model , the 3rd entry is a list of lattice generators (vectors connecting different u.c.) and the 4th entry is a list of unit cell generators (vectors connecting different impurity sites)

#================================================================================
### Lattice operators ###

# The operators are defined restrained to the impurity sites

# U (on-site Coulomb repulsion)
# Compact definition with implicit flags
model.interaction_operator('U')

# t (1st neighbor hopping)
# Format [name, vector connecting the sites, amplitude]
model.hopping_operator('t', [1,0,0], -1)	# Hopping in x
model.hopping_operator('t', [0,1,0], -1)	# Hopping in y

# t' (2nd neighbor hopping)
model.hopping_operator('t1', [1,1,0], 1)	# Hopping in x=y
model.hopping_operator('t1', [1,-1,0], 1)	# Hopping in x=-y

# Superconductivity (dx2-y2)
model.anomalous_operator('D', [1,0,0], 1, orbitals=(1,1))     # Pairing in x
model.anomalous_operator('D', [0,-1,0], -1, orbitals=(1,1))   # Pairing in y

#================================================================================
