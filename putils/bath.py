"""
Project: Read bath parameters from converged CDMFT calculation

Description:

Data filetring of the bath parameters from a converged CDMFT output file.

Authors:
- Vitor Assunção Moreira Lima (vitor.aml@hotmail.com)
-

Created: 2025-05-05
Last Modified: 2025-06-01
Version: 1.0.0

License: MIT

Dependencies:
- Python3
- pandas
- numpy
"""

import pandas as pd
import numpy as np
#================================================================================

# Indicates the converged output file and use pandas to read it
parf= 'cdmft_converged.tsv'

bp = pd.read_csv(parf, sep='\t')

# Filter the desired quantities, in this case all the tb and eb values
eb = bp.filter(like='eb').iloc[:, 0:16]
tb = bp.filter(like='tb').iloc[:, 0:16]

# Since the averages are also included in the filter, here we just remove them
eb = eb[eb.columns[::2]]
tb = tb[tb.columns[::2]]

# Gathers both variables filtered
params = pd.concat([eb,tb],axis = 1)

# Obtain the variables names and reads the values associated with the last line (usefull if more than one bath output is included in the same file)
names =np.array(params.columns,dtype = object)
values = np.array(params)[-1]

# Loop in the variables and print them in terminal
i = -1
for n in names:
    i +=1
    print (n," = ",values[i])
