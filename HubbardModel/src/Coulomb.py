#-----------------------------------------|
#   	       COULOMB TERM		  |
#-----------------------------------------|

import numpy as np

def coulomb(C, N_s):
    H_c = C * np.ones(N_s, N_s)

    return H_c
