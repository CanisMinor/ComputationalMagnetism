#-----------------------------------------|
#   	     BAND ENERGY TERM		      |
#-----------------------------------------|
import numpy as np

def elec_energy(epsi, N_s):
    return epsi * np.ones(N_s, N_s)

