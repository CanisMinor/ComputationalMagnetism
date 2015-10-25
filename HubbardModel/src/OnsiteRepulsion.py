#-----------------------------------------|
#   	       HUBBARD TERM		  |
#-----------------------------------------|
import numpy as np

def hubU(U, N_s, n_av):
    H_o = np.zeros((N_s, N_s))

    for i in range(0,N_s):
        H_o[i,i] = U[i] * n_av[i]

    return H_o


