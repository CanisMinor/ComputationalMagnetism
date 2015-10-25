import numpy as np

def neigh_1D(n_x):
    N_s = n_x   #number of states
    N_n = 2     #each site has two neighbours in one dimension

    #create neighbour list
    neighbour = np.zeros(N_n, N_s)

    #create phase list
    #integers in phase list indicate phase:
    #     0-no shift, 1-shift left, 2-shift right
    phase = np.zeros(N_n, N_s)

    #find neighbours for each atom
    for iAtom in range(0, N_s):
        neighbour[0,iAtom] = iAtom - 1
        neighbour[1,iAtom] = iAtom + 1

    #fix periodic boundary conditions
    neighbour[0, 0] = n_x - 1
    phase[0, 0] = 1
    neighbour[1, n_x - 1] = 0
    phase[1, n_x - 1] = 2

    return neighbour, phase

