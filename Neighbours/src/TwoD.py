import numpy as np

def neigh_2D(n_x, n_y):
    N_s = n_x * n_y   #number of states
    N_n = 4           #each site has four neighbours in two dimensions

    #create neighbour list
    neighbour = np.zeros((N_s, N_n))

    #create phase list
    #integers in phase list indicate phase:
    #     0-no shift, 1-shift left, 2-shift right, 3-shift up, 4-shift down
    phase = np.zeros((N_s, N_n))

    #find neighbours for each atom
    for iAtom in range(0, N_s):
        neighbour[iAtom,0] = iAtom - 1
        neighbour[iAtom,1] = iAtom + 1
        neighbour[iAtom,2] = iAtom - n_x
        neighbour[iAtom,3] = iAtom + n_x

    #fix periodic boundary conditions for left side
    for iAtom in range(0, N_s, n_x):
        neighbour[iAtom,0] = iAtom + (n_x - 1)
        phase[iAtom,0] = 1

    #fix periodic boundary conditions for right side
    for iAtom in range(n_x - 1, N_s, n_x):
        neighbour[iAtom,1] = iAtom - (n_x - 1)
        phase[iAtom,1] = 1

    #fix periodic boundary conditions for top
    for iAtom in range(0, n_x):
        neighbour[iAtom,2] = iAtom + (n_x * (n_y - 1))
        phase[iAtom,2] = 1

    #fix periodic boundary conditions for bottom
    for iAtom in range(n_x * (n_y - 1), N_s):
        neighbour[iAtom,3] = iAtom - (n_x * (n_y - 1))
        phase[iAtom,3] = 1

    return neighbour, phase





