#neighbour-related functions
import numpy as np

def ReadNeighbours(neighbfile, list_length):
    print('reading in neighbours from neighbour list ' + neighbfile + '\n')

    site_i = np.zeros(list_length)
    site_j = np.zeros(list_length)
    Jij_no = np.zeros(list_length)
    junk_stuff = np.zeros(list_length)
    
    with open(neighbfile, 'r') as nF:
        for line in nF:
            k = 3
            #read in step_i,step_j,neighbour no,tensor no
            #neighbfile >> site_i[n_line] >> site_j[n_line] >> junk_stuff[n_line] >> Jij_no[n_line];

    return site_i, site_j, Jij_no




