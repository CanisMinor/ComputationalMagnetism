#neighbour-related functions
import numpy as np

def ReadNeighbours(neighbfile, list_length):
    print('reading in neighbours from neighbour list ' + neighbfile + '\n')

    site_i = np.zeros(list_length)
    site_j = np.zeros(list_length)
    Jij_no = np.zeros(list_length)
    
    with open(neighbfile, 'r') as nF:
        j = 0
        for line in nF:
            k = 3
            qLine = line.split()
            site_i[j] = qLine[0]
            site_j[j] = qLine[1]
            #junk in third column
            Jij_no[j] = qLine[3]
            j = j + 1


    return site_i, site_j, Jij_no




