#subroutines for operations on the exchange data

import numpy as np

#count the number of exchange pairs listed in data file
def CountExch(exchFile):
    print("counting no of interactions in file " + exchFile + '\n')
    nCount = 0
    with open(exchFile, 'r') as eF:
        for line in eF:
            nCount = nCount + 1
    return nCount



#read in exchange tensors
def ReadEffExch(exchFile, n_interactions):
    Jij_array = np.zeros((n_interactions, 9))
    inform = np.zeros((n_interactions, 5))
    interat = np.zeros(n_interactions)
    j = 0
    with open(exchFile, 'r') as exchT:
        for line in exchT:
            qLine = line.split()
            inform[j,0:5] = np.array(map(float, qLine[0:5]))
            interat[j] = float(qLine[5])
            Jij_array[j, 0:9] = np.array(map(float, qLine[6:15]))
            j = j + 1
    return Jij_array, inform, interat
