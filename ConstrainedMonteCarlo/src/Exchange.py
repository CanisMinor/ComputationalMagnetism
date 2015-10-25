#subroutines for operations on the exchange data

import numpy as np

#count the number of exchange pairs listed in data file
def CountExch(exchFile):
    print("counting no of interactions in file " + s + '\n')
    nCount = 0
    with(exchFile, r) as eF:
        for line in eF:
            nCount = nCount + 1

    return nCount



#read in exchange tensors
def ReadEffExch(exchFile, n_interactions):
    Jij_array = np.zeros(n_interactions, 9)
    inform = np.zeros(n_interactions, 5)
    interat = np.zeros(n_interactions)
    with open(exchFile, 'r') as exchT:
        for line in exchT:
            q = 3 #SPLIT LINE HERE
            #   exch_templ >> inform[j][0] >> inform[j][1] >> inform[j][2] >> inform[j][3] >> inform[j][4] >> interat[j] >> Jij_array[j][0] >> Jij_array[j][1] >> Jij_array[j][2] >> Jij_array[j][3] >> Jij_array[j][4] >> Jij_array[j][5] >> Jij_array[j][6] >> Jij_array[j][7] >> Jij_array[j][8];
    return Jij_array, inform, interat




 #if(!exch_templ) print("effective exchange tensors could not be read in" +endl;

 #read in exchange parameters and pairs info from KKR data
 #for(int j=0; j<n_interactions; j++)
 #{
 #}


 #exch_templ.close();

