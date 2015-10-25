import numpy as np
import sys

sys.path.append("../../Neighbours/src/")

def neigh_3D(n_x, n_y, n_z):

    '''
    This code labels the atoms in a 3D unit cell in increasing
    integers starting along the x-axis, secondly along the y-axis
    and finally along the z-axis in a right-handed coordinate system.
    The labelling thus starts at the origin (considered close top
    left-hand corner) with atom 0 and finishes at (n_x-1,n_y-1,n_z-1)
    (the far bottom right-hand corner) with atom N_s-1.
    '''

    N_s = n_x * n_y * n_z   #number of states
    N_n = 6                 #each site has six neighbours in three dimensions

    #create neighbour list
    neighbour = np.zeros(N_n, N_s)

    #create phase list
    #integers in phase list indicate phase:
    #     0-no shift, 1-shift left, 2-shift right, 3-shift up, 4-shift down, 5-shift forward, 6-shift backward
    phase = np.zeros(N_n, N_s)

    #find neighbours for each atom
    for iAtom in range(0, N_s):
        neighbour[0,iAtom] = iAtom - 1
        neighbour[1,iAtom] = iAtom + 1
        neighbour[2,iAtom] = iAtom - n_x
        neighbour[3,iAtom] = iAtom + n_x
        neighbour[4,iAtom] = iAtom - (n_x * n_y)
        neighbour[5,iAtom] = iAtom + (n_x * n_y)

    #fix periodic boundary conditions for left face
    for iAtom in range(0, N_s, n_x):
        neighbour[0, iAtom] = iAtom + (n_x - 1)
        phase[0,iAtom] = 1

    #fix periodic boundary conditions for right face
    for iAtom in range(n_x - 1, N_s, n_x):
        neighbour[1, iAtom] = iAtom - (n_x - 1)
        phase[1,iAtom] = 1


    #fix periodic boundary conditions for top face
    iCount = 0
    for iAtom in range(0, n_x * n_z):  #there are n_x * n_z atoms on top face
        neighbour[2, iCount] = iAtom #+ (n_x * (n_y -1)) #????? CHECK
        phase[2, iCount] = 1
        if (iAtom + 1) % n_x == 0:
            iCount = iCount + (n_x * n_y)
        else:
            iCount = iCount + 1

    #fix periodic boundary conditions for bottom face
    iCount = n_x*(n_y-1)+1  #?????? check
    for iAtom in range(0, n_x * n_z):  #there are n_x * n_z atoms on bottom face
        neighbour[3, iAtom] = iAtom  # -???
        phase[3, iAtom] = 1
        if (iAtom + 1) % n_x == 0:
            iCount = iCount + (n_x * n_y)
        else:
            iCount = iCount + 1


    #fix periodic boundary conditions for front face
    for iAtom in range(0, n_x * n_y):
        neighbour[4, iAtom] = iAtom #+ ????
        phase[4,iAtom] = 2

    #fix periodic boundary conditions for bottom face
    for iAtom in range(n_x * (n_y - 1), N_s):
        neighbour[5, iAtom] = iAtom #- ????
        phase[5,iAtom] = 3

    return neighbour, phase