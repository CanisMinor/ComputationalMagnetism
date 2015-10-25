#Implementation of the Two-Dimensional Hubbard Model
import numpy as np
import math
import sys
import cmath

sys.path.append("../../Neighbours/src/")

import TwoD as neigh
import OnsiteRepulsion as onsite
import BandEnergy as band
import Coulomb as cou
import Hopping as hop


#iteration parameters
N_iter = 100
a_x = 1.0
a_y = 1.0
U_max = 0.0
t = 1.0
n_x = 4
n_y = 4


#initialise problem
N_s = n_x * n_y          #total no of sites in unit cell
k_latx = -math.pi           #initialising wave vector along x
k_laty = -math.pi           #initialising wave vector along y
n_elec =  N_s          #no of electrons (system is half-filled)
delta_U = 0.1        #Hubbard U increment size
n_do = n_elec/2.0
n_up = n_elec/2.0
n_do = 0# n_elec/2        #no of down electrons
n_up = 4 #n_elec/2        #no of up electrons
# IF(MOD(n_elec,2) .NE. 0) THEN
#    n_up = n_up +1
# END IF

print('No of sites: '+ str(N_s))
print('No of down electrons: ' + str(n_do))
print('No of up electrons: ' + str(n_up))

#initialise algorithm
tol = 0.000001   #tolerance for convergence
nkx = 100    #no of iterations over k_x
nky = nkx   #no of iterations over k_y
delta_k = (2.0 * math.pi)/float(nkx)

#initialising the k-space weighting factor# s

weights = np.ones((nky + 1, nkx + 1))
weights[0, :] = 0.5 * np.ones(nkx + 1)    #sides; share between two unit cells
weights[nky, :] = 0.5 * np.ones(nkx + 1)
weights[:, 0] = 0.5 * np.ones(nky + 1)
weights[:, nkx] = 0.5 * np.ones(nky + 1)
weights[0,0] = 0.25         #corners; share between four unit cells
weights[0, nkx] = 0.25
weights[nky, 0] = 0.25
weights[nky, nkx] = 0.25
kint = math.pi/(a_x *float(nkx))
# weights(:,:) = weights(:,:)*4*kint*kint/(4*pi*pi)

print('nkx = ' + str(nkx))
print('nky = ' + str(nky))

neighbour, phase = neigh.neigh_2D(n_x, n_y)



# allocate space to arrays
nup_tot = np.zeros(N_s)
ndown_tot = np.zeros(N_s)

#initialise Hubbard U array
U = U_max * np.ones(N_s)


#-------------------------------------------|
#              BAND OCCUPANCY               |
#-------------------------------------------|

for s in (0,N_iter):    #iterate until convergence
    # #down-spin scan through k-space
    ndown_tot = np.zeros(N_s)
    k_latx = - math.pi
    for l in range(0,nkx+1):
        k_laty = - math.pi     #=k_latx  #start at k_latx to get a wedge in BZ
        for r in range(0,nky+1):   #iterate over different wave vectors
            #DO i = 1,N_s
            #   U[i] = rand(itime)            #initialise Hubbard U
            #END DO

            #calculate Hamiltonian parts for up Hamiltonian
            H_t = hop.hopping(t, N_s,a_x,a_y,k_latx,k_laty, neighbour, phase)
            H_o = onsite.hubU(U, N_s, ndown_tot)

            #evaluate total down Hamiltonian
            H_down = np.zeros((N_s, N_s), dtype=complex)
            for i in range(0,N_s):
                for j in range(0,N_s):
                    H_down[j,i] = H_t[j,i] + complex(H_o[j,i])

            #calculate average spin-down density
            n_av_down = np.zeros(N_s)
            for j in range(0,n_do):
                for i in range(0,N_s):
                    n_av_down[i] = n_av_down[i] + abs(H_down[i,j]) * abs(H_down[i,j])

            #add to total spin-down density
            for i in range(0,N_s):
                ndown_tot[i] = ndown_tot[i] + (n_av_down[i] * weights[l,r])#/(DBLE((nkx-1)*(nky-1)))#*weights(l,r)

            k_laty = k_laty + delta_k
        k_latx = k_latx +delta_k

    for q in range(0,N_s):
        ndown_tot[q] = ndown_tot[q]/(float(nkx*nky))

    nup_tot = np.zeros(N_s)
    k_latx = -math.pi
    for l in range(0,nkx+1):
        k_laty = -math.pi
        for r in range(0, nky+1):
            #calculate Hamiltonian parts for down Hamiltonian
            H_t = hop.hopping(t,N_s,a_x,a_y, k_latx, k_laty, neighbour, phase)
            H_o = onsite.hubU(U, N_s, nup_tot)

            #evaluate total up Hamiltonian
            H_up = np.zeros((N_s, N_s), dtype=complex)
            for i in range(0,N_s):
                for j in range(0,N_s):
                    H_up[j,i] = H_t[j,i] + complex(H_o[j,i])

            #calculate average spin-up density
            n_av_up = np.zeros(N_s)
            for j in range(0, n_up):
                for i in range(0, N_s):
                    n_av_up[i] = n_av_up[i] + (abs(H_up[i,j]) * abs(H_up[i,j]))

            for i in range(0, N_s):
                nup_tot[i] = nup_tot[i] + (n_av_up[i] * weights[l,r])

            #exit loop if up and down densities have converged
            #if (s > 2) and max(change_up) < tol and max(change_down) < tol:
            #    print('exiting at iteration ' + str(r) + '\n')
            #    break;

            k_laty = k_laty +delta_k
        k_latx = k_latx + delta_k

    for q in range(0, N_s):
        nup_tot[q] = nup_tot[q]/(float(nkx*nky))


#write converged densities to file
#(loop scans site index)
with open("../data/occupancy.dat", 'w+') as outF:
    for i in range(0, N_s):
        outF.write(str(i) + ", " + str(nup_tot[i]) + ", " + str(-ndown_tot[i]) + '\n')



#-------------------------------------------|
#     BAND STRUCTURE & TOTAL ENERGY         |
#-------------------------------------------|

k_latx = - math.pi
H_o = np.zeros((N_s, N_s))
totE = np.zeros(N_s)
with open("../data/bands.dat", 'w+') as outF, open("../data/contourbands.dat", 'w+') as contourF:
    for r in range(0, nkx+1):
        k_laty = -math.pi
        for l in range(0, nky+1):
            #calculate Hamiltonian parts for down Hamiltonian
            H_t = hop.hopping(t, N_s,a_x,a_y,k_latx,k_laty, neighbour, phase)
            H_o = onsite.hubU(U,N_s,ndown_tot)

            #evaluate total up Hamiltonian
            for i in range(0, N_s):
                for j in range(0, N_s):
                    H_up[j,i] = H_t[j,i] + H_o[j,i]

            #find eigenvalues and eigenfunctions of up Hamiltonian
            w_up, v_up = np.linalg.eig(H_up)

            #calculate Hamiltonian parts for up Hamiltonian
            H_t = hop.hopping(t,N_s,a_x,a_y,k_latx,k_laty, neighbour, phase)
            H_o = onsite.hubU(U,N_s,nup_tot)

            #evaluate total down Hamiltonian
            for i in range(0,N_s):
                for j in range(0,N_s):
                    H_down[j,i] = H_t[j,i] + H_o[j,i]

            #find eigenvalues and eigenfunctions of down Hamiltonian
            w_down, v_down = np.linalg.eig(H_down)

            for i in range(0, N_s):
                outF.write(str(k_latx) + ", " + str(k_laty) + ", " + str(abs(w_up[i])) + ", " + str(abs(w_down[i])) + '\n')

            #total energy for each site from eigenenergies
            energy_onsite = np.zeros(N_s)
            for i in range(0, N_s):
                for j in range(0, n_up): #abs of eigenvalue will be approximately equal to the real part as img part are almost zero
                    energy_onsite[i] = energy_onsite[i] + abs(H_up[i,j]) * abs(H_up[i,j]) * abs(w_up[j])
                for j in range(0, n_do):  #abs of eigenvalue will be approximately equal to the real part as img part are almost zero
                    energy_onsite[i] = energy_onsite[i] + abs(H_down[i,j]) * abs(H_down[i,j]) * abs(w_down[j])

            #add double-occupancy energies
            for i in range(0, N_s):
                energy_onsite[i] = energy_onsite[i] + U[i]*ndown_tot[i]*nup_tot[i]

            #integrate total energy
            for i in range(0, N_s):
                totE[i] = totE[i] + (energy_onsite[i]*weights[r,l])
                contourF.write(str(abs(w_up[i])) + ', ')

            contourF.write('\n')  #for contour plot

            k_laty = k_laty + delta_k  #increment k_laty

        #WRITE(33,*) ' '  #for contour plot
        k_latx = k_latx + delta_k  #increment k_latx


#normalise total energy
for i in range(0, N_s):
    totE[i] = totE[i]/(float(nkx*nky))


#write total energy to file
with open("../data/totE.dat", 'w+') as eF:
    for i in range(0, N_s):
        eF.write(str(totE[i]) + '\n')





