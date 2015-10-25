#constrained Monte Carlo algorithm
#Canis Minor  Joanna.Hill@writeme.com
#13th Dec 2011

import numpy as np
import math
import random
import Parameters as pr
import Exchange as ex
import Neighbours as ne
import RotMatrix as rotate


#parameters for time- and temperature loops
N_iter = 2000         #no of time steps
N_eq = 1500           #equilibration iteration treshold

#constrained Monte Carlo variables
M_tot = 0.0  #magnetisation
M_x = 0.0    #components of magnetisation
M_y = 0.0
M_z = 0.0
summa=0.0

mz_unit= np.zeros(3)     #unit vector along z

nsp=0   #primary and secondary spin site indices
nsc=0

sp = np.zeros(3)
spp = np.zeros(3)
rsp = np.zeros(3)
rspp = np.zeros(3)
sc = np.zeros(3)
scp = np.zeros(3)
rsc = np.zeros(3)

norm=0.0            #normalisation temp variable

slump_rand = np.zeros(3)            #random number
sarg=0.0      #square root argument
deltaE=0.0    #change in energy
prand=0.0     #random number
pboltz=0.0    #Boltzmann probability
rot = np.zeros(3)
irot = np.zeros(3)    #rotation matrix and its inverse
torque_vec = np.zeros(3)   #torque vector
torque_sum = np.zeros(3)


#initialise Mersenne-Twister random number generator
random.seed(23435)      # make instance of random library

#initialise to the ferromagnetic state
Sx = np.zeros(pr.N_s)
Sy = np.zeros(pr.N_s)
Sz = np.ones(pr.N_s)

#initialise fields
H_app = np.zeros(3)         #effective, applied magnetic fields
H_ani = np.zeros(3)         #anisotropy field

#evaluate rotation matrix
rot, irot, mz_unit = rotate.RotMatrix(pr.theta, pr.phi)
print('mz_unit: (' + str(mz_unit[0]) + ' ' + str(mz_unit[1]) + ' ' + str(mz_unit[2]) + ')\n' )

#read in neighbours and exchange templates
print('counting no of neighbours\n' )
nFeFe = ex.CountExch('../data/effectiveJij.dat')
print('reading in ' + str(nFeFe) + ' Fe-neighbours for each of ' + str(pr.N_s) + ' Fe-sites' +'\n')
neigh_list_length = nFeFe * pr.N_s
site_central, site_neighbour, Jijno = ne.ReadNeighbours('../data/FeFe_neighbours.dat',neigh_list_length)
print('reading in ' + str(nFeFe) + ' effective exchange tensors from effectiveJij.dat' )
JijFeFe_eff, inform_FeFe, interat_FeFe = ex.ReadEffExch('../data/effectiveJij.dat', nFeFe)


#convert from Rydbergs to Joules/mu
for i in range(0, nFeFe):
    for j in range(0,9):
        JijFeFe_eff[i,j] = JijFeFe_eff[i,j] * pr.Ry



print('starting constrained Monte Carlo algorithm\n')


exch_energy_p=0.0   #exchange energies for primary spin, unprimed and primed
exch_energy_pp=0.0
exch_energy_c=0.0   #exchange energies for secondary spin, unprimed and primed
exch_energy_cp=0.0


with open('../data/magn-vs-T.dat', 'w+') as Tfile:
    for iT in range(1,1002,10):
        temp = float(iT)
        beta = 1.0 / (pr.kB*temp)
        beta_scaled = beta*pr.mu
        mz = 0.0
        mzp = 0.0

        print('at temperature: ' + str(temp) + '\n')
        print('beta = ' + str(beta) + '\n')

        #initialise spins
        Sx = np.zeros(pr.N_s)
        Sy = np.zeros(pr.N_s)
        Sz = np.ones(pr.N_s)

        sp = np.zeros(3)
        rsp = np.zeros(3)
        spp = np.zeros(3)
        rspp = np.zeros(3)
        sc = np.zeros(3)
        rsc = np.zeros(3)
        scp = np.zeros(3)
        rscp = np.zeros(3)

        #calculate mz
        for j in range(0, pr.N_s):
            mz = mz + ((Sx[j] * mz_unit[0]) + (Sy[j] * mz_unit[1]) + (Sz[j] * mz_unit[2]))

        print('mz = ' + str(mz) + '\n')
        summa=0.0

        #*----------------------!
        #      TIME LOOP        !
        #----------------------*/
        for t in range(0, N_iter):
            n_flipped=0
            n_unflipped=0

            for j in range(0, pr.N_s):
                #if j%1000==0:
                #    print(str(j) + ': ' + str(mz/float(pr.N_s)) + '\n')

                #choose a random spin
                nsp = random.randint(0, pr.N_s-1)   #primary spin
                nsc = random.randint(0, pr.N_s-1)   #compensation spin

                #primary spin cannot be the same as the compensation spin
                if(nsp==nsc):
                    continue

                if ((nsp >= pr.N_s) or (nsp < 0) or (nsc >= pr.N_s) or (nsc < 0)):
                    print('skips spin ' + nsp )
                    print('something is wrong with the RNO' )
                    print('terminating program ' )

                #load spins into temp variables
                sp[0] = Sx[nsp]
                sp[1] = Sy[nsp]
                sp[2] = Sz[nsp]
                sc[0] = Sx[nsc]
                sc[1] = Sy[nsc]
                sc[2] = Sz[nsc]

                #calculate mz
                mz=0.0
                for v in range(0, pr.N_s):
                    mz = mz + (Sx[v] * mz_unit[0]) + (Sy[v] * mz_unit[1]) + (Sz[v] * mz_unit[2])

                #rotate spins into z-direction
                rsp[0] = (rot[0][0] * sp[0]) + (rot[0][1] * sp[1]) + (rot[0][2] * sp[2])
                rsp[1] = (rot[1][0] * sp[0]) + (rot[1][1] * sp[1]) + (rot[1][2] * sp[2])
                rsp[2] = (rot[2][0] * sp[0]) + (rot[2][1] * sp[1]) + (rot[2][2] * sp[2])
                rsc[0] = (rot[0][0] * sc[0]) + (rot[0][1] * sc[1]) + (rot[0][2] * sc[2])
                rsc[1] = (rot[1][0] * sc[0]) + (rot[1][1] * sc[1]) + (rot[1][2] * sc[2])
                rsc[2] = (rot[2][0] * sc[0]) + (rot[2][1] * sc[1]) + (rot[2][2] * sc[2])

                #sample uniform sphere
                v1 = 0.0
                v2 = 0.0
                s = 2.0
                while s>1.0:
                    v1 = 2.0 * random.gauss(0,1) - 1.0
                    v2 = 2.0 * random.gauss(0,1) - 1.0
                    s = (v1 * v1) + (v2 * v2)
                ss = math.sqrt(1.0 - s)

                #random vector on unit sphere
                slump_rand[0] = 2.0 * v1 * ss
                slump_rand[1] = 2.0 * v2 * ss
                slump_rand[2] = 1.0 - (2.0 * s)

                #add random vector to primary spin
                rspp[0] = slump_rand[0]
                rspp[1] = slump_rand[1]
                rspp[2] = slump_rand[2]

                #normalise the primary spin
                norm = 1.0 / (math.sqrt(rspp[0]*rspp[0]+rspp[1]*rspp[1]+rspp[2]*rspp[2]))
                rspp[0]=rspp[0]*norm
                rspp[1]=rspp[1]*norm
                rspp[2]=rspp[2]*norm

                #adjust the compensation spin's x and y components to preserve mx=my=0
                rscp[0]=rsc[0]+rsp[0]-rspp[0]
                rscp[1]=rsc[1]+rsp[1]-rspp[1]

                #check whether to make the move or not (square root condn)
                sarg=1.0-rscp[0]*rscp[0]-rscp[1]*rscp[1] #-rscp[2]*rscp[2]
                if sarg<0:
                    continue

                #z-component sign function
                if rsc[2]<0:
                    rscp[2] = -math.sqrt(sarg)
                else:
                    rscp[2] = math.sqrt(sarg)

                #new z-component of magnetisation
                mzp=(mz+rspp[2]+rscp[2]-rsp[2]-rsc[2])

                if (mzp<0):
                    continue

                #rotate spins back into non-rotated frame
                spp[0] = (irot[0][0] * rspp[0]) + (irot[0][1] * rspp[1]) + (irot[0][2] * rspp[2])
                spp[1] = (irot[1][0] * rspp[0]) + (irot[1][1] * rspp[1]) + (irot[1][2] * rspp[2])
                spp[2] = (irot[2][0] * rspp[0]) + (irot[2][1] * rspp[1]) + (irot[2][2] * rspp[2])
                scp[0] = (irot[0][0] * rscp[0]) + (irot[0][1] * rscp[1]) + (irot[0][2] * rscp[2])
                scp[1] = (irot[1][0] * rscp[0]) + (irot[1][1] * rscp[1]) + (irot[1][2] * rscp[2])
                scp[2] = (irot[2][0] * rscp[0]) + (irot[2][1] * rscp[1]) + (irot[2][2] * rscp[2])

                #*------------------------!
                # CALCULATING OLD ENERGY  !
                #------------------------*/
                e_p=0.0
                e_pp=0.0
                e_c=0.0
                e_cp=0.0
                #primary exchange
                exch_energy_p=0.0
                for k in range(0, nFeFe):
                    neighb_site=site_neighbour[nsp*nFeFe+k]
                    Jij_line=Jijno[nsp*nFeFe+k]-1
                    exch_energy_p = (exch_energy_p
                        + JijFeFe_eff[Jij_line,0] * Sx[neighb_site]*sp[0]
                        + JijFeFe_eff[Jij_line,3] * Sy[neighb_site]*sp[0]
                        + JijFeFe_eff[Jij_line,6] * Sz[neighb_site]*sp[0]
                        + JijFeFe_eff[Jij_line,1] * Sx[neighb_site]*sp[1]
                        + JijFeFe_eff[Jij_line,4] * Sy[neighb_site]*sp[1]
                        + JijFeFe_eff[Jij_line,7] * Sz[neighb_site]*sp[1]
                        + JijFeFe_eff[Jij_line,2] * Sx[neighb_site]*sp[2]
                        + JijFeFe_eff[Jij_line,5] * Sy[neighb_site]*sp[2]
                        + JijFeFe_eff[Jij_line,8] * Sz[neighb_site]*sp[2])

                    if JijFeFe_eff[Jij_line,8]==0:
                        print('Have you remembered to copy over the correct neighbour list and exchange template?\n' )
                        print('Check headers in neighbour list and exchange template.\n' )

                e_p=e_p-exch_energy_p

                #compensation exchange
                exch_energy_c=0.0
                for k in range(0, nFeFe):
                    neighb_site=site_neighbour[nsc*nFeFe+k]
                    Jij_line=Jijno[nsc*nFeFe+k]-1
                    exch_energy_c = (exch_energy_c
                        + JijFeFe_eff[Jij_line,0]*Sx[neighb_site] * sc[0]
                        + JijFeFe_eff[Jij_line,3]*Sy[neighb_site] * sc[0]
                        + JijFeFe_eff[Jij_line,6]*Sz[neighb_site] * sc[0]
                        + JijFeFe_eff[Jij_line,1]*Sx[neighb_site] * sc[1]
                        + JijFeFe_eff[Jij_line,4]*Sy[neighb_site] * sc[1]
                        + JijFeFe_eff[Jij_line,7]*Sz[neighb_site] * sc[1]
                        + JijFeFe_eff[Jij_line,2]*Sx[neighb_site] * sc[2]
                        + JijFeFe_eff[Jij_line,5]*Sy[neighb_site] * sc[2]
                        + JijFeFe_eff[Jij_line,8]*Sz[neighb_site] * sc[2])

                e_c = e_c - exch_energy_c

                #primary Zeeman
                e_p = e_p - (H_app[0] * sp[0]) - (H_app[1] * sp[1]) - (H_app[2] * sp[2])

                #secondary zeeman
                e_c = e_c - (H_app[0] * sc[0]) - (H_app[1]*sc[1]) - (H_app[2]*sc[2])

                #primary uniaxial anisotropy
                e_p = e_p - (H_ani[2] * sp[2] * sp[2])

                #secondary uniaxial anisotropy
                e_c = e_c - (H_ani[2] * sp[2] * sp[2])

                #----------------------------!
                #     PRIMARY SPIN MOVE      !
                #----------------------------!
                Sx[nsp]=spp[0]
                Sy[nsp]=spp[1]
                Sz[nsp]=spp[2]

                #----------------------------!
                #         NEW ENERGY         !
                #----------------------------!
                #primary exchange
                exch_energy_pp=0.0
                for k in range(0, nFeFe):
                    neighb_site=site_neighbour[nsp*nFeFe+k]
                    Jij_line=Jijno[nsp*nFeFe+k]-1
                    exch_energy_pp = (exch_energy_pp
                        + JijFeFe_eff[Jij_line,0] * Sx[neighb_site] * spp[0]
                        + JijFeFe_eff[Jij_line,3] * Sy[neighb_site] * spp[0]
                        + JijFeFe_eff[Jij_line,6] * Sz[neighb_site] * spp[0]
                        + JijFeFe_eff[Jij_line,1] * Sx[neighb_site] * spp[1]
                        + JijFeFe_eff[Jij_line,4] * Sy[neighb_site] * spp[1]
                        + JijFeFe_eff[Jij_line,7] * Sz[neighb_site] * spp[1]
                        + JijFeFe_eff[Jij_line,2] * Sx[neighb_site] * spp[2]
                        + JijFeFe_eff[Jij_line,5] * Sy[neighb_site] * spp[2]
                        + JijFeFe_eff[Jij_line,8] * Sz[neighb_site] * spp[2])
                e_pp=e_pp-exch_energy_pp

                #primary zeeman
                e_pp=e_pp-(H_app[0]*spp[0])-(H_app[1]*spp[1])-(H_app[2]*spp[2])

                #primary uniaxial anisotropy
                e_pp=e_pp-(H_ani[2]*spp[2]*spp[2])

                #-------------------------------!
                #    COMPENSATION SPIN MOVE     !
                #-------------------------------!
                Sx[nsc]=scp[0]
                Sy[nsc]=scp[1]
                Sz[nsc]=scp[2]

                #------------------------------!
                #          NEW ENERGY          !
                #------------------------------!
                #compensation exchange
                exch_energy_cp=0.0
                for k in range(0, nFeFe):
                    neighb_site=site_neighbour[nsp*nFeFe+k]
                    Jij_line=Jijno[nsp*nFeFe+k]-1
                    exch_energy_cp = (exch_energy_cp
                        + JijFeFe_eff[Jij_line,0] * Sx[neighb_site] * scp[0]
                        + JijFeFe_eff[Jij_line,3] * Sy[neighb_site] * scp[0]
                        + JijFeFe_eff[Jij_line,6] * Sz[neighb_site] * scp[0]
                        + JijFeFe_eff[Jij_line,1] * Sx[neighb_site] * scp[1]
                        + JijFeFe_eff[Jij_line,4] * Sy[neighb_site] * scp[1]
                        + JijFeFe_eff[Jij_line,7] * Sz[neighb_site] * scp[1]
                        + JijFeFe_eff[Jij_line,2] * Sx[neighb_site] * scp[2]
                        + JijFeFe_eff[Jij_line,5] * Sy[neighb_site] * scp[2]
                        + JijFeFe_eff[Jij_line,8] * Sz[neighb_site] * scp[2])
                e_cp=e_cp-exch_energy_cp

                #compensation Zeeman
                e_cp=e_cp-(H_app[0]*scp[0])-(H_app[1]*scp[1])-(H_app[2]*scp[2])

                #compensation uniaxial anisotropy
                e_cp=e_cp-(H_ani[2]*scp[2]*scp[2])

                #---------------------------------!
                #   CALCULATE ENERGY DIFFERENCE   !
                #---------------------------------!
                deltaE = e_pp - e_p + e_cp - e_c
                #print('printing: ' + deltaE + '   ' + exch_energy_p + '   ' + exch_energy_pp )

                #if energy has increased, calculate Boltzmann probability
                if deltaE>0:
                    pboltz = (mzp/mz) * (mzp/mz) * math.fabs(rsc[2]/rscp[2]) * math.exp(-deltaE * beta_scaled)
                    prand = random.uniform(0,1)

                    #if random no > Boltzmann prob, revert spins to original state
                    if prand>pboltz:
                        Sx[nsp]=sp[0]
                        Sy[nsp]=sp[1]
                        Sz[nsp]=sp[2]

                        Sx[nsc]=sc[0]
                        Sy[nsc]=sc[1]
                        Sz[nsc]=sc[2]

                        n_unflipped += 1

                n_flipped += 1
            #end of trial-loop

            #find total magnetisation
            M_x=0.0
            M_y=0.0
            M_z=0.0
            for v in range(0, pr.N_s):
                M_x = M_x + Sx[v]
                M_y = M_y + Sy[v]
                M_z = M_z + Sz[v]

            if t > N_eq:
                summa = summa + (math.sqrt((M_x*M_x)+(M_y*M_y)+(M_z*M_z)) / (float(pr.N_s)))

        #end of time-loop

        print(str(pboltz) + ' ' +  str(prand) + ' ' + str(mz) + ' ' +  str(mzp) + ' ' + str(rsc[2]) + ' ' + str(rscp[2]) + ' ' + str(deltaE) + ' ' + str(beta_scaled) + '\n')

        #M_tot=M_tot/dble(N_iter-N_eq)
        summa = summa / (float(N_iter - N_eq))

        #calculate mz
        mz=0.0
        for j in range(0, pr.N_s):
            mz = mz + ((Sx[j] * mz_unit[0]) + (Sy[j] * mz_unit[1]) + (Sz[j] * mz_unit[2])) # /pr.N_s

        Tfile.write(str(temp) + ', ' + str(summa) + ', ' + str(mz) + '\n')
        summa=0.0
    #end of temperature loop


print('constrained Monte Carlo algorithm finished, exiting program\n' )
  

with open('../data/spin-configuration.dat', 'w+') as spipr:
    for i in range(0, pr.N_s):
        spipr.write(str(i) + ', ' + str(Sz[i]))


