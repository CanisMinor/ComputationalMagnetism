#parameters for one-dimensional Hubbard model
n_x = 2                #no of sites along x-direction
N_s = n_x              #total no of sites in unit cell
t = 1                #nearest-neighbour kinetic energy integral
U = 0                #Hubbard U
a_x = 1              #lattice parameter
k_latx = 0           #initialising wave vector along x
epsi = 0             #energy of electron occupying site
N_e = 1                #no of electrons (band is half-filled)
C = 0                #Coulomb repulsion integral
delta_U = 0.1        #Hubbard U increment size
delta_k = 0.05       #wave vector increment size
n_av_up = 0          #average spin-up density initialised to zero
n_av_down = 0        #average spin-down density initialised to zero
n_av_prod = 0        #product of average spin-up and spin-down densities initialised to zero