#-----------------------------------------|
#   	       HOPPING TERM		  |
#-----------------------------------------|
import numpy as np
import cmath

def hopping(t, N_s, a_x, a_y, k_latx, k_laty, neighbs, phases):
    H_t = np.zeros((N_s, N_s), dtype=complex)

    imag = complex(0.0, 1.0)

    for iSite in range(0, N_s):
        currentNeighs = neighbs[iSite,:]
        currentPhases = phases[iSite,:]
        H_t[iSite, currentNeighs[0]] = H_t[iSite, currentNeighs[0]] + t * (currentPhases[0] - 1)
        H_t[iSite, currentNeighs[1]] = H_t[iSite, currentNeighs[1]] + t * (currentPhases[1] - 1)
        H_t[iSite, currentNeighs[2]] = H_t[iSite, currentNeighs[2]] + t * (currentPhases[2] - 1)
        H_t[iSite, currentNeighs[3]] = H_t[iSite, currentNeighs[3]] + t * (currentPhases[3] - 1)
        H_t[iSite, currentNeighs[0]] = H_t[iSite, currentNeighs[0]] - (t * currentPhases[0] * cmath.exp(imag * k_latx * a_x))
        H_t[iSite, currentNeighs[1]] = H_t[iSite, currentNeighs[1]] - (t * currentPhases[1] * cmath.exp(-imag * k_latx * a_x))
        H_t[iSite, currentNeighs[2]] = H_t[iSite, currentNeighs[2]] - (t * currentPhases[2] * cmath.exp(imag * k_laty * a_y))
        H_t[iSite, currentNeighs[3]] = H_t[iSite, currentNeighs[3]] - (t * currentPhases[3] * cmath.exp(-imag * k_laty * a_y))
    return H_t


