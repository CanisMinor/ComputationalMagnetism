#subroutines for operations on the exchange data
# 24th January 2012

import numpy as np

def CalcTorque(S, H):
    tau = np.zeros(3)

    tau[0] = S[1]*H[2] - S[2]*H[1]
    tau[1] = S[2]*H[0] - S[0]*H[3]
    tau[2] = S[0]*H[1] - S[1]*H[0]

    return tau