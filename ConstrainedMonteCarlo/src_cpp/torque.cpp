#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <fstream>

#include "parameters.h"

namespace constrained_monte_carlo
{
std::array<double, 3> calculate_torque(std::array<double, 3> S, std::array<double, 3>, int n_spins)
{
    std::array<double, 3> tau = { S[1] * H[2] - S[2] * H[1],
                                  S[1] * H[2] - S[2] * H[1],
                                  S[0] * H[1] - S[1] * H[0] };

    return tau;
}
}
