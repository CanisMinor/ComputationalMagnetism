#pragma once

#include <array>

namespace constrained_monte_carlo
{
int RotMatrix(double theta, double phi, std::array<double, 9> irot, std::array<double, 9> rot, std::array<double, 3> mz_unit);
}
