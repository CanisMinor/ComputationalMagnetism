#include <iostream>
#include <cmath>
#include "Globals.h"

std::array<double, 9> get_rotation_matrix(double theta, double phi, std::array<double, 9>& irot, 
                                          std::array<double, 9>& rot, std::array<double, 3>& mz_unit)
{
    // convert angles to radians
    const double half_turn_degrees = 180.0;
    const double theta_rad = pi * theta / half_turn_degrees;
    const double phi_rad = pi * phi / half_turn_degrees;

    // calculate trig functions
    const double cost = cos(theta_rad);
    const double sint = sin(theta_rad);
    const double cosp = cos(phi_rad);
    const double sinp = sin(phi_rad);

    // spin components (normalised)
    mz_unit[0] = cost * sinp;
    mz_unit[1] = sint * sinp;
    mz_unit[2] = cosp;

    // rotation matrix
    rot[0] = cost * cosp;
    rot[1] = -sint;
    rot[2] = cost * sinp;
    rot[3] = sint * cosp;
    rot[4] = cost;
    rot[5] = sint * sinp;
    rot[6] = -sinp;
    rot[7] = 0.0;
    rot[8] = cosp;

    // inverse rotation matrix
    irot[0] = cost * cosp;
    irot[1] = sint * cosp;
    irot[2] = -sinp;
    irot[3] = -sint;
    irot[4] = cost;
    irot[5] = 0.0;
    irot[6] = cost * sinp;
    irot[7] = sint * sinp;
    irot[8] = cosp;

    return 0;
}
