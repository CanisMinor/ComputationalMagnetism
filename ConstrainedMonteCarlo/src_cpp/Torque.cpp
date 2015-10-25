//subroutines for operations on the exchange data
// 24th January 2012

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <fstream>


using namespace std;

#include "Globals.h"

//calculte torque for given N, S, H
vector <double> CalcTorque(double S[3], double H[3], int n_spins)
{
 vector <double> tau(3);

 tau[0] = S[1]*H[2] - S[2]H[1]
 tau[1] = S[2]*H[0] - S[0]H[3]
 tau[2] = S[0]*H[1] - S[1]H[0]

 return tau;
}

