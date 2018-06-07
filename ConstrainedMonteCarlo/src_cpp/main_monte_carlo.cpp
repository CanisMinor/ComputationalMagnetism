/*----------------------------------------------------------------------------!
!                                                                             !
! KKR-Atomistic Model of FePt                                                 !
!                                                                             !
! This program performs the Monte Carlo algorithm for                         !
! an FePt system using KKR exchange tensor input.                             !
!                                                                             !
! Author: Canis Minor                                                         !
! Brehult, December 2011                                                      !
!                                                                             !
!----------------------------------------------------------------------------*/

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <fstream>
#include <cstdio>
#include <time.h>

#include "randomc.h"
#include "stocc.h"
#include "Globals.h"
#include "RotMat.h"
#include "ExchangeReadIn.h"
#include "NeighbReadIn.h"
#include "marsaglia_algorithm.hpp"

#ifndef MULTIFILE_PROJECT
#include "mersenne.cpp"  // code for random number generator
#include "stoc1.cpp"     // random library source code
#include "userintf.cpp"  // define system specific user interface
#endif

using namespace std;

namespace constrained_monte_carlo
{
double generate_random_number(double a, double b)
{
    static std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(a, b);
    return distribution(generator);
}

double generate_random_number(int a, int b)
{
    static std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(a, b);
    return distribution(generator);
}

int main()
{
    // parameters for time- and temperature loops
    const int N_iter = 20000;  // no of time steps
    const int N_eq = 15000;    // equilibration iteration treshold

    // constrained Monte Carlo variables
    double mz_unit[3];  // unit vector along z
    double rot[3][3], irot[3][3];      // rotation matrix and its inverse

    // initialise magnetic fields
    std::array<double, 3> H_appplied = { 0.0, 0.0, 0.0 };  // effective, applied magnetic fields
    std::array<double, 3> H_anisotropic = { 0.0, 0.0, 0.0 };

    if (RotMatrix(theta, phi, irot, rot, mz_unit) != 0)
    {
        cout << "error in evaluation of rotation matrix" << endl;
        cout << "terminating program" << endl;
        return 1;
    }

    const int nFeFe = CountExch("effectiveJij.dat");
    const int neigh_list_length = (nFeFe)*N_s;

    std::vector<int> site_central[neigh_list_length];
    std::vector<int> site_neighbour[neigh_list_length];
    std::vector<int> Jijno[neigh_list_length];

    ReadNeighbours("FeFe_neighbours.dat", neigh_list_length, site_central, site_neighbour, Jijno);

    std::vector<double> interat_FeFe(N_s);
    int inform_FeFe[N_s][5];
    double JijFeFe_eff[nFeFe][9];

    ReadEffExch("effectiveJij.dat", JijFeFe_eff, inform_FeFe, interat_FeFe, nFeFe);

    // convert from Rydbergs to Joules/mu
    for (int i_iron_site = 0; i_iron_site < nFeFe; ++i_iron_site)
    {
        for (int i_dimension = 0; i_dimension < 9; ++i_dimension)
        {
            JijFeFe_eff[i_iron_site][i_dimension] *= Ry;
        }
    }

    ofstream magn_temp;
    magn_temp.open("magn-vs-T.dat");
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_int_distribution<int> uniform_distribution(0, N_s);

    const int start_temperature = 300;
    const int stop_temperature = 700;
    const int step_temperature = 5;
    for (int i_temp = start_temperature; i_temp < stop_temperature; i_temp += step_temperature)
    {
        const double temp = (double)(i_temp);
        const double beta = 1.0 / (kB * temp);
        const double beta_scaled = beta * mu;

        // initialise to the ferromagnetic state
        std::vector<double> Sx(N_s, 0.0);
        std::vector<double> Sy(N_s, 0.0);
        std::vector<double> Sz(N_s, 1.0);

        // calculate mz
        for (int j = 0; j < N_s; j++)
        {
            mz = mz + ((Sx[j] * mz_unit[0]) + (Sy[j] * mz_unit[1]) + (Sz[j] * mz_unit[2]));
        }

        /*----------------------!
        !      TIME LOOP        !
        !----------------------*/
        double summa = 0.0;
        for (int i_time = 1; i_time <= N_iter; ++i_time)
        {
            int n_flipped = 0;
            int n_unflipped = 0;
            for (int j = 0; j < N_s; j++)
            {
                // choose a random primary spin
                int primary_spin_index = (int)(uniform_distribution(generator));

                // if nsp is outside the system, cycle to the next random spin
                if (primary_spin_index >= N_s || primary_spin_index < 0)
                {
                    std::cerr << "skips spin " << primary_spin_index << std::endl;
                    std::cerr << "something is wrong with the RNG" << std::endl;
                    std::cerr << "terminating program " << std::endl;
                    return 1;
                }

                // load spins into temp variables
                std::array<double, 3> spin_old = { Sx[nsp], Sy[nsp], Sz[nsp] };

                // sample uniform sphere (Marsaglia)
                std::array<double, 3> slump_rand = sample_random_unit_vector();

                // move random spin to Marsaglia direction
                std::array<double, 3> spin_new = { slump_rand[0], slump_rand[1], slump_rand[2] };

                const double norm =
                    1.0 / sqrt((spin_new[0] * spin_new[0]) + (spin_new[1] * spin_new[1]) + (spin_new[2] * spin_new[2]));
                spin_new[0] /= norm;
                spin_new[1] /= norm;
                spin_new[2] /= norm;

                /*------------------------!
                ! CALCULATING OLD ENERGY  !
                !------------------------*/
                double energy_old = 0.0;
                // primary exchange; calculate exchange field by summing over contributions from each neighbour
                double exch_energy_old = 0.0;
                for (int k = 0; k < nFeFe; k++)
                {
                    int neighb_site = site_neighbour[nsp * nFeFe + k];
                    int Jij_line = Jijno[nsp * nFeFe + k] - 1;
                    exch_energy_old += JijFeFe_eff[Jij_line][8] * Sz[neighb_site] * spin_old[2];
                }

                energy_old -= exch_energy_old;

                // primary Zeeman energy (old)
                energy_old -= (H_app[0] * spin_old[0] + H_app[1] * spin_old[1] + H_app[2] * spin_old[2]);

                // primary spin uniaxial anisotropy (old)
                energy_old -= H_ani[2] * spin_old[2] * spin_old[2];

                // move primary spin
                Sx[nsp] = spin_new[0];
                Sy[nsp] = spin_new[1];
                Sz[nsp] = spin_new[2];

                // new primary exchange; calculate exchange field by summing over contributions from each neighbour
                double energy_new = 0.0;
                double exch_energy_new = 0.0;
                for (int k = 0; k < nFeFe; ++k)
                {
                    int neighb_site = site_neighbour[nsp * nFeFe + k];
                    int Jij_line = Jijno[nsp * nFeFe + k] - 1;
                    exch_energy_new += JijFeFe_eff[Jij_line][0] * Sx[neighb_site] * spin_new[0]
                                     + JijFeFe_eff[Jij_line][3] * Sy[neighb_site] * spin_new[0]
                                     + JijFeFe_eff[Jij_line][6] * Sz[neighb_site] * spin_new[0]
                                     + JijFeFe_eff[Jij_line][1] * Sx[neighb_site] * spin_new[1]
                                     + JijFeFe_eff[Jij_line][4] * Sy[neighb_site] * spin_new[1]
                                     + JijFeFe_eff[Jij_line][7] * Sz[neighb_site] * spin_new[1]
                                     + JijFeFe_eff[Jij_line][2] * Sx[neighb_site] * spin_new[2]
                                     + JijFeFe_eff[Jij_line][5] * Sy[neighb_site] * spin_new[2]
                                     + JijFeFe_eff[Jij_line][8] * Sz[neighb_site] * spin_new[2];
                }

                energy_new -= exch_energy_new;

                // new primary Zeeman energy
                energy_new -= ((H_app[0] * spin_new[0]) + (H_app[1] * spin_new[1]) + (H_app[2] * spin_new[2]));

                // new primary uniaxial anisotropy
                energy_new -= H_ani[2] * spin_new[2] * spin_new[2];

                // calculate total change in energy
                double deltaE = energy_new - energy_old;

                // calculate Boltzmann probability
                const double max_probability = 1.00;
                double pboltz = std::min(exp(-deltaE * beta_scaled), max_probability);
                double prand = sto.Random();

                // if random no < Boltzmann prob, revert spins to original state
                if (prand > pboltz)
                {
                    Sx[nsp] = sp[0];
                    Sy[nsp] = sp[1];
                    Sz[nsp] = sp[2];
                    n_unflipped++;
                    continue;
                }

                n_flipped++;
            }

            std::array<double, 3> magnetisation = { 0.0, 0.0, 0.0 };
            magnetisation[0] += std::accumulate(Sx.begin(), Sx.end(), 0);  ///(N_s)
            magnetisation[1] += std::accumulate(Sy.begin(), Sy.end(), 0);
            magnetisation[2] += std::accumulate(Sz.begin(), Sz.end(), 0);

            magnetisation[0] /= (double)N_s;
            magnetisation[1] /= (double)N_s;
            magnetisation[2] /= (double)N_s;

            if (i_time > N_eq)
            {
                for (int v = 0; v < N_s; v++)
                {
                    double magnitude = 0.0;
                    for (int dimension = 0; dimension < 3; ++dimension)
                    {
                        magnitude += magnetisation[dimension] * magnetisation[dimension];
                    }

                    summa += std::sqrt(magnitude);
                }
            }
        }  // time-loop

        summa /= ((double)(N_iter - N_eq));

        // calculate mz
        double mz = 0.0;
        for (int j = 0; j < N_s; j++)
        {
            mz += ((Sx[j] * mz_unit[0]) + (Sy[j] * mz_unit[1]) + (Sz[j] * mz_unit[2]));  //!/N_s
        }

        magn_temp << temp << ", " << summa << ", " << M_tot << ", " << M_x << ", " << M_y << ", " << M_z << ", "
                  << M_z / N_s << ", " << mz / N_s << endl;
    }
    magn_temp.close();

    ofstream spin_structure;
    spin_structure.open("spin-configuration.dat");
    for (int i_spin = 0; i_spin < N_s; ++i_spin)
    {
        spin_structure << i_spin << "\t" << Sz[i_spin] << endl;
    }

    spin_structure.close();

    return 0;
}
}
