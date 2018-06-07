// constrained Monte Carlo algorithm
// Canis Minor
// York 13th Dec 2011

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <fstream>
#include <cstdio>
#include <time.h>
#include <vector>

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
int main()
{
    // parameters for time- and temperature loops
    const int N_iter = 2000;  // no of time steps
    const int N_eq = 1500;    // equilibration iteration treshold

    // constrained Monte Carlo variables
    double mz_unit[3];     // unit vector along z
    double slump_rand[3];              // random number
    double sarg = 0.0;                 // square root argument
    double deltaE = 0.0;               // change in energy
    double prand = 0.0, pboltz = 0.0;  // random number, Boltzmann probability
    double rot[3][3], irot[3][3];      // rotation matrix and its inverse

    // initialise Mersenne-Twister random number generator
    int seed = (int)time(0);   // random seed
    StochasticLib1 sto(seed);  // make instance of random library

    // initialise fields
    std::array<double, 3> H_app = {0.0, 0.0, 0.0};  // effective, applied magnetic fields
    std::array<double, 3> H_ani = {0.0, 0.0, 0.0};  //anisotropic field

    if (RotMatrix(theta, phi, irot, rot, mz_unit) != 0)
    {
        cout << "error in evaluation of rotation matrix" << endl;
        cout << "terminating program" << endl;
        return 1;
    }

    const int nFeFe = CountExch("../data/effectiveJij.dat");
    const int neigh_list_length = nFeFe * N_s;

    std::vector<int> site_central(neigh_list_length, 0);
    std::vector<int> site_neighbour(neigh_list_length, 0);
    std::vector<int> Jijno(neigh_list_length, 0);

    ReadNeighbours("../data/FeFe_neighbours.dat", neigh_list_length, site_central, site_neighbour, Jijno);
    std::vector<double> interatomic_FeFe(N_s, 0.0);
    std::vector<std::array<int, 5>> info_FeFe(N_s, {0.0, 0.0, 0.0, 0.0, 0.0});
    std::vector<std::array<double, 9>> JijFeFe_eff(nFeFe, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});

    ReadEffExch("../data/effectiveJij.dat", JijFeFe_eff, inform_FeFe, interat_FeFe, nFeFe);

    // convert from Rydbergs to Joules/mu
    for (int i_iron_pair = 0; i_iron_pair < nFeFe; ++i_iron_pair)
    {
        for (int dimension = 0; dimension < 9; ++dimension)
        {
            JijFeFe_eff[i_iron_pair][dimension]*= Ry;
        }
    }

    ofstream magn_temp;
    magn_temp.open("magn-vs-T.dat");
    double exch_energy_p = 0.0, exch_energy_pp = 0.0;  // exchange energies for primary spin, unprimed and primed
    double exch_energy_c = 0.0, exch_energy_cp = 0.0;  // exchange energies for secondary spin, unprimed and primed
    double beta = 0.0;
    double beta_scaled = 0.0;
    double mz = 0.0;
    double mzp = 0.0;
    double temp = 0.0;

    for (int iT = 1; iT < 1002; iT += 10)
    {
        temp = (double)(iT);
        beta = 1.0 / (kB * temp);
        beta_scaled = beta * mu;
        mz = 0.0;
        mzp = 0.0;

        cout << "at temperature: " << temp << endl;
        cout << "beta = " << beta << endl;

        // initialise spins to ferromagnetic state
        std::vector<double> Sx(N_s, 0.0);
        std::vector<double> Sy(N_s, 0.0);
        std::vector<double> Sz(N_s, 1.0);

        for (int i = 0; i < 3; i++)
            sp[i] = 0.0;
        for (int i = 0; i < 3; i++)
            rsp[i] = 0.0;
        for (int i = 0; i < 3; i++)
            spp[i] = 0.0;
        for (int i = 0; i < 3; i++)
            rspp[i] = 0.0;
        for (int i = 0; i < 3; i++)
            sc[i] = 0.0;
        for (int i = 0; i < 3; i++)
            rsc[i] = 0.0;
        for (int i = 0; i < 3; i++)
            scp[i] = 0.0;
        for (int i = 0; i < 3; i++)
            rscp[i] = 0.0;

        // calculate mz
        for (int j = 0; j < N_s; j++)
        {
            mz += ((Sx[j] * mz_unit[0]) + (Sy[j] * mz_unit[1]) + (Sz[j] * mz_unit[2]));
        }

        double summa = 0.0;
        for (int t = 1; t <= N_iter; t++)
        {
            int n_flipped = 0;
            int n_unflipped = 0;

            for (int j = 0; j < N_s; j++)
            {
                // choose a random spin
                int i_primary_spin = (int)(sto.Random() * N_s);  // primary spin
                int i_compensation_spin = (int)(sto.Random() * N_s);  // compensation spin

                // primary spin cannot be the same as the compensation spin
                if (i_primary_spin == i_compensation_spin)
                {
                    continue;
                }

                // check that the spin indices are within the number of spins
                if (i_primary_spin >= N_s || i_primary_spin < 0 || i_compensation_spin >= N_s || i_compensation_spin < 0)
                {
                    cout << "skips spin " << primary_spin << endl;
                    cout << "something is wrong with the RNO" << endl;
                    cout << "terminating program " << endl;
                    return 1;
                }

                // load spins into temp variables
                sp[0] = Sx[i_primary_spin];
                sp[1] = Sy[i_primary_spin];
                sp[2] = Sz[i_primary_spin];
                sc[0] = Sx[i_compensation_spin];
                sc[1] = Sy[i_compensation_spin];
                sc[2] = Sz[i_compensation_spin];

                // calculate mz
                double mz = 0.0;
                for (int v = 0; v < N_s; v++)
                {
                    mz += (Sx[v] * mz_unit[0]) + (Sy[v] * mz_unit[1]) + (Sz[v] * mz_unit[2]);
                }

                // rotate spins into z-direction
                rsp[0] = (rot[0][0] * sp[0]) + (rot[0][1] * sp[1]) + (rot[0][2] * sp[2]);
                rsp[1] = (rot[1][0] * sp[0]) + (rot[1][1] * sp[1]) + (rot[1][2] * sp[2]);
                rsp[2] = (rot[2][0] * sp[0]) + (rot[2][1] * sp[1]) + (rot[2][2] * sp[2]);
                rsc[0] = (rot[0][0] * sc[0]) + (rot[0][1] * sc[1]) + (rot[0][2] * sc[2]);
                rsc[1] = (rot[1][0] * sc[0]) + (rot[1][1] * sc[1]) + (rot[1][2] * sc[2]);
                rsc[2] = (rot[2][0] * sc[0]) + (rot[2][1] * sc[1]) + (rot[2][2] * sc[2]);

                // sample uniform sphere (Marsaglia)
                std::array<double, 3> slump_rand = sample_random_unit_vector();

                // add random vector to primary spin
                rspp[0] = slump_rand[0];
                rspp[1] = slump_rand[1];
                rspp[2] = slump_rand[2];

                // normalise the primary spin
                norm = 1.0 / (sqrt(rspp[0] * rspp[0] + rspp[1] * rspp[1] + rspp[2] * rspp[2]));
                rspp[0] = rspp[0] * norm;
                rspp[1] = rspp[1] * norm;
                rspp[2] = rspp[2] * norm;

                // adjust the compensation spin's x and y components to preserve mx=my=0
                rscp[0] = rsc[0] + rsp[0] - rspp[0];
                rscp[1] = rsc[1] + rsp[1] - rspp[1];

                // check whether to make the move or not (square root condn)
                sarg = 1.0 - rscp[0] * rscp[0] - rscp[1] * rscp[1];  //-rscp[2]*rscp[2];
                if (sarg < 0)
                {
                    continue;
                }

                // z-component sign function
                if (rsc[2] < 0)
                {
                    rscp[2] = -sqrt(sarg);
                }
                else
                {
                    rscp[2] = sqrt(sarg);
                }

                // new z-component of magnetisation
                double mz_new = (mz_old + rspp[2] + rscp[2] - rsp[2] - rsc[2]);
                if (mz_new < 0)
                {
                    continue;
                }

                // rotate spins back into non-rotated frame
                spp[0] = (irot[0][0] * rspp[0]) + (irot[0][1] * rspp[1]) + (irot[0][2] * rspp[2]);
                spp[1] = (irot[1][0] * rspp[0]) + (irot[1][1] * rspp[1]) + (irot[1][2] * rspp[2]);
                spp[2] = (irot[2][0] * rspp[0]) + (irot[2][1] * rspp[1]) + (irot[2][2] * rspp[2]);
                scp[0] = (irot[0][0] * rscp[0]) + (irot[0][1] * rscp[1]) + (irot[0][2] * rscp[2]);
                scp[1] = (irot[1][0] * rscp[0]) + (irot[1][1] * rscp[1]) + (irot[1][2] * rscp[2]);
                scp[2] = (irot[2][0] * rscp[0]) + (irot[2][1] * rscp[1]) + (irot[2][2] * rscp[2]);

                /*------------------------!
                // CALCULATING OLD ENERGY  !
                //------------------------*/
                e_p = 0.0;
                e_pp = 0.0;
                e_c = 0.0;
                e_cp = 0.0;
                // primary exchange
                exch_energy_p = 0.0;
                for (int k = 0; k < nFeFe; k++)
                {
                    int neighb_site = site_neighbour[primary_spin * nFeFe + k];
                    int Jij_line = Jijno[primary_spin * nFeFe + k] - 1;
                    exch_energy_p = exch_energy_p + JijFeFe_eff[Jij_line][0] * Sx[neighb_site] * sp[0] +
                                    JijFeFe_eff[Jij_line][3] * Sy[neighb_site] * sp[0] +
                                    JijFeFe_eff[Jij_line][6] * Sz[neighb_site] * sp[0] +
                                    JijFeFe_eff[Jij_line][1] * Sx[neighb_site] * sp[1] +
                                    JijFeFe_eff[Jij_line][4] * Sy[neighb_site] * sp[1] +
                                    JijFeFe_eff[Jij_line][7] * Sz[neighb_site] * sp[1] +
                                    JijFeFe_eff[Jij_line][2] * Sx[neighb_site] * sp[2] +
                                    JijFeFe_eff[Jij_line][5] * Sy[neighb_site] * sp[2] +
                                    JijFeFe_eff[Jij_line][8] * Sz[neighb_site] * sp[2];

                    if (JijFeFe_eff[Jij_line][8] == 0)
                    {
                        cout << "error: neighb_site, primary_spin, j, Jij_line are " << neighb_site << ", " << primary_spin << ", " << j
                             << ", " << Jij_line << ", " << primary_spin * nFeFe + k << endl;
                        cout << "Have you remembered to copy over the correct neighbour list and exchange template? "
                             << endl;
                        cout << "Check headers in neighbour list and exchange template." << endl;
                        return 1;
                    }
                }

                e_p = e_p - exch_energy_p;

                // compensation exchange
                exch_energy_c = 0.0;
                for (int k = 0; k < nFeFe; k++)
                {
                    int neighb_site = site_neighbour[i_compensation_spin * nFeFe + k];
                    int Jij_line = Jijno[i_compensation_spin * nFeFe + k] - 1;
                    exch_energy_c = exch_energy_c + JijFeFe_eff[Jij_line][0] * Sx[neighb_site] * sc[0] +
                                    JijFeFe_eff[Jij_line][3] * Sy[neighb_site] * sc[0] +
                                    JijFeFe_eff[Jij_line][6] * Sz[neighb_site] * sc[0] +
                                    JijFeFe_eff[Jij_line][1] * Sx[neighb_site] * sc[1] +
                                    JijFeFe_eff[Jij_line][4] * Sy[neighb_site] * sc[1] +
                                    JijFeFe_eff[Jij_line][7] * Sz[neighb_site] * sc[1] +
                                    JijFeFe_eff[Jij_line][2] * Sx[neighb_site] * sc[2] +
                                    JijFeFe_eff[Jij_line][5] * Sy[neighb_site] * sc[2] +
                                    JijFeFe_eff[Jij_line][8] * Sz[neighb_site] * sc[2];
                }
                e_c = e_c - exch_energy_c;

                // primary zeeman
                e_p = e_p - H_app[0] * sp[0] - H_app[1] * sp[1] - H_app[2] * sp[2];

                // compensation zeeman
                e_c = e_c - H_app[0] * sc[0] - H_app[1] * sc[1] - H_app[2] * sc[2];

                // primary uniaxial anisotropy
                e_p = e_p - (H_ani[2] * sp[2] * sp[2]);

                // compensation uniaxial anisotropy
                e_c = e_c - (H_ani[2] * sp[2] * sp[2]);

                //----------------------------!
                //     PRIMARY SPIN MOVE      !
                //----------------------------!
                Sx[primary_spin] = spp[0];
                Sy[primary_spin] = spp[1];
                Sz[primary_spin] = spp[2];

                //----------------------------!
                //         NEW ENERGY         !
                //----------------------------!
                // primary exchange
                exch_energy_pp = 0.0;
                for (int k = 0; k < nFeFe; k++)
                {
                    int neighb_site = site_neighbour[i_primary_spin * nFeFe + k];
                    int Jij_line = Jijno[i_primary_spin * nFeFe + k] - 1;
                    exch_energy_pp = exch_energy_pp + JijFeFe_eff[Jij_line][0] * Sx[neighb_site] * spp[0] +
                                     JijFeFe_eff[Jij_line][3] * Sy[neighb_site] * spp[0] +
                                     JijFeFe_eff[Jij_line][6] * Sz[neighb_site] * spp[0] +
                                     JijFeFe_eff[Jij_line][1] * Sx[neighb_site] * spp[1] +
                                     JijFeFe_eff[Jij_line][4] * Sy[neighb_site] * spp[1] +
                                     JijFeFe_eff[Jij_line][7] * Sz[neighb_site] * spp[1] +
                                     JijFeFe_eff[Jij_line][2] * Sx[neighb_site] * spp[2] +
                                     JijFeFe_eff[Jij_line][5] * Sy[neighb_site] * spp[2] +
                                     JijFeFe_eff[Jij_line][8] * Sz[neighb_site] * spp[2];
                }
                e_pp = e_pp - exch_energy_pp;

                // primary zeeman
                e_pp = e_pp - (H_app[0] * spp[0]) - (H_app[1] * spp[1]) - (H_app[2] * spp[2]);

                // primary uniaxial anisotropy
                e_pp = e_pp - (H_ani[2] * spp[2] * spp[2]);

                //-------------------------------!
                //    COMPENSATION SPIN MOVE     !
                //-------------------------------!
                Sx[i_compensation_spin] = scp[0];
                Sy[i_compensation_spin] = scp[1];
                Sz[i_compensation_spin] = scp[2];

                //------------------------------!
                //          NEW ENERGY          !
                //------------------------------!
                // compensation exchange
                exch_energy_cp = 0.0;
                for (int k = 0; k < nFeFe; k++)
                {
                    int neighb_site = site_neighbour[primary_spin * nFeFe + k];
                    int Jij_line = Jijno[primary_spin * nFeFe + k] - 1;
                    exch_energy_cp = exch_energy_cp + JijFeFe_eff[Jij_line][0] * Sx[neighb_site] * scp[0] +
                                     JijFeFe_eff[Jij_line][3] * Sy[neighb_site] * scp[0] +
                                     JijFeFe_eff[Jij_line][6] * Sz[neighb_site] * scp[0] +
                                     JijFeFe_eff[Jij_line][1] * Sx[neighb_site] * scp[1] +
                                     JijFeFe_eff[Jij_line][4] * Sy[neighb_site] * scp[1] +
                                     JijFeFe_eff[Jij_line][7] * Sz[neighb_site] * scp[1] +
                                     JijFeFe_eff[Jij_line][2] * Sx[neighb_site] * scp[2] +
                                     JijFeFe_eff[Jij_line][5] * Sy[neighb_site] * scp[2] +
                                     JijFeFe_eff[Jij_line][8] * Sz[neighb_site] * scp[2];
                }
                e_cp = e_cp - exch_energy_cp;

                // compensation Zeeman
                e_cp = e_cp - (H_app[0] * scp[0]) - (H_app[1] * scp[1]) - (H_app[2] * scp[2]);

                // compensation uniaxial anisotropy
                e_cp = e_cp - (H_ani[2] * scp[2] * scp[2]);
                
                // if energy has increased, calculate Boltzmann probability
                const double deltaE = e_pp - e_p + e_cp - e_c;
                if (deltaE > 0)
                {
                    const double pboltz = (mz_new / mz_old) * (mz_new / mz_old) * std::abs(rsc[2] / rscp[2]) * exp(-deltaE * beta_scaled);
                    const double prand = sto.Random();

                    // if random no > Boltzmann prob, revert spins to original state
                    if (prand > pboltz)
                    {
                        Sx[primary_spin] = sp[0];
                        Sy[primary_spin] = sp[1];
                        Sz[primary_spin] = sp[2];

                        Sx[compensation_spin] = sc[0];
                        Sy[compensation_spin] = sc[1];
                        Sz[compensation_spin] = sc[2];

                        n_unflipped++;
                    }
                }
                n_flipped++;
            }  // trial-loop

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

        summa = summa / (double)(N_iter - N_eq);

        // calculate mz
        mz = 0.0;
        for (int j = 0; j < N_s; j++)
        {
            mz = mz + ((Sx[j] * mz_unit[0]) + (Sy[j] * mz_unit[1]) + (Sz[j] * mz_unit[2]));  // /N_s;
        }

        magn_temp << temp << ", " << summa << ", " << M_x << ", " << M_y << ", " << M_z << ", " 
                  << M_z / N_s << ", " << mz / N_s << endl;  // pboltz,prand,deltaE !,beta_scaled!mz_ratio,sc_ratio
    }  // temperature loop
    magn_temp.close();
    cout << "constrained Monte Carlo algorithm finished, exiting program" << endl;

    ofstream spin_structure;
    spin_structure.open("spin-configuration.dat");
    for (int i = 0; i < N_s; i++)
    {
        spin_structure << i << "\t" << Sz[i] << endl;
    }
    spin_structure.close();

    return 0;
}
}  // namespace constrained_monte_carlo
