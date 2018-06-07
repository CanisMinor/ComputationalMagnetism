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
#include "marsaglia_algorithm.hpp"
#include "parameters.h"
#include "reader_exchange_constants.h"
#include "reader_neighbour_lists.h"
#include "rotation_matrix.h"


namespace constrained_monte_carlo
{
std::array<double, 3> rotate_spin(std::array<double, 3> s, std::array<double, 9> rot)
{
    std::array<double, 3> rotated_spin = { (rot[0] * s[0]) + (rot[1] * s[1]) + (rot[2] * s[2]),
                                           (rot[3] * s[0]) + (rot[4] * s[1]) + (rot[5] * s[2]),
                                           (rot[6] * s[0]) + (rot[7] * s[1]) + (rot[8] * s[2]) };

    return rotated_spin;
}

double calculate_total_energy(std::array<double, 3> H_app, std::array<double, 3> H_ani, std::array<double, 3> spin, std::array<double, 3> neighbour_spin, std::array<double, 9> JijFeFe)
{
    double total_energy = 0.0;
    double exchange_energy = 0.0;
    for (int k = 0; k < nFeFe; k++)
    {
        int neighb_site = site_neighbour[primary_spin * nFeFe + k];
        int Jij_line = Jijno[primary_spin * nFeFe + k] - 1;
        exchange_energy += JijFeFe[0] * neighbour_spin[0] * spin[0]
                         + JijFeFe[3] * neighbour_spin[1] * spin[0]
                         + JijFeFe[6] * neighbour_spin[2] * spin[0]
                         + JijFeFe[1] * neighbour_spin[0] * spin[1]
                         + JijFeFe[4] * neighbour_spin[1] * spin[1]
                         + JijFeFe[7] * neighbour_spin[2] * spin[1]
                         + JijFeFe[2] * neighbour_spin[0] * spin[2]
                         + JijFeFe[5] * neighbour_spin[1] * spin[2]
                         + JijFeFe[8] * neighbour_spin[2] * spin[2];

        if (JijFeFe_eff[Jij_line][8] == 0)
        {
            cout << "error: neighb_site, primary_spin, j, Jij_line are " << neighb_site << ", " << primary_spin << ", "
                 << j << ", " << Jij_line << ", " << primary_spin * nFeFe + k << endl;
            cout << "Have you remembered to copy over the correct neighbour list and exchange template? " << endl;
            cout << "Check headers in neighbour list and exchange template." << endl;
            return 1;
        }
    }

    total_energy -= exchange_energy;  // exchange energy
    total_energy -= (H_app[0] * primary_spin_old[0] + H_app[1] * primary_spin_old[1] +
                     H_app[2] * primary_spin_old[2]);                        // Zeeman energy
    total_energy -= (H_ani[2] * primary_spin_old[2] * primary_spin_old[2]);  // uniaxial anisotropy energy

    return total_energy;
}

int constrained_monte_carlo()
{
    // parameters for time- and temperature loops
    const int N_iter = 2000;  // no of time steps
    const int N_eq = 1500;    // equilibration iteration treshold

    // rotation matrix
    std::array<double, 9> rot;
    std::array<double, 9> irot;
    std::array<double, 3> mz_unit;
    get_rotation_matrix(theta, phi, irot, rot, mz_unit);

    // initialise Mersenne-Twister random number generator
    int seed = (int)time(0);   // random seed
    StochasticLib1 sto(seed);  // make instance of random library

    // initialise fields
    std::array<double, 3> H_app = { 0.0, 0.0, 0.0 };  // effective, applied magnetic fields
    std::array<double, 3> H_ani = { 0.0, 0.0, 0.0 };  // anisotropic field

    const int nFeFe = CountExch("../data/effectiveJij.dat");
    const int neigh_list_length = nFeFe * N_s;

    std::vector<int> site_central(neigh_list_length, 0);
    std::vector<int> site_neighbour(neigh_list_length, 0);
    std::vector<int> Jijno(neigh_list_length, 0);

    ReadNeighbours("../data/FeFe_neighbours.dat", neigh_list_length, site_central, site_neighbour, Jijno);
    std::vector<double> interatomic_FeFe(N_s, 0.0);
    std::vector<std::array<int, 5>> info_FeFe(N_s, { 0.0, 0.0, 0.0, 0.0, 0.0 });
    std::vector<std::array<double, 9>> JijFeFe_eff(nFeFe, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 });

    ReadEffExch("../data/effectiveJij.dat", JijFeFe_eff, inform_FeFe, interat_FeFe, nFeFe);

    // convert from Rydbergs to Joules/mu
    for (int i_iron_pair = 0; i_iron_pair < nFeFe; ++i_iron_pair)
    {
        for (int dimension = 0; dimension < 9; ++dimension)
        {
            JijFeFe_eff[i_iron_pair][dimension] *= Ry;
        }
    }

    ofstream magn_temp;
    magn_temp.open("magn-vs-T.dat");
    for (int iT = 1; iT < 1002; iT += 10)
    {
        const double temp = (double)(iT);
        const double beta = 1.0 / (kB * temp);
        const double beta_scaled = beta * mu;

        cout << "at temperature: " << temp << endl;
        cout << "beta = " << beta << endl;

        // initialise spins to ferromagnetic state
        std::vector<double> Sx(N_s, 0.0);
        std::vector<double> Sy(N_s, 0.0);
        std::vector<double> Sz(N_s, 1.0);

        std::array<double, 3> primary_spin_old = { 0.0, 0.0, 0.0 };
        std::array<double, 3> rotated_primary_spin_old = { 0.0, 0.0, 0.0 };
        std::array<double, 3> primary_spin_new = { 0.0, 0.0, 0.0 };
        std::array<double, 3> rotated_primary_spin_new = { 0.0, 0.0, 0.0 };
        std::array<double, 3> compensation_spin_old = { 0.0, 0.0, 0.0 };
        std::array<double, 3> rotated_compensation_spin_old = { 0.0, 0.0, 0.0 };
        std::array<double, 3> compensation_spin_new = { 0.0, 0.0, 0.0 };
        std::array<double, 3> rotated_compensation_spin_new = { 0.0, 0.0, 0.0 };

        // calculate mz
        double mz_old = 0.0;
        for (int j = 0; j < N_s; j++)
        {
            mz_old += ((Sx[j] * mz_unit[0]) + (Sy[j] * mz_unit[1]) + (Sz[j] * mz_unit[2]));
        }

        double summa = 0.0;
        for (int t = 1; t <= N_iter; t++)
        {
            int n_flipped = 0;
            int n_unflipped = 0;

            for (int j = 0; j < N_s; j++)
            {
                // choose a random spin
                int i_primary_spin = (int)(sto.Random() * N_s);       // primary spin
                int i_compensation_spin = (int)(sto.Random() * N_s);  // compensation spin

                // primary spin cannot be the same as the compensation spin
                if (i_primary_spin == i_compensation_spin)
                {
                    continue;
                }

                // check that the spin indices are within the number of spins
                if (i_primary_spin >= N_s || i_primary_spin < 0 || i_compensation_spin >= N_s ||
                    i_compensation_spin < 0)
                {
                    cout << "skips spin " << primary_spin << endl;
                    cout << "something is wrong with the RNO" << endl;
                    cout << "terminating program " << endl;
                    return 1;
                }

                // load spins into temp variables
                primary_spin[0] = Sx[i_primary_spin];
                primary_spin[1] = Sy[i_primary_spin];
                primary_spin[2] = Sz[i_primary_spin];
                compensation_spin[0] = Sx[i_compensation_spin];
                compensation_spin[1] = Sy[i_compensation_spin];
                compensation_spin[2] = Sz[i_compensation_spin];

                // calculate mz
                double mz = 0.0;
                for (int v = 0; v < N_s; v++)
                {
                    mz += (Sx[v] * mz_unit[0]) + (Sy[v] * mz_unit[1]) + (Sz[v] * mz_unit[2]);
                }

                // rotate spins into z-direction
                rotated_primary_spin = rotate_spin(primary_spin, rot);
                rotated_compensation_spin = rotate_spin(compensation_spin, rot);

                // sample uniform sphere (Marsaglia)
                std::array<double, 3> slump_rand = sample_random_unit_vector();

                // add random vector to primary spin
                rotated_primary_spin_new[0] = slump_rand[0];
                rotated_primary_spin_new[1] = slump_rand[1];
                rotated_primary_spin_new[2] = slump_rand[2];

                // normalise the primary spin
                norm = 1.0 / (std::sqrt(rotated_primary_spin_new[0] * rotated_primary_spin_new[0] +
                                        rotated_primary_spin_new[1] * rotated_primary_spin_new[1] +
                                        rotated_primary_spin_new[2] * rotated_primary_spin_new[2]));
                rotated_primary_spin_new[0] *= norm;
                rotated_primary_spin_new[1] *= norm;
                rotated_primary_spin_new[2] *= norm;

                // adjust the compensation spin's x and y components to preserve mx=my=0
                rotated_compensation_spin_new[0] =
                    rotated_compensation_spin_old[0] + rotated_primary_spin_old[0] - rotated_primary_spin_new[0];
                rotated_compensation_spin_new[1] =
                    rotated_compensation_spin_old[1] + rotated_primary_spin_old[1] - rotated_primary_spin_new[1];

                // check whether to make the move or not (square root condn)
                double square_root_argument =
                    1.0 - rotated_compensation_spin_new[0] * rotated_compensation_spin_new[0] -
                    rotated_compensation_spin_new[1] * rotated_compensation_spin_new[1];  //-rscp[2]*rscp[2];
                if (square_root_argument < 0)
                {
                    continue;
                }

                // z-component sign function
                if (rotated_compensation_spin[2] < 0)
                {
                    rotated_compensation_spin_new[2] = -std::sqrt(square_root_argument);
                }
                else
                {
                    rotated_compensation_spin_new[2] = std::sqrt(square_root_argument);
                }

                // new z-component of magnetisation
                double mz_new = (mz_old + rotated_primary_spin_new[2] + rotated_compensation_spin_new[2] -
                                 rotated_primary_spin_old[2] - rotated_compensation_spin_old[2]);
                if (mz_new < 0)
                {
                    continue;
                }

                // rotate spins back into non-rotated frame
                primary_spin_new = rotate_spin(rotated_primary_spin_new, irot);
                compensation_spin_new = rotate_spin(rotated_compensation_spin_new, irot);

                // Calculate old energies for primary spin
                double primary_energy_old = 0.0;
                exch_energy_p = 0.0;
                for (int k = 0; k < nFeFe; k++)
                {
                    int neighb_site = site_neighbour[primary_spin * nFeFe + k];
                    int Jij_line = Jijno[primary_spin * nFeFe + k] - 1;
                    exch_energy_p += JijFeFe_eff[Jij_line][0] * Sx[neighb_site] * sp[0] +
                                     JijFeFe_eff[Jij_line][3] * Sy[neighb_site] * primary_spin_old[0] +
                                     JijFeFe_eff[Jij_line][6] * Sz[neighb_site] * primary_spin_old[0] +
                                     JijFeFe_eff[Jij_line][1] * Sx[neighb_site] * primary_spin_old[1] +
                                     JijFeFe_eff[Jij_line][4] * Sy[neighb_site] * primary_spin_old[1] +
                                     JijFeFe_eff[Jij_line][7] * Sz[neighb_site] * primary_spin_old[1] +
                                     JijFeFe_eff[Jij_line][2] * Sx[neighb_site] * primary_spin_old[2] +
                                     JijFeFe_eff[Jij_line][5] * Sy[neighb_site] * primary_spin_old[2] +
                                     JijFeFe_eff[Jij_line][8] * Sz[neighb_site] * primary_spin_old[2];

                    if (JijFeFe_eff[Jij_line][8] == 0)
                    {
                        cout << "error: neighb_site, primary_spin, j, Jij_line are " << neighb_site << ", "
                             << primary_spin << ", " << j << ", " << Jij_line << ", " << primary_spin * nFeFe + k
                             << endl;
                        cout << "Have you remembered to copy over the correct neighbour list and exchange template? "
                             << endl;
                        cout << "Check headers in neighbour list and exchange template." << endl;
                        return 1;
                    }
                }

                primary_energy_old -= exch_energy_p;  // exchange energy
                primary_energy_old -= (H_app[0] * primary_spin_old[0] + H_app[1] * primary_spin_old[1] +
                                       H_app[2] * primary_spin_old[2]);  // Zeeman energy
                primary_energy_old -=
                    (H_ani[2] * primary_spin_old[2] * primary_spin_old[2]);  // uniaxial anisotropy energy

                // Calculate old energy for compensation spin
                double compensation_energy_old = 0.0;
                double exch_energy_c = 0.0;
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

                compensation_energy_old -= exch_energy_c;
                compensation_energy_old -= (H_app[0] * compensation_spin_old[0] + H_app[1] * compensation_spin_old[1] +
                                            H_app[2] * compensation_spin_old[2]);  // Zeeman energy
                compensation_energy_old -= (H_ani[2] * compensation_spin_old[2] * compensation_spin_old[2]);

                // move primary spin
                Sx[primary_spin] = primary_spin_new[0];
                Sy[primary_spin] = primary_spin_new[1];
                Sz[primary_spin] = primary_spin_new[2];

                // calculate new energy for primary spin
                double primary_energy_new = 0.0;
                double exchange_energy_primary_new = 0.0;
                for (int k = 0; k < nFeFe; k++)
                {
                    int neighb_site = site_neighbour[i_primary_spin * nFeFe + k];
                    int Jij_line = Jijno[i_primary_spin * nFeFe + k] - 1;
                    exch_enexchange_energy_primary_newergy_pp += JijFeFe_eff[Jij_line][0] * Sx[neighb_site] * primary_spin_new[0] +
                                                                 JijFeFe_eff[Jij_line][3] * Sy[neighb_site] * primary_spin_new[0] +
                                                                 JijFeFe_eff[Jij_line][6] * Sz[neighb_site] * primary_spin_new[0] +
                                                                 JijFeFe_eff[Jij_line][1] * Sx[neighb_site] * primary_spin_new[1] +
                                                                 JijFeFe_eff[Jij_line][4] * Sy[neighb_site] * primary_spin_new[1] +
                                                                 JijFeFe_eff[Jij_line][7] * Sz[neighb_site] * primary_spin_new[1] +
                                                                 JijFeFe_eff[Jij_line][2] * Sx[neighb_site] * primary_spin_new[2] +
                                                                 JijFeFe_eff[Jij_line][5] * Sy[neighb_site] * primary_spin_new[2] +
                                                                 JijFeFe_eff[Jij_line][8] * Sz[neighb_site] * primary_spin_new[2];
                }

                primary_energy_new -= exchange_energy_primary_new;
                primary_energy_new -= ((H_app[0] * primary_spin_new[0]) + (H_app[1] * primary_spin_new[1]) + (H_app[2] * primary_spin_new[2]));  // Zeeman
                primary_energy_new -= (H_ani[2] * primary_spin_new[2] * primary_spin_new[2]);  // uniaxial anisotropy


                // move compensation spin
                Sx[i_compensation_spin] = compensation_spin_new[0];
                Sy[i_compensation_spin] = compensation_spin_new[1];
                Sz[i_compensation_spin] = compensation_spin_new[2];

                // calculate new compensation energy
                double compensation_energy_new = 0.0;
                double exchange_energy_compensation_new = 0.0;
                for (int k = 0; k < nFeFe; k++)
                {
                    int neighb_site = site_neighbour[primary_spin * nFeFe + k];
                    int Jij_line = Jijno[primary_spin * nFeFe + k] - 1;
                    exchange_energy_compensation_new += JijFeFe_eff[Jij_line][0] * Sx[neighb_site] * compensation_spin_new[0] +
                                                        JijFeFe_eff[Jij_line][3] * Sy[neighb_site] * compensation_spin_new[0] +
                                                        JijFeFe_eff[Jij_line][6] * Sz[neighb_site] * compensation_spin_new[0] +
                                                        JijFeFe_eff[Jij_line][1] * Sx[neighb_site] * compensation_spin_new[1] +
                                                        JijFeFe_eff[Jij_line][4] * Sy[neighb_site] * compensation_spin_new[1] +
                                                        JijFeFe_eff[Jij_line][7] * Sz[neighb_site] * compensation_spin_new[1] +
                                                        JijFeFe_eff[Jij_line][2] * Sx[neighb_site] * compensation_spin_new[2] +
                                                        JijFeFe_eff[Jij_line][5] * Sy[neighb_site] * compensation_spin_new[2] +
                                                        JijFeFe_eff[Jij_line][8] * Sz[neighb_site] * compensation_spin_new[2];
                }

                compensation_energy_new -= exchange_energy_compensation_new;
                compensation_energy_new -= ((H_app[0] * compensation_spin_new[0]) - (H_app[1] * compensation_spin_new[1]) - (H_app[2] * compensation_spin_new[2]));  // Zeeman
                compensation_energy_new -= (H_ani[2] * compensation_spin_new[2] * compensation_spin_new[2]);  // uniaxial anisotropy

                // if energy has increased, calculate Boltzmann probability
                const double deltaE = compensation_energy_new - compensation_energy_old + primary_energy_new - primary_energy_old;
                if (deltaE > 0)
                {
                    const double pboltz =
                        (mz_new / mz_old) * (mz_new / mz_old) * std::abs(rotated_compensation_spin_old[2] / rotated_compensation_spin_new[2]) * exp(-deltaE * beta_scaled);
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

        summa = summa / (double)(N_iter - N_eq);

        // calculate mz
        double mz = 0.0;
        for (int j = 0; j < N_s; j++)
        {
            mz += ((Sx[j] * mz_unit[0]) + (Sy[j] * mz_unit[1]) + (Sz[j] * mz_unit[2]));  // /N_s;
        }

        magn_temp << temp << ", " << summa << ", " << M_x << ", " << M_y << ", " << M_z << ", " << M_z / N_s << ", "
                  << mz / N_s << endl;  // pboltz,prand,deltaE !,beta_scaled!mz_ratio,sc_ratio
    }                                   // temperature loop

    magn_temp.close();

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
