namespace constrained_monte_carlo
{
const int n_x = 9;
const int n_y = 9;
const int n_z = 9;
const int N_s = n_x * n_y * n_z;
const int theta = 0.0;
const int phi = 0.0;
const double alpha = 0.1;                   // Gilbert damping parameter
const double gamm = 1.76E11;                // gyromagnetic ratio
const double pi = 3.14159265358979;         // pi=3.14159265358979...
const double kB = 1.38E-23;                 // Boltzmann constant
const double mu = 9.274E-24;                // Bohr magneton
const double Ry = 13.6 * 1.6 * 1E-19 / mu;  // Rydberg in Joules, divided by the Bohr magneton
const double dz = 0.6;                      // anisotropy
const double mPt = 0.3234180633;            // Pt moment from KKR calc
}
