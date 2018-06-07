#include <cmath>
#include <random>
#include <vector>

namespace constrained_monte_carlo
{
std::array<double, 3> sample_random_unit_vector()
{
    const double mean = 0.0;
    const double variance = 1.0;
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(mean, variance);
    double s = 2.0;
    double v1 = 0.0;
    double v2 = 0.0;
    while (s > 1.0)
    {
        v1 = 2.0 * distribution(generator) - 1.0;
        v2 = 2.0 * distribution(generator) - 1.0;
        s = v1 * v1 + v2 * v2;
    }

    const double sqrt_s = std::sqrt(1.0 - s);

    // random vector on unit sphere
    std::array<double, 3> random_unit_vector = {2.0 * v1 * sqrt_s, 2.0 * v2 * sqrt_s, 1.0 - 2.0 * s};
    return random_unit_vector;
}
}
