#pragma once

namespace constrained_monte_carlo
{
int ReadEffExch(std::string file_name, std::vector<std::array<double, 9>> Jij_array, std::vector<std::array<int, 5>> inform, std::vector<double> interat, int n_interactions);
}
