#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <array>
#include <fstream>
#include <vector>
#include "parameters.h"

namespace constrained_monte_carlo
{
int read_exchange_constants(std::string& file_name, std::vector<std::array<double, 9>>& Jij_array,
                            std::vector<std::array<int, 5>>& information, std::vector<double>& inter_atomic,
                            int n_interactions)
{
    std::ifstream exchange_constants_input(file_name.c_str());
    if (!exchange_constants_input)
    {
        std::cout << "effective exchange tensors could not be read in" << std::endl;
    }

    // read in exchange parameters and pairs info from KKR data
    std::string input_line;
    while (std::getline(std::cin, input_line))
    {
        std::array<int, 5> information_line;
        double interatomic_distance;
        std::array<double, 9> Jij_array_line;
        input_line >> information_line[0] >> information_line[1] >> information_line[2] >> information_line[3] >>
            information_line[4] >> interatomic_distance >> Jij_array_line[0] >> Jij_array_line[1] >>
            Jij_array_line[2] >> Jij_array_line[3] >> Jij_array_line[4] >> Jij_array_line[5] >> Jij_array_line[6] >>
            Jij_array_line[7] >> Jij_array_line[8];

        Jij_array.push_back(Jij_array_line);
        inter_atomic.push_back(interatomic_distance);
        information.push_back(information_line);
    }

    exchange_constants_input.close();

    return 0;
}
}  // namespace constrained_monte_carlo
