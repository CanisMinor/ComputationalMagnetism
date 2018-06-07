// neighbour-related functions

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <fstream>
#include "parameters.h"

using namespace std;

void read_neighbour_list(const std::string& file_name, std::vector<int> site_i, std::vector<int> site_j, std::vector<int> Jij_no)
{
    std::ifstream neighbour_file(file_name.c_str());
    while (std::getline(neighbour_file, input_line))
    {
        // read in step_i,step_j,neighbour no,tensor no
        int site_i_val;
        int site_j_val;
        double junk_stuff;
        int Jij_no_val;
        input_line >> site_i_val >> site_j_val >> junk_stuff >> Jij_no_val;

        site_i.push_back(site_i_val);
        site_j.push_back(site_j_val);
        Jij_no.push_back(Jij_no_val);
    }

    neighbour_file.close();
}
