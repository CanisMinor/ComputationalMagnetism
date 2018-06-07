// neighbour-related functions

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <fstream>
#include "parameters.h"

using namespace std;

int ReadNeighbours(const char* s, int list_length, int site_i[], int site_j[], int Jij_no[])
{
    cout << "reading in neighbours from neighbour list " << s << '\n';
    ifstream neighbfile;
    int junk_stuff[list_length];
    neighbfile.open(s);
    if (!neighbfile)
    {
        cout << "neighbour list couldn't be read in ReadNeighbours subroutine" << endl;
    }

    for (int n_line = 0; n_line < list_length; n_line++)
    {
        // read in step_i,step_j,neighbour no,tensor no
        neighbfile >> site_i[n_line] >> site_j[n_line] >> junk_stuff[n_line] >> Jij_no[n_line];
    }

    neighbfile.close();

    return 0;
}
