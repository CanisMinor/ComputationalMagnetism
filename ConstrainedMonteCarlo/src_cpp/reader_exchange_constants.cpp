//subroutines for operations on the exchange data

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <fstream>

using namespace std;

#include "Globals.h"

//count the number of exchange pairs listed in data file
int CountExch(const char *s)
{
 cout << "counting no of interactions in file " << s << '\n';

 ifstream exchfile;   //exchfile is a variable of type "input file"
 exchfile.open(s);
 if(!exchfile) return 0;

 string junk;
 int ncount=-1;
 while ( !exchfile.eof() )
 {                                                                          
    ncount++;                                                               
    getline(exchfile, junk);
 }
 exchfile.close();

 return ncount;
}



//read in exchange tensors
int ReadEffExch(const char *s, double Jij_array[][9], int inform[][5], double interat[], int n_interactions)
{
 ifstream exch_templ;
 exch_templ.open(s);
 if(!exch_templ) cout << "effective exchange tensors could not be read in" << endl;
 
 //read in exchange parameters and pairs info from KKR data
 for(int j=0; j<n_interactions; j++)
 {
    exch_templ >> inform[j][0] >> inform[j][1] >> inform[j][2] >> inform[j][3] >> inform[j][4] >> interat[j] >> Jij_array[j][0] >> Jij_array[j][1] >> Jij_array[j][2] >> Jij_array[j][3] >> Jij_array[j][4] >> Jij_array[j][5] >> Jij_array[j][6] >> Jij_array[j][7] >> Jij_array[j][8]; 
 } 


 exch_templ.close();
 return 0;
}

