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


#ifndef MULTIFILE_PROJECT
   #include "mersenne.cpp"             // code for random number generator
   #include "stoc1.cpp"                // random library source code
   #include "userintf.cpp"             // define system specific user interface
#endif

using namespace std;

int main()
{
 //parameters for time- and temperature loops 
 const int N_iter=20000;         //no of time steps
 const int N_eq=15000;           //equilibration iteration treshold 

 //constrained Monte Carlo variables
 double M_x=0.0, M_y=0.0, M_z=0.0;    //components of magnetisation    
 double summa=0.0;
 double mz_unit[3];     //unit vector along z
 int nsp=0;      //primary spin site index
 double sp[3], spp[3]; 
 double e_p=0.0, e_pp=0.0;
 double v1=0.0, v2=0.0, ss=0.0, s=2.0;          
 double slump_rand[3];          //random number
 double deltaE=0.0;        //change in energy
 double prand=0.0, pboltz=0.0;   //random number, Boltzmann probability
 double rot[3][3], irot[3][3];      //rotation matrix and its inverse

 //-x, -y and -z components for each spin 
 double Sx[N_s];
 double Sy[N_s];
 double Sz[N_s];

 //initialise Mersenne-Twister random number generator
 int seed = (int)time(0);            // random seed
 StochasticLib1 sto(seed);           // make instance of random library

 //initialise to the ferromagnetic state
 for(int i=0; i<N_s; i++) Sx[i]=0.0;
 for(int i=0; i<N_s; i++) Sy[i]=0.0;
 for(int i=0; i<N_s; i++) Sz[i]=1.0;

 //initialise fields
 double H_app[3];         //effective, applied magnetic fields
 H_app[0]=0.0;
 H_app[1]=0.0;
 H_app[2]=0.0;
 double H_ani[3];         //anisotropy field
 H_ani[0]=0.0;
 H_ani[1]=0.0;
 H_ani[2]=0.0;


 cout << endl; 
 cout << "****************************************************************************" << endl;
 cout << "*        Monte Carlo Simulation for an FePt System with KKR Jij Input      *" << endl;
 cout << "****************************************************************************" << endl;
 cout << endl; 
 cout << endl; 

 if(RotMatrix(theta,phi,irot,rot,mz_unit)!=0)
 {
    cout << "error in evaluation of rotation matrix" << endl;
    cout << "terminating program" << endl;
    return 1;
 }
 cout << endl;

 cout << "mz_unit: " << mz_unit[0] << "\t" << mz_unit[1] << "\t" << mz_unit[2] << endl;
 
 cout << "counting no of neighbours" << endl;
 int nFeFe=CountExch("effectiveJij.dat");
 cout << "reading in " << nFeFe << " Fe-neighbours for each of " << N_s << " Fe-sites" << endl;
 int neigh_list_length=(nFeFe)*N_s;

 int site_central[neigh_list_length];
 int site_neighbour[neigh_list_length];
 int Jijno[neigh_list_length];

 ReadNeighbours("FeFe_neighbours.dat",neigh_list_length,site_central,site_neighbour,Jijno);
 cout << endl;

 double interat_FeFe[N_s];
 int inform_FeFe[N_s][5];
 double JijFeFe_eff[nFeFe][9];

 cout << "reading in " << nFeFe << " effective exchange tensors from effectiveJij.dat" << endl; 
 ReadEffExch("effectiveJij.dat", JijFeFe_eff, inform_FeFe, interat_FeFe, nFeFe);
 cout << endl;

 //convert from Rydbergs to Joules/mu
 for(int i=0; i<nFeFe; i++)
 {
   JijFeFe_eff[i][0]=JijFeFe_eff[i][0]*Ry;
   JijFeFe_eff[i][1]=JijFeFe_eff[i][1]*Ry;
   JijFeFe_eff[i][2]=JijFeFe_eff[i][2]*Ry;
   JijFeFe_eff[i][3]=JijFeFe_eff[i][3]*Ry;
   JijFeFe_eff[i][4]=JijFeFe_eff[i][4]*Ry;
   JijFeFe_eff[i][5]=JijFeFe_eff[i][5]*Ry;
   JijFeFe_eff[i][6]=JijFeFe_eff[i][6]*Ry;
   JijFeFe_eff[i][7]=JijFeFe_eff[i][7]*Ry;
   JijFeFe_eff[i][8]=JijFeFe_eff[i][8]*Ry;
 }


 cout << "starting Monte Carlo algorithm" << endl;

// ofstream test;
// test.open("test.dat");
 ofstream magn_temp;
 magn_temp.open("magn-vs-T.dat");

 double exch_energy_p=0.0, exch_energy_pp=0.0;    //exchange energies for primary spin, unprimed and primed
 double beta=0.0; 
 double beta_scaled=0.0;
 double mz=0.0;
 double mzp=0.0;
 double temp=0.0;

                                                                                         
 for(int iT=300; iT < 700; iT=iT+5)
 {                                                                         
    temp=(double)(iT);
    beta=1.0/(kB*temp);  
    beta_scaled=beta*mu;
    mz=0.0;
    mzp=0.0;

    cout << "at temperature: " << temp << endl;
    cout << "beta = " << beta << endl;
 
    //initialise spins                                                                    
    for(int i=0; i<N_s; i++) Sx[i]=0.0;
    for(int i=0; i<N_s; i++) Sy[i]=0.0;
    for(int i=0; i<N_s; i++) Sz[i]=1.0;

    for(int i=0; i<3; i++) sp[i]=0.0;
    for(int i=0; i<3; i++) spp[i]=0.0;
       
    //calculate mz
    for(int j=0; j<N_s; j++)
    {
       mz=mz+((Sx[j]*mz_unit[0])+(Sy[j]*mz_unit[1])+(Sz[j]*mz_unit[2]));
    }

    cout << "mz = " << mz << endl;
    summa=0.0;

    /*----------------------!
    !      TIME LOOP        !
    !----------------------*/
    for(int t=1; t<=N_iter; t++)
    {
       int n_flipped=0;
       int n_unflipped=0;

       for(int j=0; j<N_s; j++)
       {
          if(j%1000==0) 
          {
             //cout << j << "   " << mz/N_s << endl;
          }
    
          //choose a random primary spin
          nsp=(int)(sto.Random()*N_s); 

          //if nsp is outside the system, cycle to the next random spin
          if ((nsp >= N_s) || (nsp < 0))
          {
             cout << "skips spin " << nsp << endl;
             cout << "something is wrong with the RNO" << endl;
             cout << "terminating program " << endl;
             return 1;
          }

          //load spins into temp variables
          sp[0]=Sx[nsp];  
          sp[1]=Sy[nsp];
          sp[2]=Sz[nsp];
        
          //sample uniform sphere (Marsaglia)
          v1=0.0;
          v2=0.0;
          s=2.0;
          ss=0.0;
          while(s>1.0)
          {
             v1=2.0*sto.Normal(0,1)-1.0;
             v2=2.0*sto.Normal(0,1)-1.0;
             s=v1*v1 + v2*v2;
          }
          ss=sqrt(1.0-s);
  
          //random vector on unit sphere
          slump_rand[0]=2.0*v1*ss;
          slump_rand[1]=2.0*v2*ss;
          slump_rand[2]=1.0-2.0*s;
          //test << "slump rand: " << slump_rand[0] << "\t" << slump_rand[1] << "\t" << slump_rand[2] << endl;
  
          //move random spin to Marsaglia direction
          spp[0]=slump_rand[0];
          spp[1]=slump_rand[1];
          spp[2]=slump_rand[2];

          double norm=1.0/sqrt((spp[0]*spp[0])+(spp[1]*spp[1])+(spp[2]*spp[2]));

          //test << "spp: " << spp[0] << "\t" << spp[1] << "\t" << spp[2] << " norm: " << norm << endl;

          spp[0]=spp[0]/norm;
          spp[1]=spp[1]/norm;
          spp[2]=spp[2]/norm;

 
          /*------------------------!
          ! CALCULATING OLD ENERGY  !
          !------------------------*/
          e_p=0.0;
          //primary exchange; calculate exchange field by summing over contributions from each neighbour    
          exch_energy_p=0.0;
          for(int k=0; k<nFeFe; k++)
          {                                                                 
             int neighb_site=site_neighbour[nsp*nFeFe+k];
             int Jij_line=Jijno[nsp*nFeFe+k]-1;
             exch_energy_p=exch_energy_p
             /* +JijFeFe_eff[Jij_line][0]*Sx[neighb_site]*sp[0]
              +JijFeFe_eff[Jij_line][3]*Sy[neighb_site]*sp[0]
              +JijFeFe_eff[Jij_line][6]*Sz[neighb_site]*sp[0]
              +JijFeFe_eff[Jij_line][1]*Sx[neighb_site]*sp[1]
              +JijFeFe_eff[Jij_line][4]*Sy[neighb_site]*sp[1]
              +JijFeFe_eff[Jij_line][7]*Sz[neighb_site]*sp[1]
              +JijFeFe_eff[Jij_line][2]*Sx[neighb_site]*sp[2]
              +JijFeFe_eff[Jij_line][5]*Sy[neighb_site]*sp[2]*/
              +JijFeFe_eff[Jij_line][8]*Sz[neighb_site]*sp[2];
              //test << "iteration " << k << " neighbour " << neighb_site << " exchange tensor no " << Jij_line << 
              //"\t" << JijFeFe_eff[Jij_line][8] << "\t" << JijFeFe_eff[Jij_line][8]*Sz[neighb_site]*sp[2] << endl;
          }
          e_p=e_p-exch_energy_p;

          //test << "exchange energy 1: " << e_p << endl;
        
          //primary Zeeman energy (old)
          e_p=e_p - H_app[0]*sp[0] - H_app[1]*sp[1] - H_app[2]*sp[2];
        
          //primary spin uniaxial anisotropy (old)
          e_p=e_p-(H_ani[2]*sp[2]*sp[2]);
     
          //move primary spin
          Sx[nsp]=spp[0];
          Sy[nsp]=spp[1];
          Sz[nsp]=spp[2];
        
          //new primary exchange; calculate exchange field by summing over contributions from each neighbour    
          e_pp=0.0;
          exch_energy_pp=0.0;
          for(int k=0; k<nFeFe; k++)
          {                                                                 
             int neighb_site=site_neighbour[nsp*nFeFe+k];
             int Jij_line=Jijno[nsp*nFeFe+k]-1;
             exch_energy_pp=exch_energy_pp
              +JijFeFe_eff[Jij_line][0]*Sx[neighb_site]*spp[0]
              +JijFeFe_eff[Jij_line][3]*Sy[neighb_site]*spp[0]
              +JijFeFe_eff[Jij_line][6]*Sz[neighb_site]*spp[0]
              +JijFeFe_eff[Jij_line][1]*Sx[neighb_site]*spp[1]
              +JijFeFe_eff[Jij_line][4]*Sy[neighb_site]*spp[1]
              +JijFeFe_eff[Jij_line][7]*Sz[neighb_site]*spp[1]
              +JijFeFe_eff[Jij_line][2]*Sx[neighb_site]*spp[2]
              +JijFeFe_eff[Jij_line][5]*Sy[neighb_site]*spp[2]
              +JijFeFe_eff[Jij_line][8]*Sz[neighb_site]*spp[2];
          }
          e_pp=e_pp-exch_energy_pp;
          //test << "exchange energy 2: " << e_pp << "\t" << spp[0] << "\t" << spp[1] << "\t" << spp[2] << endl;
        
          //new primary Zeeman energy
          e_pp=e_pp-(H_app[0]*spp[0])-(H_app[1]*spp[1])-(H_app[2]*spp[2]);
        
          //new primary uniaxial anisotropy
          e_pp=e_pp-(H_ani[2]*spp[2]*spp[2]);  
        
          //calculate total change in energy
          deltaE=(e_pp-e_p);
          //test << "out: " << nsp << "\t" << deltaE << "\t" << beta_scaled << "\t" << prand << "\t" << pboltz << endl;

          //calculate Boltzmann probability
          pboltz=min(exp(-deltaE*beta_scaled),1.00);
          prand=sto.Random();
          //test << "in: " << nsp << "\t" << deltaE << "\t" << beta_scaled << "\t" << prand << "\t" << pboltz << endl;
          
          //if random no < Boltzmann prob, revert spins to original state
          if (prand>pboltz)
          {
             Sx[nsp]=sp[0];
             Sy[nsp]=sp[1];
             Sz[nsp]=sp[2];
             n_unflipped++;
            // test << "reverts " << t << "   " << j << "   " << endl; 
            // test << "DOES NOT FLIP" << endl;
             continue;
          }
          //test << "FLIPS" << endl;
          n_flipped++;
       }
       M_x=0.0;
       M_y=0.0;
       M_z=0.0;
     
       for(int v=0; v<N_s; v++)
       {
          M_x=M_x+Sx[v];///(N_s)
       }
       M_x=M_x/N_s;
     
       for(int v=0; v<N_s; v++)
       {
          M_y=M_y+Sy[v];///(N_s)
       }
       M_y=M_y/N_s;
     
       for(int v=0; v<N_s; v++)
       {
          M_z=M_z+Sz[v];///(N_s)
       }
       M_z=M_z/N_s;

       if(t > N_eq) 
       {
          for(int v=0; v<N_s; v++)
          {
             summa=summa+sqrt((M_x*M_x)+(M_y*M_y)+(M_z*M_z)); 
          }
       }
    //test << "time step: " << t << " n_flipped: " << n_flipped << " n_unflipped: " << n_unflipped << endl;
    }  //time-loop 

    summa=summa/((double)(N_iter-N_eq));
  
    //calculate mz
    mz=0.0;
    for(int j=0; j<N_s; j++)
    {
       mz=mz+((Sx[j]*mz_unit[0])+(Sy[j]*mz_unit[1])+(Sz[j]*mz_unit[2])); //!/N_s
    }
 
    magn_temp << temp << "  " << summa << "   " <<  M_x << "   " <<  M_y << "   " <<  M_z << "   " <<  M_z/N_s << "   " <<  mz/N_s << endl; 
    summa=0.0;
 }   //temperature-loop
 magn_temp.close();
 //test.close();

 cout << "Monte Carlo algorithm finished, exiting program" << endl;

 ofstream spin_structure;
 spin_structure.open("spin-configuration.dat");
 for(int i=0; i<N_s; i++)
 {
    spin_structure << i << "\t" << Sz[i] << endl;
 }
 spin_structure.close();

 return 0;
}
