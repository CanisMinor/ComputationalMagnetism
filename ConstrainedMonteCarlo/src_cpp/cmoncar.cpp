//constrained Monte Carlo algorithm
//Canis Minor
//York 13th Dec 2011

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


#ifndef MULTIFILE_PROJECT
   #include "mersenne.cpp"             // code for random number generator
   #include "stoc1.cpp"                // random library source code
   #include "userintf.cpp"             // define system specific user interface
#endif

using namespace std;

int main()
{
 //parameters for time- and temperature loops 
 const int N_iter=2000;         //no of time steps
 const int N_eq=1500;           //equilibration iteration treshold 

 //constrained Monte Carlo variables
 double M_tot=0.0;      //magnetisation
 double M_x=0.0, M_y=0.0, M_z=0.0;    //components of magnetisation    
 double summa=0.0;
 double mz_unit[3];     //unit vector along z
 int nsp=0, nsc=0;      //primary and secondary spin site indices
 double sp[3], spp[3], rsp[3], rspp[3]; 
 double sc[3], scp[3], rsc[3], rscp[3]; 
 double norm=0.0;            //normalisation temp variable
 double e_p=0.0, e_pp=0.0;
 double e_c=0.0, e_cp=0.0;
 double v1=0.0, v2=0.0, ss=0.0, s=2.0;          
 double slump_rand[3];            //random number
 double sarg=0.0;                 //square root argument
 double deltaE=0.0;               //change in energy
 double prand=0.0, pboltz=0.0;    //random number, Boltzmann probability
 double rot[3][3], irot[3][3];    //rotation matrix and its inverse
 vector <double> torque_vec(3);   //torque vector
 vector <double> torque_sum(3);

 //-x, -y and -z components for each spin 
 double Sx[N_s];
 double Sy[N_s];
 double Sz[N_s];

 //initialise Mersenne-Twister random number generator
 int seed = (int)time(0);       // random seed
 StochasticLib1 sto(seed);      // make instance of random library

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
 cout << "* Constrained Monte Carlo Simulation for an FePt System with KKR Jij Input *" << endl;
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
 int nFeFe=CountExch("../data/effectiveJij.dat");
 cout << "reading in " << nFeFe << " Fe-neighbours for each of " << N_s << " Fe-sites" << endl;
 int neigh_list_length=(nFeFe)*N_s;

 int site_central[neigh_list_length];
 int site_neighbour[neigh_list_length];
 int Jijno[neigh_list_length];

 ReadNeighbours("../data/FeFe_neighbours.dat",neigh_list_length,site_central,site_neighbour,Jijno);
 cout << endl;

 double interat_FeFe[N_s];
 int inform_FeFe[N_s][5];
 double JijFeFe_eff[nFeFe][9];

 cout << "reading in " << nFeFe << " effective exchange tensors from effectiveJij.dat" << endl; 
 ReadEffExch("../data/effectiveJij.dat", JijFeFe_eff, inform_FeFe, interat_FeFe, nFeFe);
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
/////////////////////////////////////////////////////////////!
 cout << "starting constrained Monte Carlo algorithm" << endl; 

 ofstream magn_temp;
 magn_temp.open("magn-vs-T.dat");


 double exch_energy_p=0.0, exch_energy_pp=0.0;    //exchange energies for primary spin, unprimed and primed
 double exch_energy_c=0.0, exch_energy_cp=0.0;    //exchange energies for secondary spin, unprimed and primed
 double beta=0.0; 
 double beta_scaled=0.0;
 double mz=0.0;
 double mzp=0.0;
 double temp=0.0;

                                                                                      
 for(int iT=1; iT<1002; iT+=10)
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
    for(int i=0; i<3; i++) rsp[i]=0.0;
    for(int i=0; i<3; i++) spp[i]=0.0;
    for(int i=0; i<3; i++) rspp[i]=0.0;
    for(int i=0; i<3; i++) sc[i]=0.0;
    for(int i=0; i<3; i++) rsc[i]=0.0;
    for(int i=0; i<3; i++) scp[i]=0.0;
    for(int i=0; i<3; i++) rscp[i]=0.0;

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

          //choose a random spin
          nsp=(int)(sto.Random()*N_s);   //primary spin
          nsc=(int)(sto.Random()*N_s);   //compensation spin
  
          //primary spin cannot be the same as the compensation spin
          if(nsp==nsc)
          {
             continue;
          }
          if ((nsp >= N_s) || (nsp < 0) || (nsc >= N_s) || (nsc < 0))
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
          sc[0]=Sx[nsc];
          sc[1]=Sy[nsc];
          sc[2]=Sz[nsc];

          //calculate mz
          mz=0.0;
          for(int v=0; v<N_s; v++)
          {
             mz=mz+(Sx[v]*mz_unit[0])+(Sy[v]*mz_unit[1])+(Sz[v]*mz_unit[2]);
          }
          
          //rotate spins into z-direction
          rsp[0]=(rot[0][0]*sp[0])+(rot[0][1]*sp[1])+(rot[0][2]*sp[2]);
          rsp[1]=(rot[1][0]*sp[0])+(rot[1][1]*sp[1])+(rot[1][2]*sp[2]);
          rsp[2]=(rot[2][0]*sp[0])+(rot[2][1]*sp[1])+(rot[2][2]*sp[2]);
          rsc[0]=(rot[0][0]*sc[0])+(rot[0][1]*sc[1])+(rot[0][2]*sc[2]);
          rsc[1]=(rot[1][0]*sc[0])+(rot[1][1]*sc[1])+(rot[1][2]*sc[2]);
          rsc[2]=(rot[2][0]*sc[0])+(rot[2][1]*sc[1])+(rot[2][2]*sc[2]);
          
          //sample uniform sphere
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
  
          //add random vector to primary spin
          rspp[0]=slump_rand[0];
          rspp[1]=slump_rand[1];
          rspp[2]=slump_rand[2];
  
          //normalise the primary spin
          norm=1.0/(sqrt(rspp[0]*rspp[0]+rspp[1]*rspp[1]+rspp[2]*rspp[2]));
          rspp[0]=rspp[0]*norm;
          rspp[1]=rspp[1]*norm;
          rspp[2]=rspp[2]*norm;
  
          //adjust the compensation spin's x and y components to preserve mx=my=0
          rscp[0]=rsc[0]+rsp[0]-rspp[0];
          rscp[1]=rsc[1]+rsp[1]-rspp[1];
  
          //check whether to make the move or not (square root condn)
          sarg=1.0-rscp[0]*rscp[0]-rscp[1]*rscp[1]; //-rscp[2]*rscp[2];
          if(sarg<0)
          {
             continue;
          }

          //z-component sign function
          if(rsc[2]<0) 
          {
             rscp[2]=-sqrt(sarg);
          }
          else
          {
             rscp[2]=sqrt(sarg);
          }
  
          //new z-component of magnetisation 
          mzp=(mz+rspp[2]+rscp[2]-rsp[2]-rsc[2]);
  
          if (mzp<0)
          {
             continue;
          }
          //rotate spins back into non-rotated frame
          spp[0]=(irot[0][0]*rspp[0])+(irot[0][1]*rspp[1])+(irot[0][2]*rspp[2]);
          spp[1]=(irot[1][0]*rspp[0])+(irot[1][1]*rspp[1])+(irot[1][2]*rspp[2]);
          spp[2]=(irot[2][0]*rspp[0])+(irot[2][1]*rspp[1])+(irot[2][2]*rspp[2]);
          scp[0]=(irot[0][0]*rscp[0])+(irot[0][1]*rscp[1])+(irot[0][2]*rscp[2]);
          scp[1]=(irot[1][0]*rscp[0])+(irot[1][1]*rscp[1])+(irot[1][2]*rscp[2]);
          scp[2]=(irot[2][0]*rscp[0])+(irot[2][1]*rscp[1])+(irot[2][2]*rscp[2]);
 
          /*------------------------!
          // CALCULATING OLD ENERGY  !
          //------------------------*/
          e_p=0.0;
          e_pp=0.0; 
          e_c=0.0;
          e_cp=0.0;
          //primary exchange
          exch_energy_p=0.0;
          for(int k=0; k<nFeFe; k++)
          {                                                                 
             int neighb_site=site_neighbour[nsp*nFeFe+k];
             int Jij_line=Jijno[nsp*nFeFe+k]-1;
             exch_energy_p=exch_energy_p
              +JijFeFe_eff[Jij_line][0]*Sx[neighb_site]*sp[0]
              +JijFeFe_eff[Jij_line][3]*Sy[neighb_site]*sp[0]
              +JijFeFe_eff[Jij_line][6]*Sz[neighb_site]*sp[0]
              +JijFeFe_eff[Jij_line][1]*Sx[neighb_site]*sp[1]
              +JijFeFe_eff[Jij_line][4]*Sy[neighb_site]*sp[1]
              +JijFeFe_eff[Jij_line][7]*Sz[neighb_site]*sp[1]
              +JijFeFe_eff[Jij_line][2]*Sx[neighb_site]*sp[2]
              +JijFeFe_eff[Jij_line][5]*Sy[neighb_site]*sp[2]
              +JijFeFe_eff[Jij_line][8]*Sz[neighb_site]*sp[2];
              //test << "iteration " << k << " neighbour " << neighb_site << " exchange tensor no " << Jij_line << 
              //"\t" << JijFeFe_eff[Jij_line][8] << "\t" << JijFeFe_eff[Jij_line][8]*Sz[neighb_site]*sp[2] << endl;
             if(JijFeFe_eff[Jij_line][8]==0) 
             { 
                cout << "error: neighb_site, nsp, j, Jij_line are " << neighb_site << ", " << nsp << ", " << j << ", " << Jij_line << ", " << nsp*nFeFe+k << endl;
                cout << "Have you remembered to copy over the correct neighbour list and exchange template? " << endl;
                cout << "Check headers in neighbour list and exchange template." << endl;
                return 1;
             }
          }
          e_p=e_p-exch_energy_p;

          //compensation exchange
          exch_energy_c=0.0;
          for(int k=0; k<nFeFe; k++)
          {                                                                 
             int neighb_site=site_neighbour[nsc*nFeFe+k];
             int Jij_line=Jijno[nsc*nFeFe+k]-1;
             exch_energy_c=exch_energy_c
              +JijFeFe_eff[Jij_line][0]*Sx[neighb_site]*sc[0]
              +JijFeFe_eff[Jij_line][3]*Sy[neighb_site]*sc[0]
              +JijFeFe_eff[Jij_line][6]*Sz[neighb_site]*sc[0]
              +JijFeFe_eff[Jij_line][1]*Sx[neighb_site]*sc[1]
              +JijFeFe_eff[Jij_line][4]*Sy[neighb_site]*sc[1]
              +JijFeFe_eff[Jij_line][7]*Sz[neighb_site]*sc[1]
              +JijFeFe_eff[Jij_line][2]*Sx[neighb_site]*sc[2]
              +JijFeFe_eff[Jij_line][5]*Sy[neighb_site]*sc[2]
              +JijFeFe_eff[Jij_line][8]*Sz[neighb_site]*sc[2];
          }
          e_c=e_c-exch_energy_c;
 
          //primary zeeman
          e_p=e_p-H_app[0]*sp[0] - H_app[1]*sp[1] - H_app[2]*sp[2];
  
          //secondary zeeman
          e_c=e_c-H_app[0]*sc[0] - H_app[1]*sc[1] - H_app[2]*sc[2];
  
          //primary uniaxial anisotropy
          e_p=e_p-(H_ani[2]*sp[2]*sp[2]);
  
          //secondary uniaxial anisotropy
          e_c=e_c-(H_ani[2]*sp[2]*sp[2]);
        
          //----------------------------!
          //     PRIMARY SPIN MOVE      !
          //----------------------------!
          Sx[nsp]=spp[0];
          Sy[nsp]=spp[1];
          Sz[nsp]=spp[2];
  
          //----------------------------!
          //         NEW ENERGY         !
          //----------------------------!
          //primary exchange
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
 
          //primary zeeman
          e_pp=e_pp-(H_app[0]*spp[0])-(H_app[1]*spp[1])-(H_app[2]*spp[2]);
  
          //primary uniaxial anisotropy
          e_pp=e_pp-(H_ani[2]*spp[2]*spp[2]);
     
          //-------------------------------!
          //    COMPENSATION SPIN MOVE     !
          //-------------------------------!
          Sx[nsc]=scp[0];
          Sy[nsc]=scp[1];
          Sz[nsc]=scp[2];
  
          //------------------------------!
          //          NEW ENERGY          !
          //------------------------------!
          //compensation exchange
          exch_energy_cp=0.0;
          for(int k=0; k<nFeFe; k++)
          {                                                                
             int neighb_site=site_neighbour[nsp*nFeFe+k];
             int Jij_line=Jijno[nsp*nFeFe+k]-1;
             exch_energy_cp=exch_energy_cp
              +JijFeFe_eff[Jij_line][0]*Sx[neighb_site]*scp[0]
              +JijFeFe_eff[Jij_line][3]*Sy[neighb_site]*scp[0]
              +JijFeFe_eff[Jij_line][6]*Sz[neighb_site]*scp[0]
              +JijFeFe_eff[Jij_line][1]*Sx[neighb_site]*scp[1]
              +JijFeFe_eff[Jij_line][4]*Sy[neighb_site]*scp[1]
              +JijFeFe_eff[Jij_line][7]*Sz[neighb_site]*scp[1]
              +JijFeFe_eff[Jij_line][2]*Sx[neighb_site]*scp[2]
              +JijFeFe_eff[Jij_line][5]*Sy[neighb_site]*scp[2]
              +JijFeFe_eff[Jij_line][8]*Sz[neighb_site]*scp[2];
          }
          e_cp=e_cp-exch_energy_cp;        
         
          //compensation Zeeman
          e_cp=e_cp-(H_app[0]*scp[0])-(H_app[1]*scp[1])-(H_app[2]*scp[2]);
  
          //compensation uniaxial anisotropy
          e_cp=e_cp-(H_ani[2]*scp[2]*scp[2]); 
     
          //---------------------------------!
          //   CALCULATE ENERGY DIFFERENCE   !
          //---------------------------------!
          deltaE=e_pp-e_p+e_cp-e_c;
          //cout << "printing: " << deltaE << "   " << exch_energy_p << "   " << exch_energy_pp << endl;
  
          //if energy has increased, calculate Boltzmann probability
          if(deltaE>0) 
          {
             pboltz=(mzp/mz)*(mzp/mz)*abs(rsc[2]/rscp[2])*exp(-deltaE*beta_scaled);
             prand=sto.Random();
         
             //if random no > Boltzmann prob, revert spins to original state
             if (prand>pboltz) 
             {
                Sx[nsp]=sp[0];
                Sy[nsp]=sp[1];
                Sz[nsp]=sp[2];
  
                Sx[nsc]=sc[0];
                Sy[nsc]=sc[1];
                Sz[nsc]=sc[2];

                n_unflipped++;
             }
          }
          n_flipped++;
       }   //trial-loop
       M_x=0.0;
       M_y=0.0;
       M_z=0.0;
  
       for(int v=0; v<N_s; v++) M_x=M_x+Sx[v];
       for(int v=0; v<N_s; v++) M_y=M_y+Sy[v];
       for(int v=0; v<N_s; v++) M_z=M_z+Sz[v];
  
       if(t > N_eq) 
       {
          summa=summa+sqrt((M_x*M_x)+(M_y*M_y)+(M_z*M_z))/(double(N_s));
       }
    }    //time-loop

    cout << pboltz << " " <<  prand << " " << mz << " " <<  mzp << " " << rsc[2] << " " << rscp[2] << " " << deltaE << " " << beta_scaled << endl;

    //M_tot=M_tot/dble(N_iter-N_eq) 
    summa=summa/(double)(N_iter-N_eq);

    //calculate mz
    mz=0.0;
    for(int j=0; j<N_s; j++)
    {
       mz=mz+((Sx[j]*mz_unit[0])+(Sy[j]*mz_unit[1])+(Sz[j]*mz_unit[2])); // /N_s;
    }
    magn_temp << temp << "  " << summa << "   " << M_tot << "   " <<  M_x << "   " <<  M_y << "   " <<  M_z << "   " <<  M_z/N_s << "   " <<  mz/N_s << endl; // pboltz,prand,deltaE !,beta_scaled!mz_ratio,sc_ratio
    summa=0.0;
 } //temperature loop
 magn_temp.close();
 cout << "constrained Monte Carlo algorithm finished, exiting program" << endl;
  
 ofstream spin_structure;
 spin_structure.open("spin-configuration.dat");
 for(int i=0; i<N_s; i++)
 {
    spin_structure << i << "\t" << Sz[i] << endl;
 }
 spin_structure.close();

 return 0;
}
