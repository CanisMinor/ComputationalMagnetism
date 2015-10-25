
#include <iostream>
#include <cmath>

#include "Globals.h"

using namespace std;


//-------------------------------------------------!
//        Construction of Rotation Matrix          !
//-------------------------------------------------!

int RotMatrix(double theta, double phi, double irot[3][3], double rot[3][3], double mz_unit[3])
{ 
 //convert angles to radians
 double theta_rad=pi*theta/180.0;
 double phi_rad=pi*phi/180.0;

 //calculate trig functions
 double cost=cos(theta_rad);
 double sint=sin(theta_rad);
 double cosp=cos(phi_rad);
 double sinp=sin(phi_rad);

 //spin components (normalised)
 mz_unit[0]=cost*sinp;
 mz_unit[1]=sint*sinp;
 mz_unit[2]=cosp;

 cout << "constraint direction: " << "<" << mz_unit[0] << " " << mz_unit[1] << " " << mz_unit[2] << ">" << endl;
 cout << "evaluating rotation matrix (and its inverse) for theta = " << theta << ", phi = " << phi << endl;
 
 //rotation matrix
 rot[0][0]=cost*cosp;
 rot[0][1]=-sint;
 rot[0][2]=cost*sinp;
 rot[1][0]=sint*cosp;
 rot[1][1]=cost;
 rot[1][2]=sint*sinp;
 rot[2][0]=-sinp;
 rot[2][1]=0.0;
 rot[2][2]=cosp;

 cout << "rotation matrix: " << endl;
 cout << rot[0][0] << " " << rot[0][1] << " " << rot[0][2] << endl;
 cout << rot[1][0] << " " << rot[1][1] << " " << rot[1][2] << endl;
 cout << rot[2][0] << " " << rot[2][1] << " " << rot[2][2] << endl;


 //inverse rotation matrix
 irot[0][0]=cost*cosp;
 irot[0][1]=sint*cosp;
 irot[0][2]=-sinp;
 irot[1][0]=-sint;
 irot[1][1]=cost;
 irot[1][2]=0.0;
 irot[2][0]=cost*sinp;
 irot[2][1]=sint*sinp;
 irot[2][2]=cosp;

 cout << "inverse rotation matrix: " << endl;
 cout << irot[0][0] << " " << irot[0][1] << " " << irot[0][2] << endl;
 cout << irot[1][0] << " " << irot[1][1] << " " << irot[1][2] << endl;
 cout << irot[2][0] << " " << irot[2][1] << " " << irot[2][2] << endl;

 return 0;
}
