
#-------------------------------------------------!
#        Construction of Rotation Matrix          !
#-------------------------------------------------!
import math
import numpy as np

def RotMatrix(theta, phi):

    #initialise matrices
    rot = np.zeros(3,3)
    irot = np.zeros(3,3)
    mz_unit = np.zeros(3)

    #convert angles to radians
    theta_rad = math.pi * theta /180.0
    phi_rad = math.pi * phi /180.0

    #calculate trig functions
    cost=math.cos(theta_rad)
    sint=math.sin(theta_rad)
    cosp=math.cos(phi_rad)
    sinp=math.sin(phi_rad)

    #spin components (normalised)
    mz_unit[0]=cost*sinp
    mz_unit[1]=sint*sinp
    mz_unit[2]=cosp

    print('constraint direction: ' + '<' + str(mz_unit[0]) + ' ' + str(mz_unit[1]) + ' ' + str(mz_unit[2]) + '>' + '\n')
    print('evaluating rotation matrix (and its inverse) for theta = ' + str(theta) + ', phi = ' + str(phi) + '\n')
 
    #rotation matrix
    rot[0,0]=cost*cosp
    rot[0,1]=-sint
    rot[0,2]=cost*sinp
    rot[1,0]=sint*cosp
    rot[1,1]=cost
    rot[1,2]=sint*sinp
    rot[2,0]=-sinp
    rot[2,1]=0.0
    rot[2,2]=cosp

    print("rotation matrix: " + '\n')
    print(rot[0,0] + ' ' + rot[0,1] + ' ' + rot[0,2] + '\n')
    print(rot[1,0] + ' ' + rot[1,1] + ' ' + rot[1,2] + '\n')
    print(rot[2,0] + ' ' + rot[2,1] + ' ' + rot[2,2] + '\n')


    #inverse rotation matrix
    irot[0,0]=cost*cosp
    irot[0,1]=sint*cosp
    irot[0,2]=-sinp
    irot[1,0]=-sint
    irot[1,1]=cost
    irot[1,2]=0.0
    irot[2,0]=cost*sinp
    irot[2,1]=sint*sinp
    irot[2,2]=cosp

    print( "inverse rotation matrix: " + '\n')
    print( irot[0,0] + ' ' + irot[0,1] + ' ' + irot[0,2] + '\n')
    print( irot[1,0] + ' ' + irot[1,1] + ' ' + irot[1,2] + '\n')
    print( irot[2,0] + ' ' + irot[2,1] + ' ' + irot[2,2] + '\n')

    return rot, irot, mz_unit
