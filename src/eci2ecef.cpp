/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Jun 2016
*  LAST MODIFIED:    Jun 2016
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Convert states from inertial frame to body frame
*
* INPUT:
*    t   -- time (s)
*    X   -- ECI Position (km)
*    V   -- ECI Velocity (km/s)
*
* OUTPUTS:
*    xB  -- ECEF Position (km)
*    vB  -- ECEF Velocity (km/s)
*/
#include <math.h>
#include "eci2ecef.h"
#include "const.h"

void eci2ecef(double t, double* X, double* V, double* xB, double* vB){
  //FIXME: Transform based on the desired frames. Don't assume near Earth
  double th       = t*C_omega;
  double cos_th   = cos(th);
  double sin_th   = sin(th);
  double vB_temp[3];
  // Convert to Position to Rotating Frame
  xB[0]   =  cos_th*X[0] + sin_th*X[1];
  xB[1]   = -sin_th*X[0] + cos_th*X[1];
  xB[2]   = X[2];

  // find velocity in rotating frame (aligned with inertial axes)
  if (V[0]+V[1]+V[2] != 0.0){
    vB_temp[0]   = V[0] + C_omega*X[1];
    vB_temp[1]   = V[1] - C_omega*X[0];
    vB_temp[2]   = V[2];
  }
  if (V[0]+V[1]+V[2] == 0.0){
    vB_temp[0] = 0.0;
    vB_temp[1] = 0.0;
    vB_temp[2] = 0.0;
  }
  // Now rotate the relative velocity into the body frame
  vB[0]   =  cos_th*vB_temp[0] + sin_th*vB_temp[1];
  vB[1]   = -sin_th*vB_temp[0] + cos_th*vB_temp[1];
  vB[2]   = vB_temp[2];
  return;
}
