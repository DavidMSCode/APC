/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Jun 2022
*  LAST MODIFIED:    Jun 2022
*  AFFILIATION:      Department of Aerospace Engineering, University of Illinois, Urbana-Champaign
*  DESCRIPTION:      Convert modified equinoctial elements to cartesion position and velocity
*
* INPUT:
*    p, f, g, h, k, L -- Modified Equinoctial Elements
*
* OUTPUTS:
*    r_eci -- Position (km)
*    v_eci -- Velocity (km/s)
*/

#include <math.h>
// #include <stdio.h>
// #include <string.h>
// #include <stdlib.h>
#include "const.h"
#include "c_functions.h"

void mee2rv(double* mee, double* r_eci, double* v_eci ){
    
    double p = mee[0];
    double f = mee[1];
    double g = mee[2];
    double h = mee[3];
    double k = mee[4];
    double L = mee[5];

  // Radius
  double radius = p/(1.0 + f*cos(L) + g*sin(L));

  // Common terms
  double alpha2 = pow(h,2.0) - pow(k,2.0);
  double tani2s = pow(h,2.0) + pow(k,2.0);
  double s2     = 1.0 + tani2s;

  // ECI Position and Velocity
  r_eci[0] = radius*(cos(L) + alpha2*cos(L) + 2.0*h*k*sin(L))/s2;
  r_eci[1] = radius*(sin(L) - alpha2*sin(L) + 2.0*h*k*cos(L))/s2;
  r_eci[2] = 2.0*radius*(h*sin(L) - k*cos(L))/s2;
  v_eci[0] = -sqrt(C_MU/p)*(sin(L) + alpha2*sin(L) - 2.0*h*k*cos(L) + g - 2.0*f*h*k + alpha2*g)/s2;
  v_eci[1] = -sqrt(C_MU/p)*(-cos(L) + alpha2*cos(L) + 2.0*h*k*sin(L) - f + 2.0*g*h*k + alpha2*f)/s2;
  v_eci[2] = 2.0*sqrt(C_MU/p)*(h*cos(L) + k*sin(L) + f*h + g*k)/s2;

}
