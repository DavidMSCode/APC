/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Jun 2022
*  LAST MODIFIED:    Jun 2022
*  AFFILIATION:      Department of Aerospace Engineering, University of Illinois, Urbana-Champaign
*  DESCRIPTION:      Converts classical orbital elements to modified equinoctial elements
*
* INPUT:
*    a, e, inc, Om, w, nu -- Classical orbital elements
*
* OUTPUTS:
*    mee  -- Modified Equinoctial Elements
*/

#include <math.h>
// #include <stdio.h>
// #include <string.h>
// #include <stdlib.h>

void classical2mee( double* coe, double* mee ){

  double a = coe[0];
  double e = coe[1];
  double inc = coe[2];
  double Om = coe[3];
  double w = coe[4];
  double nu = coe[5];
      
  mee[0] = a*(1.0 - pow(e,2));
  mee[1] = e*cos(w + Om);
  mee[2] = e*sin(w + Om);
  mee[3] = tan(inc/2)*cos(Om);
  mee[4] = tan(inc/2)*sin(Om);
  mee[5] = Om + w + nu;

}
