/*
*  AUTHORS:          David Stanley (DavidMS4@illinois.edu) based on code by Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     August 2022
*  LAST MODIFIED:    August 2022
*  AFFILIATION:      Department of Aerospace Engineering, University of Illinois Urbana Champaign
*  DESCRIPTION:      Computes gravity using pre computed gravity perturbations of a nearby reference orbit
*
* INPUT:
*    t     -- Time (s)
*    Xo    -- State (position and velocity)
*    err   -- Picard iteration relative error
*    i     -- index along array for current segment
*    M     -- Number of sample points
*    deg   -- Gravity degree
*    hot   -- Hot start switch condition
*    tol   -- Tolerance
*    itr   -- Picard iteration counter
*    Feval -- Function evaluation counter
*    del_G -- Array containing precomputed gravity perturbations
*
* OUTPUTS:
*    G   -- Gravitational acceleration
*
* REFERENCES:
* 1. Macomber, B., Probe, A., Woollands, R., Read, J., and Junkins, J., "Enhancements of
*    Modified Chebyshev Picard Iteration for Perturbed Orbit Propagation", Computational
*    Modelling in Engineering & Sciences, Vol. 111, pp, 29-64, 2016.
*/
#include <math.h>
#include <cstdio>

#include "c_functions.h"
#include "const.h"
#include "EGM2008.h"
#include "matrix_loader.h"
#include "radial_gravity.h"
#include "perturbed_gravity.h"
#include "offset_gravity.h"

#define debug_offset 0
#define debug_offset_itr 0

 void offset_gravity(double t, double* Xo, double err, int i, int M, double deg, int hot, double* G, double tol, int* itr, double* Feval, IterCounters& ITRs, double* del_G){
  double Gapprox[3] = {0.0};
  //retrieve iteration counter values
  int ITR1=ITRs.ITR1;
  int ITR2=ITRs.ITR2;
  int ITR3=ITRs.ITR3;
  int ITR4=ITRs.ITR4;
  int MODEL=ITRs.MODEL;
   // Initialization
   if (*itr == 0 && hot == 0){
     ITR1     = 0;
     ITR2     = 0;
     ITR3     = 0;
     ITR4     = 0;
     MODEL    = 0;
   }

   // Initialization with hot start
   if (*itr == 0 && hot == 1){
     ITR1     = -1;
     ITR2     = -1;
     ITR3     = -1;
     ITR4     = -1;
     MODEL    = 0;
   }

if (err>tol)
{
    // Approximate Gravity with precomputed del_G
    if (debug_offset == 1){
      if (i==1){
      printf("Offset Gravity\n");
      }
    }
    Grav_Approx(t,Xo,Gapprox,Feval);
    for (int j=0; j<=2; j++){
      G[j] = Gapprox[j] + del_G[ID2(i,j+1,Nmax+1)];
    }
    
}
if (i==M+1){
  *itr = *itr + 1;
}
//store iteration counters in struct
 ITRs.ITR1=ITR1;
 ITRs.ITR2=ITR2;
 ITRs.ITR3=ITR3;
 ITRs.ITR4=ITR4;
 ITRs.MODEL=MODEL;
 return;
 }

 