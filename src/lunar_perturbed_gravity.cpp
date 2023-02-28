/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     June 2016
*  LAST MODIFIED:    April 2022
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Computes gravity using the variable fidelity force approximations
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

#include "lunar_perturbed_gravity.h"
#include "const.h"
#include "GRGM1200b.h"
#include "c_functions.h"
#include "matrix_loader.h"
#include "radial_gravity.h"

#define debug_grav 0
#define debug_grav_itr 0

 void lunar_perturbed_gravity(double t, double* Xo, double err, int i, int M, double deg, int hot, double* G, double tol, int* itr, double* Feval, IterCounters& ITRs, double* del_G){

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


   //////////////////////////////// J2-J6 ///////////////////////////////
   if (err > 1.0e-1){
     // J2-J6
     if (debug_grav == 1){
       if (i==1){
        printf("J2-J6\n");
       }
     }
     lunar_Grav_Approx(t,Xo,G,Feval);
     if (i==M+1){
        ITR1  = *itr;
        ITR2  = *itr;
        ITR3  = *itr;
        ITR4  = *itr;
        if (debug_grav_itr == 1){
           printf("ITR1 %i\n",ITR1);
        }
     }
   }


   //////////////////////////////// 1e-1 1e-4 ///////////////////////////////
   else if (err <= 1.0e-1 && err > 1.0e-4 && ITR1 == *itr - 1){ // 1e-1 1e-4
     // FULL Gravity
     if (debug_grav == 1){
       if (i==1){
        printf("Full Gravity 1\n");
       }
     }
     lunar_Grav_Full(t,Xo,G,tol,deg,Feval);
     lunar_Grav_Approx(t,Xo,Gapprox,Feval);
     for (int j=0; j<=2; j++){
       del_G[ID2(i,j+1,Nmax+1)] = G[j] - Gapprox[j];
     }
     MODEL = MODEL + 3;
   } else if (err <= 1.0e-1 && err > 1.0e-4){ // 1e-1 1e-4
     // Approximate Gravity
     if (debug_grav == 1){
       if (i==1){
        printf("Approx Gravity 1\n");
       }
     }
     lunar_Grav_Approx(t,Xo,Gapprox,Feval);
     for (int j=0; j<=2; j++){
        G[j] = Gapprox[j] + del_G[ID2(i,j+1,Nmax+1)];
     }
     if (i==M+1){
       ITR2 = *itr;
       if (debug_grav_itr == 1){
         printf("ITR2 %i\n",ITR2);
       }
     }
   }


   //////////////////////////////// 1e-4 1e-7 ///////////////////////////////
   // else if (err <= 1.0e-4 && err > 1.0e-7 && ITR2 == *itr - 1){ // 1e-4 1e-7
   //   // FULL Gravity
   //   if (debug_grav == 1){
   //     if (i==1){
   //      printf("Full Gravity 2\n");
   //     }
   //   }
   //   lunar_Grav_Full(t,Xo,G,tol,deg,Feval);
   //   Grav_Approx(t,Xo,Gapprox,Feval);
   //   for (int j=0; j<=2; j++){
   //      del_G[ID2(i,j+1,Nmax+1)] = G[j] - Gapprox[j];
   //   }
   // }
   else if (err <= 1.0e-4 && err > 1.0e-7){ // 1e-4 1e-7
     // Approximate Gravity
     if (debug_grav == 1){
       if (i==1){
        printf("Approx Gravity 2\n");
       }
     }
     lunar_Grav_Approx(t,Xo,Gapprox,Feval);
     for (int j=0; j<=2; j++){
       G[j] = Gapprox[j] + del_G[ID2(i,j+1,Nmax+1)];
     }
     if (i==M+1){
        ITR3 = *itr;
        if (debug_grav_itr == 1){
          printf("ITR3 %i\n",ITR3);
        }
     }
   }


   //////////////////////////////// 1e-7 1e-10 ///////////////////////////////
   else if (err <= 1.0e-7 && err > 1.0e-10 && ITR3 == *itr - 1){ // 1e-7 1e-10
     // FULL Gravity
     if (debug_grav == 1){
       if (i==1){
        printf("Full Gravity 3\n");
       }
     }
     lunar_Grav_Full(t,Xo,G,tol,deg,Feval);
     lunar_Grav_Approx(t,Xo,Gapprox,Feval);
     for (int j=0; j<=2; j++){
        del_G[ID2(i,j+1,Nmax+1)] = G[j] - Gapprox[j];
     }
     MODEL = MODEL + 3;
   }
   else if (err <= 1.0e-7 && err > 1.0e-10){ // 1e-7 1e-10
     // Approximate Gravity
     if (debug_grav == 1){
       if (i==1){
        printf("Approx Gravity 3\n");
       }
     }
     lunar_Grav_Approx(t,Xo,Gapprox,Feval);
     for (int j=0; j<=2; j++){
       G[j] = Gapprox[j] + del_G[ID2(i,j+1,Nmax+1)];
     }
     if (i==M+1){
        ITR4 = *itr;
        if (debug_grav_itr == 1){
          printf("ITR4 %i",ITR4);
        }
     }
   }


   //////////////////////////////// 1e-10 1e-12 ///////////////////////////////
   // else if (err <= 1.0e-10 && err > 1.0e-12 && ITR4 == *itr - 1){ // 1e-10 1e-12
   else if (MODEL == 3 && err > tol){ // 1e-10 1e-12
     // FULL Gravity
     if (debug_grav == 1){
       if (i==1){
        printf("Full Gravity 4\n");
       }
     }
     lunar_Grav_Full(t,Xo,G,tol,deg,Feval);
     lunar_Grav_Approx(t,Xo,Gapprox,Feval);
     for (int j=0; j<=2; j++){
        del_G[ID2(i,j+1,Nmax+1)] = G[j] - Gapprox[j];
     }
     MODEL = MODEL + 3;
   }
   else if (err > tol){
     // Approximate Gravity
     if (debug_grav == 1){
       if (i==1){
        printf("Approx Gravity 4\n");
       }
     }
     lunar_Grav_Approx(t,Xo,Gapprox,Feval);
     for (int j=0; j<=2; j++){
       G[j] = Gapprox[j] + del_G[ID2(i,j+1,Nmax+1)];
     }
   }

   if (i==M+1){
      *itr = *itr + 1;
   }

//store iteration coutners in struct
 ITRs.ITR1=ITR1;
 ITRs.ITR2=ITR2;
 ITRs.ITR3=ITR3;
 ITRs.ITR4=ITR4;
 ITRs.MODEL=MODEL;
 return;
 }

//FIXME: Need to change to lunar spherical harmonic coefficicents that best approximate masscons
 void lunar_Grav_Approx(double t, double* X, double* dX, double* Feval){
  int deg = 4;
  GRGM1200b(X, dX, 4);

   Feval[1] = Feval[1] + 1.0;
   return;
 }

 void lunar_Grav_Full(double t, double* Xo, double* acc, double tol, double deg, double* Feval){

   double state[6]  = {0.0};
   double dstate[6] = {0.0};

   for (int jj=0; jj<=5; jj++){
     if (jj>2){
       state[jj] = 0.0;
     }
     else{
       state[jj] = Xo[jj];
     }
   }

   double grav = 0.0;
   //FIXME: Radial gravity is written for Earth model assuming that atmospheric drag perturbation dominates at low orbit
   radial_gravity(Xo,tol,deg,&grav);
   GRGM1200b(state, &dstate[3], grav);
   Feval[0] = Feval[0] + pow(grav,2)/pow(deg,2);

   for (int jj=0; jj<=2; jj++){
     acc[jj] = dstate[jj+3];
   }

   return;
 }
