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

#include "perturbed_gravity.h"
#include "const.h"
#include "EGM2008.h"
#include "c_functions.h"
#include "matrix_loader.h"
#include "radial_gravity.h"

#define debug_grav 0
#define debug_grav_itr 0

void perturbed_gravity_error(double t, double *Xo, double err, int i, int M, double deg, int hot, double *G, double tol, int *itr, double *Feval, IterCounters &ITRs, double *del_G)
{
  double& SF = ITRs.SF;
  double Gapprox[3] = {0.0};
  if(err<tol/ITRs.TOL_SCALE){
    //adjust tolerance
    if(ITRs.TOL_SCALE*SF<1){
      ITRs.TOL_SCALE = ITRs.TOL_SCALE*SF;
    }
    else{
      ITRs.TOL_SCALE = 1;
    }
    ITRs.DID_FULL = true;
    Grav_Full(t, Xo, G, tol, deg, Feval);
    Grav_Approx(t, Xo, Gapprox, Feval);
    for (int j = 0; j <= 2; j++)
    {
      del_G[ID2(i, j + 1, Nmax + 1)] = G[j] - Gapprox[j];
    }
  }
  else{
    ITRs.DID_FULL = false;
    Grav_Approx(t, Xo, Gapprox, Feval);
    for (int j = 0; j <= 2; j++)
    {
      G[j] = Gapprox[j] + del_G[ID2(i, j + 1, Nmax + 1)];
    }
  }
  return;
}

void perturbed_gravity(double t, double *Xo, double err, int i, int M, double deg, int hot, double *G, double tol, int *itr, double *Feval, IterCounters &ITRs, double *del_G)
{

  double Gapprox[3] = {0.0};
  // retrieve iteration counter values
  int ITR1 = ITRs.ITR1;
  int ITR2 = ITRs.ITR2;
  int ITR3 = ITRs.ITR3;
  int ITR4 = ITRs.ITR4;
  int MODEL = ITRs.MODEL;
  bool NEXT_FULL = ITRs.NEXT_FULL;
  // Initialization
  if (*itr == 0 && hot == 0)
  {
    ITR1 = 0;
    ITR2 = 0;
    ITR3 = 0;
    ITR4 = 0;
    MODEL = 0;
  }

  // Initialization with hot start
  if (*itr == 0 && hot == 1)
  {
    ITR1 = -1;
    ITR2 = -1;
    ITR3 = -1;
    ITR4 = -1;
    MODEL = 0;
  }

  if(NEXT_FULL){
    NEXT_FULL = false;
    { // 1e-1 1e-4
    // FULL Gravity
    if (debug_grav == 1)
    {
      if (i == 1)
      {
        printf("Forced Full Gravity 1\n");
      }
    }
    Grav_Full(t, Xo, G, tol, deg, Feval);
    Grav_Approx(t, Xo, Gapprox, Feval);
    for (int j = 0; j <= 2; j++)
    {
      del_G[ID2(i, j + 1, Nmax + 1)] = G[j] - Gapprox[j];
    }
    MODEL = MODEL + 3;
  }
  }
  //////////////////////////////// J2-J6 ///////////////////////////////
  else if (err > 1.0e-1)
  {
    // J2-J6
    if (debug_grav == 1)
    {
      if (i == 1)
      {
        printf("J2-J6\n");
      }
    }
    Grav_Approx(t, Xo, G, Feval);
    if (i == M + 1)
    {
      ITR1 = *itr;
      ITR2 = *itr;
      ITR3 = *itr;
      ITR4 = *itr;
      if (debug_grav_itr == 1)
      {
        printf("ITR1 %i\n", ITR1);
      }
    }
  }

  //////////////////////////////// 1e-1 1e-4 ///////////////////////////////
  else if (err <= 1.0e-1 && err > 1.0e-4 && ITR1 == *itr - 1)
  { // 1e-1 1e-4
    // FULL Gravity
    if (debug_grav == 1)
    {
      if (i == 1)
      {
        printf("Full Gravity 1\n");
      }
    }
    Grav_Full(t, Xo, G, tol, deg, Feval);
    Grav_Approx(t, Xo, Gapprox, Feval);
    for (int j = 0; j <= 2; j++)
    {
      del_G[ID2(i, j + 1, Nmax + 1)] = G[j] - Gapprox[j];
    }
    MODEL = MODEL + 3;
  }
  else if (err <= 1.0e-1 && err > 1.0e-4)
  { // 1e-1 1e-4
    // Approximate Gravity
    if (debug_grav == 1)
    {
      if (i == 1)
      {
        printf("Approx Gravity 1\n");
      }
    }
    Grav_Approx(t, Xo, Gapprox, Feval);
    for (int j = 0; j <= 2; j++)
    {
      G[j] = Gapprox[j] + del_G[ID2(i, j + 1, Nmax + 1)];
    }
    if (i == M + 1)
    {
      ITR2 = *itr;
      if (debug_grav_itr == 1)
      {
        printf("ITR2 %i\n", ITR2);
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
  //   Grav_Full(t,Xo,G,tol,deg,Feval);
  //   Grav_Approx(t,Xo,Gapprox,Feval);
  //   for (int j=0; j<=2; j++){
  //      del_G[ID2(i,j+1,Nmax+1)] = G[j] - Gapprox[j];
  //   }
  // }
  else if (err <= 1.0e-4 && err > 1.0e-7)
  { // 1e-4 1e-7
    // Approximate Gravity
    if (debug_grav == 1)
    {
      if (i == 1)
      {
        printf("Approx Gravity 2\n");
      }
    }
    Grav_Approx(t, Xo, Gapprox, Feval);
    for (int j = 0; j <= 2; j++)
    {
      G[j] = Gapprox[j] + del_G[ID2(i, j + 1, Nmax + 1)];
    }
    if (i == M + 1)
    {
      ITR3 = *itr;
      if (debug_grav_itr == 1)
      {
        printf("ITR3 %i\n", ITR3);
      }
    }
  }

  //////////////////////////////// 1e-7 1e-10 ///////////////////////////////
  else if (err <= 1.0e-7 && err > 1.0e-10 && ITR3 == *itr - 1)
  { // 1e-7 1e-10
    // FULL Gravity
    if (debug_grav == 1)
    {
      if (i == 1)
      {
        printf("Full Gravity 3\n");
      }
    }
    Grav_Full(t, Xo, G, tol, deg, Feval);
    Grav_Approx(t, Xo, Gapprox, Feval);
    for (int j = 0; j <= 2; j++)
    {
      del_G[ID2(i, j + 1, Nmax + 1)] = G[j] - Gapprox[j];
    }
    MODEL = MODEL + 3;
  }
  else if (err <= 1.0e-7 && err > 1.0e-10)
  { // 1e-7 1e-10
    // Approximate Gravity
    if (debug_grav == 1)
    {
      if (i == 1)
      {
        printf("Approx Gravity 3\n");
      }
    }
    Grav_Approx(t, Xo, Gapprox, Feval);
    for (int j = 0; j <= 2; j++)
    {
      G[j] = Gapprox[j] + del_G[ID2(i, j + 1, Nmax + 1)];
    }
    if (i == M + 1)
    {
      ITR4 = *itr;
      if (debug_grav_itr == 1)
      {
        printf("ITR4 %i", ITR4);
      }
    }
  }

  //////////////////////////////// 1e-10 1e-12 ///////////////////////////////
  // else if (err <= 1.0e-10 && err > 1.0e-12 && ITR4 == *itr - 1){ // 1e-10 1e-12
  else if (MODEL == 3 && err > tol)
  { // 1e-10 1e-12
    // FULL Gravity
    if (debug_grav == 1)
    {
      if (i == 1)
      {
        printf("Full Gravity 4\n");
      }
    }
    Grav_Full(t, Xo, G, tol, deg, Feval);
    Grav_Approx(t, Xo, Gapprox, Feval);
    for (int j = 0; j <= 2; j++)
    {
      del_G[ID2(i, j + 1, Nmax + 1)] = G[j] - Gapprox[j];
    }
    MODEL = MODEL + 3;
  }
  else if (err > tol)
  {
    // Approximate Gravity
    if (debug_grav == 1)
    {
      if (i == 1)
      {
        printf("Approx Gravity 4\n");
      }
    }
    Grav_Approx(t, Xo, Gapprox, Feval);
    for (int j = 0; j <= 2; j++)
    {
      G[j] = Gapprox[j] + del_G[ID2(i, j + 1, Nmax + 1)];
    }
  }

  if (i == M + 1)
  {
    *itr = *itr + 1;
  }

  // store iteration coutners in struct
  ITRs.ITR1 = ITR1;
  ITRs.ITR2 = ITR2;
  ITRs.ITR3 = ITR3;
  ITRs.ITR4 = ITR4;
  ITRs.MODEL = MODEL;
  return;
}

void Grav_EarthJ2J6(double *X, double *dX)
{
  // Return the J2-J6 Earth Gravity Evaluation with X in ECEF frame.
  double r = sqrt(pow(X[0], 2) + pow(X[1], 2) + pow(X[2], 2));

  double J2 = 1082.63e-6;
  double J3 = -2.52e-6;
  double J4 = -1.61e-6;
  double J5 = -0.15e-6;
  double J6 = 0.57e-6;

  double x_r_1;
  double y_r_1;
  double z_r_1;
  double z_r_2;
  double z_r_3;
  double z_r_4;
  double z_r_5;
  double z_r_6;

  x_r_1 = X[0] / r;
  y_r_1 = X[1] / r;
  z_r_1 = X[2] / r;
  z_r_2 = pow((X[2] / r), 2);
  z_r_3 = pow((X[2] / r), 3);
  z_r_4 = pow((X[2] / r), 4);
  z_r_5 = pow((X[2] / r), 5);
  z_r_6 = pow((X[2] / r), 6);

  double aTB[3] = {0.0};
  double aJ2[3] = {0.0};
  double aJ3[3] = {0.0};
  double aJ4[3] = {0.0};
  double aJ5[3] = {0.0};
  double aJ6[3] = {0.0};
  double aJ2J6[3] = {0.0};

  for (int i = 0; i <= 2; i++)
  {
    aTB[i] = -(C_MU / pow(r, 3)) * X[i];
  }

  aJ2[0] = -3.0 / 2.0 * J2 * (C_MU / pow(r, 2)) * pow((C_Req / r), 2) * (1.0 - 5.0 * z_r_2) * x_r_1;
  aJ2[1] = -3.0 / 2.0 * J2 * (C_MU / pow(r, 2)) * pow((C_Req / r), 2) * (1.0 - 5.0 * z_r_2) * y_r_1;
  aJ2[2] = -3.0 / 2.0 * J2 * (C_MU / pow(r, 2)) * pow((C_Req / r), 2) * (3.0 - 5.0 * z_r_2) * z_r_1;

  aJ3[0] = (1.0 / 2.0 * J3 * (C_MU / pow(r, 2)) * pow((C_Req / r), 3)) * 5.0 * (7.0 * z_r_3 - 3.0 * z_r_1) * x_r_1;
  aJ3[1] = (1.0 / 2.0 * J3 * (C_MU / pow(r, 2)) * pow((C_Req / r), 3)) * 5.0 * (7.0 * z_r_3 - 3.0 * z_r_1) * y_r_1;
  aJ3[2] = (1.0 / 2.0 * J3 * (C_MU / pow(r, 2)) * pow((C_Req / r), 3)) * 3.0 * (1.0 - 10.0 * z_r_2 + 35.0 / 3.0 * z_r_4);

  aJ4[0] = (5.0 / 8.0 * J4 * (C_MU / pow(r, 2)) * pow((C_Req / r), 4)) * (3.0 - 42.0 * z_r_2 + 63.0 * z_r_4) * x_r_1;
  aJ4[1] = (5.0 / 8.0 * J4 * (C_MU / pow(r, 2)) * pow((C_Req / r), 4)) * (3.0 - 42.0 * z_r_2 + 63.0 * z_r_4) * y_r_1;
  aJ4[2] = (5.0 / 8.0 * J4 * (C_MU / pow(r, 2)) * pow((C_Req / r), 4)) * (15.0 - 70.0 * z_r_2 + 63.0 * z_r_4) * z_r_1;

  aJ5[0] = (J5 / 8.0 * (C_MU / pow(r, 2)) * pow((C_Req / r), 5)) * 3.0 * (35.0 * z_r_1 - 210.0 * z_r_3 + 231.0 * z_r_5) * x_r_1;
  aJ5[1] = (J5 / 8.0 * (C_MU / pow(r, 2)) * pow((C_Req / r), 5)) * 3.0 * (35.0 * z_r_1 - 210.0 * z_r_3 + 231.0 * z_r_5) * y_r_1;
  aJ5[2] = (J5 / 8.0 * (C_MU / pow(r, 2)) * pow((C_Req / r), 5)) * (693.0 * z_r_6 - 945.0 * z_r_4 + 315.0 * z_r_2 - 15.0);

  aJ6[0] = (-J6 / 16.0 * (C_MU / pow(r, 2)) * pow((C_Req / r), 6)) * (35.0 - 945.0 * z_r_2 + 3465.0 * z_r_4 - 3003.0 * z_r_6) * x_r_1;
  aJ6[1] = (-J6 / 16.0 * (C_MU / pow(r, 2)) * pow((C_Req / r), 6)) * (35.0 - 945.0 * z_r_2 + 3465.0 * z_r_4 - 3003.0 * z_r_6) * y_r_1;
  aJ6[2] = (-J6 / 16.0 * (C_MU / pow(r, 2)) * pow((C_Req / r), 6)) * (245.0 - 2205.0 * z_r_2 + 4851.0 * z_r_4 - 3003.0 * z_r_6) * z_r_1;

  for (int i = 0; i <= 2; i++)
  {
    dX[i] = aTB[i] + aJ2[i] + aJ3[i] + aJ4[i] + aJ5[i] + aJ6[i];
  }
  return;
}

void Grav_Approx(double t, double *X, double *dX, double *Feval)
{
  // Return the J2-J6 Earth Gravity Evaluation with X in ECEF frame and increase func eval counter.
  Grav_EarthJ2J6(X, dX);
  Feval[1] = Feval[1] + 1.0;
  return;
}

void Grav_Full(double t, double *Xo, double *acc, double tol, double deg, double *Feval)
{

  double state[6] = {0.0};
  double dstate[6] = {0.0};

  for (int jj = 0; jj <= 5; jj++)
  {
    if (jj > 2)
    {
      state[jj] = 0.0;
    }
    else
    {
      state[jj] = Xo[jj];
    }
  }

  double grav = deg;
  // radial_gravity(Xo, tol, deg, &grav);
  EGM2008(state, &dstate[3], grav);
  Feval[0] = Feval[0] + pow(grav, 2) / pow(deg, 2);

  for (int jj = 0; jj <= 2; jj++)
  {
    acc[jj] = dstate[jj + 3];
  }

  return;
}

bool fullgravswitch(double err, double tol, int hot, int itr, IterCounters &ITRs)
{

  double Gapprox[3] = {0.0};
  // retrieve iteration counter values
  int &ITR1 = ITRs.ITR1;
  int &ITR2 = ITRs.ITR2;
  int &ITR3 = ITRs.ITR3;
  int &ITR4 = ITRs.ITR4;
  int &MODEL = ITRs.MODEL;
  bool grav_switch = false; // flag to switch to full gravity

  //////////////////////////////// J2-J6 ///////////////////////////////
  if (err > 1.0e-1)
  {
    // J2-J6
    grav_switch = false;
    if (debug_grav == 1)
    {
      printf("J2-J6\n");
    }
    ITR1 = itr;
    ITR2 = itr;
    ITR3 = itr;
    ITR4 = itr;
    if (debug_grav_itr == 1)
    {
      printf("ITR1 %i\n", ITR1);
    }
  }

  //////////////////////////////// 1e-1 1e-4 ///////////////////////////////
  else if (err <= 1.0e-1 && err > 1.0e-4 && ITR1 == itr - 1)
  //If the error is less than 1e-1 and greater than 1e-4 and the previous iteration was an approximate gravity iteration then switch to full gravity
  { // 1e-1 1e-4
    // FULL Gravity
    grav_switch = true;
    if (debug_grav == 1)
    {
      printf("Full Gravity 1\n");
    }
    MODEL = MODEL + 3;
  }
  else if (err <= 1.0e-1 && err > 1.0e-4)
  { // 1e-1 1e-4
  //If the error is less than 1e-1 and greater than 1e-4 and the first full gravity iteration has been performed then switch to approximate gravity
    // Approximate Gravity
    grav_switch = false;
    if (debug_grav == 1)
    {
      printf("Approx Gravity 1\n");
    }
    {
      ITR2 = itr;
      if (debug_grav_itr == 1)
      {
        printf("ITR2 %i\n", ITR2);
      }
    }
  }

  //////////////////////////////// 1e-4 1e-7 ///////////////////////////////
  // else if (err <= 1.0e-4 && err > 1.0e-7 && ITR2 == itr - 1){ // 1e-4 1e-7
  //   // FULL Gravity
  //   grav_switch = true;
  //   if (debug_grav == 1){
  //     printf("Full Gravity 2\n");

  //   }
  // }
  else if (err <= 1.0e-4 && err > 1.0e-7)
  { // 1e-4 1e-7
  //If the error is less than 1e-4 and greater than 1e-7 and the second full gravity iteration has been performed  (or skipped) then switch to approximate gravity
    // Approximate Gravity
    grav_switch = false;
    if (debug_grav == 1)
    {
      printf("Approx Gravity 2\n");
    }

    ITR3 = itr;
    if (debug_grav_itr == 1)
    {
      printf("ITR3 %i\n", ITR3);
    }
  }

  //////////////////////////////// 1e-7 1e-10 ///////////////////////////////
  else if (err <= 1.0e-7 && err > 1.0e-10 && ITR3 == itr - 1)
  //If the error is less than 1e-7 and greater than 1e-10 and the previous iteration was an approximate gravity iteration then switch to full gravity
  { // 1e-7 1e-10
    // FULL Gravity
    grav_switch = true;
    if (debug_grav == 1)
    {
    printf("Full Gravity 3\n");
    }
    MODEL = MODEL + 3;
  }

  else if (err <= 1.0e-7 && err > 1.0e-10)
  { // 1e-7 1e-10
  //If the error is less than 1e-7 and greater than 1e-10 and the third full gravity iteration has been performed then switch to approximate gravity
    // Approximate Gravity
    grav_switch = false;
    if (debug_grav == 1)
    {

      printf("Approx Gravity 3\n");
    }

    ITR4 = itr;
    if (debug_grav_itr == 1)
    {
      printf("ITR4 %i\n", ITR4);
    }
  }

  //////////////////////////////// 1e-10 1e-12 ///////////////////////////////
  // else if (err <= 1.0e-10 && err > 1.0e-12 && ITR4 == *itr - 1){ // 1e-10 1e-12
  else if (MODEL == 3 && err > tol)
  { // 1e-10 1e-12
  //If the error is less than 1e-10 and the previous iteration was an approximate gravity iteration then switch to full gravity
    // FULL Gravity
    grav_switch = true;
    if (debug_grav == 1)
    {

      printf("Full Gravity 4\n");
    }
    MODEL = MODEL + 3;
  }
  else if (err > tol)
  {
  //If the error is less than 1e-10 and the previous iteration was an approximate gravity iteration then switch to full gravity
    // Approximate Gravity
    grav_switch = false;
    if (debug_grav == 1)
    {

      printf("Approx Gravity 4\n");
    }
  }

  return grav_switch;
}

void variableGrav(double t, double *Xo, double *G, double *del_G, double tol, int i, double deg, double *Feval, bool fullgravswitch, string primary){
//Calculate the approximate or full gravity depending on the fullgravswitch flag
//set ode functions
  void (*ode_full)(double t, double *Xo, double *G, double tol, double deg, double *Feval);
  void (*ode_approx)(double t, double *Xo, double *G, double *Feval);


  if (primary=="EARTH"){
      ode_full=&Grav_Full;
      ode_approx=&Grav_Approx;

  }
  else if (primary=="MOON"){
    printf("function variableGrav: Primary Body not recognized\n");
      // ode_full=&lunar_Grav_Full;
      // ode_approx=&lunar_Grav_Approx;
  }
  else{
    printf("function variableGrav: Primary Body not recognized\n");
  }

  double Gapprox[3] = {0.0};
  if (fullgravswitch)
  {
    ode_full(t, Xo, G, tol, deg, Feval);
    ode_approx(t, Xo, Gapprox, Feval);
    for (int j = 0; j <= 2; j++)
    {
      del_G[ID2(i, j + 1, Nmax + 1)] = G[j] - Gapprox[j];
    }
  }
  else
  {
    ode_approx(t, Xo, Gapprox, Feval);
    for (int j = 0; j <= 2; j++)
    {
      G[j] = Gapprox[j] + del_G[ID2(i, j + 1, Nmax + 1)];
    }
    }
  return;
}