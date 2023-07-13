/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     May 2017
*  LAST MODIFIED:    May 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Prepares arrays for storing iteration computations
*
* INPUT:
*    r0      -- Initial position (km)
*    v0      -- Initial velocity (km/s)
*    t0      -- Initial time (s)
*    t_final -- Final time (s)
*    dt      -- Time interval (s)
*    tp      -- Approximate time of perigee passage (s)
*    tol     -- Tolerance
*    N       -- Polynomial degree
*    M       -- Sample points
*    seg     -- Segmentss per orbit
*    prep_HS -- Hot start switch condition
*
* OUTPUTS:
*    tvec    -- Segment start and end times (s)
*    t_orig  -- Segment start and end times for first segment (s)
*    P1      -- First integration operator
*    P2      -- Second integration operator
*    T1      -- Chebyshev velocity matrix
*    T2      -- Chebyshev position matrix
*    A       -- Least squares operator
*    Ta      -- Chebyshev acceleration matrix
*
* COMMENTS:
*
*/
#include <math.h>

#include <vector>
#include <iostream>

#include <omp.h>

#include "prepare_propagator.h"
#include "const.h"
#include "matrix_loader.h"
#include "rv2elm.h"
#include "c_functions.h"
#include "flags.h"
#include "Orbit.h"

using namespace std;


void prepare_propagator(double tol, int* prep_HS, Orbit &orbit){
  //Initialize arrays
  int N  = orbit.N;
  int M = orbit.N;                // # sample points = polynomial degree
  vector<double> T2((M+1)*(N+1),0.0);   // [(M+1)x(N+1)]
  vector<double> P2((N+1)*N,0.0);       // [(N+1)xN]
  vector<double> T1((M+1)*N,0.0);       // [(M+1)xN]
  vector<double> P1(N*(N-1),0.0);       // [Nx(N-1)]
  vector<double> Ta((M+1)*(N-1),0.0);   // [(M+1)x(N-1)]
  vector<double> A((N-1)*(M+1),0.0);    // [(N-1)x(M+1)]
  vector<double> t_orig(orbit.seg+1,0.0);
  vector<double> tvec(orbit.seg+1,0.0);

  //Get values from orbit class
  int seg = orbit.seg;
  double* r0 = &orbit._In_r0[0];
  double* v0 = &orbit._In_v0[0];
  double t0 = orbit._Integrator_t0;
  double t_final = orbit._Integrator_tf;
  double dt = orbit._dt;
  double tp = orbit.tp;
  // Compute Keplerian Orbit Period
  double a, e, Period, n;
  double elm[10] = {0.0};
  double mu = orbit.GetPrimaryGravitationalParameter();
  rv2elm(r0,v0,tol,mu,elm);
  a      = elm[1];
  e      = elm[2];
  Period = 2.0*C_PI*sqrt(pow(a,3)/mu);
  n      = 2.0*C_PI/Period;

  // Compute Time Vector (based on Keplerian true anomaly segments)
  double f, E, MA, df;
  t_orig[0] = 0.0;
  tvec[0]   = 0.0;
  df        = 2.0*C_PI/seg;
  f         = 0.0;
  E         = 0.0;
  MA        = 0.0;
  for (int i=1; i<=seg; i++){
    f = f + df;
    E = 2.0*atan2(tan(0.5*f)*sqrt(1.0-e),sqrt(1.0+e));
    if (E < 0){
      E = 2.0*C_PI + E;
    }
    MA        = E - e*sin(E);
    t_orig[i] = MA/n;
    tvec[i]   = MA/n;
  }

  // Short first segment if user specified IC's do not coincide with segment break
  double ts;
  ts = Period - tp;

  if (fabs(tp) > 1.0e-5){
    for (int i=0; i<=seg; i++){
      if (t_orig[i] < ts){
        tvec[i] = 0.0;
      }
      if (t_orig[i] >= ts){
        tvec[i] = t_orig[i]-ts;
      }
    }
    *prep_HS = 0;
  }

  // User specified time vector for output
  int len;
  len = ceil(t_final/dt);
  std::vector<double> time_out(len,0.0);
  //memset( time_out, 0.0, (len*sizeof(double)));
  time_out[0] = t0;
  for (int i=1; i<len; i++){
    time_out[i] = time_out[i-1] + dt;
  }

  // LOAD PRECOMPUTED MATRICES
    double* temp1;
    double* temp2;
    double* temp3;
    double* temp4;
    double* temp5;
    double* temp6;

  #pragma omp critical(matrixloader)    //ensure only one thread loads the Picard iteration matrices
  {
    if (!g_MATRICES_LOADED)
    {
      //matrices aren't loaded so load them now
      matrix_loader();
      
    }
  }
    int idN;
    idN = (N-10);
    // Retrive Data from Storage Arrays

    temp1 = &arr_T2[idN][0];
    temp2 = &arr_P2[idN][0];
    temp3 = &arr_T1[idN][0];
    temp4 = &arr_P1[idN][0];
    temp5 = &arr_Ta[idN][0];
    temp6 = &arr_A[idN][0];
  
  // BUILD MATRICES
  for (int j=1; j<=M+1; j++){
    for (int k=1; k<=N+1; k++){
      T2[ID2(j,k,M+1)] = temp1[ID2(j,k,Nmax+1)];  // Chebyshev Position Matrix
    }
  }
  for (int j=1; j<=N+1; j++){
    for (int k=1; k<=N; k++){
      P2[ID2(j,k,N+1)] = temp2[ID2(j,k,Nmax+1)];  // Second Integration Operator (Velocity to Position)
    }
  }
  for (int j=1; j<=M+1; j++){
    for (int k=1; k<=N; k++){
      T1[ID2(j,k,M+1)] = temp3[ID2(j,k,Nmax+1)];  // Chebyshev Velocity Matrix
    }
  }
  for (int j=1; j<=N; j++){
    for (int k=1; k<=N-1; k++){
      P1[ID2(j,k,N)] = temp4[ID2(j,k,Nmax+1)];    // First Integration Operator (Acceleration to Velocity)
    }
  }
  for (int j=1; j<=M+1; j++){
    for (int k=1; k<=N-1; k++){
      Ta[ID2(j,k,M+1)] = temp5[ID2(j,k,Nmax+1)];  // Chebyshev Acceleration Matrix
    }
  }
  for (int j=1; j<=N-1; j++){
    for (int k=1; k<=M+1; k++){
      A[ID2(j,k,N-1)] = temp6[ID2(j,k,Nmax+1)];   // Least Squares Operator
    }
  }
//std::cout << "finished loading matrices";
orbit.T2 = T2;   // [(M+1)x(N+1)]
orbit.P2 = P2;       // [(N+1)xN]
orbit.T1 = T1;       // [(M+1)xN]
orbit.P1 = P1;       // [Nx(N-1)]
orbit.Ta = Ta;   // [(M+1)x(N-1)]
orbit.A = A;    // [(N-1)x(M+1)]
orbit.t_orig = t_orig;
orbit.tvec = tvec;
}
