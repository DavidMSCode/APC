/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     May 2017
*  LAST MODIFIED:    May 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Picard Iteration
*
* INPUT:
*    Xinit   -- Initial position (km)
*    Vinit   -- Initial velocity (km/s)
*    X       -- Position warm start for current segment (km)
*    V       -- Velocity warm start for current segment (km/s)
*    times   -- Time array for current segment (s)
*    N       -- Degree of Chebyshev polynomial
*    M       -- Number of sample points
*    hot     -- Hot start on/off switch condition
*    tol     -- Tolerance
*    P1      -- First integration operator (Acceleration to Velocity)
*    P2      -- Second integration operator (Velocity to Position)
*    T1      -- First Chebyshev matrix
*    T2      -- Second Chebyshev matrix
*    A       -- Least squares operator
*    Feval   -- Function evaluation counter
*
* OUTPUTS:
*    X       -- Position solution for current segment (km)
*    V       -- Velocity solution for current segment (km/s)
*    Alpha   -- Position coefficients for current segment
*    Beta    -- Velocity coefficients for current segment
*
* REFERENCES:
* 1. Junkins, J.L., and Woollands, R., "Nonlinear Differential Equations Solvers via Adaptive Picard-Chebyshev Iteration: Applications in Astrodynamics",
*    AAS/AIAA Astrodynamics Specialist Conference, Stevenson, WA, 2017.
* 2. Junkins, J.L., and Woollands, R., "Adaptive-Picard-Chebyshev for Propagating Perturbed Two-Body Orbits",
*    JGCD, submitted 2017.
*/
#include <math.h>

#include <vector>

#include "picard_iteration.h"
#include "const.h"
#include "c_functions.h"
#include "FandG.h"
#include "eci2ecef.h"
#include "ecef2eci.h"
#include "lunar_perturbed_gravity.h"
#include "picard_error_feedback.h"
#include "perturbations.h"
#include "Orbit.h"
#include "EGM2008.h"
#include "matrix_loader.h"
#include "Ephemeris.hpp"


void picard_iteration(double* Xint, double* Vint, std::vector<double> &X, std::vector<double> &V, std::vector<double> &times, int N, int M, double deg, int hot, double tol,
  std::vector<double> &P1, std::vector<double> &P2, std::vector<double> &T1, std::vector<double> &T2, std::vector<double> &A, double* Feval, std::vector<double> &Alpha, std::vector<double> &Beta, Orbit &orbit, EphemerisManager &ephem){
  
  // Initialization
  bool suborbital = false;
  double alt = 0.0;
  double xI[3]        = {0.0};
  double vI[3]        = {0.0};
  double xECEF[3]     = {0.0};
  double vECEF[3]     = {0.0};
  double aECEF[3]     = {0.0};
  double aECI[3]      = {0.0};
  double del_X[3]     = {0.0};
  double del_aECEF[3] = {0.0};
  double del_aECI[3]  = {0.0};
  double drag_aECEF[3] = {0.0};
  double SRP_aECI[3] = {0.0};
  double third_body_aECI[3] = {0.0};

  std::vector<double> G((M+1)*3,0.0);
  std::vector<double> beta(N*3,0.0);
  std::vector<double> gamma(N*3,0.0);
  std::vector<double> alpha((N+1)*3,0.0);
  std::vector<double> kappa((N+1)*3,0.0);
  std::vector<double> Xorig;
  std::vector<double> Vorig;
  std::vector<double> Xnew;
  std::vector<double> Vnew;
  std::vector<double> xECEFp((M+1)*3,0.0);
  std::vector<double> xECIp((M+1)*3,0.0);
  std::vector<double> del_a((M+1)*3,0.0);
  //Perturbed Gravity iteration storage
  IterCounters ITRs;
  double del_G[3*(Nmax+1)];

  int itr, MaxIt;
  double err, w2;
  itr   = 0;
  MaxIt = 30;
  err   = 10.0;
  w2    = (times[M]-times[0])/2.0;

  if (hot == 1){
    err = 1e-2; // Prevents low fidelity J2 to J6 computations
  }

  while(err > tol){
    suborbital = false;
    for (int i=1; i<=M+1; i++){

      for (int j=1; j<=3; j++){
        xI[j-1] = X[ID2(i,j,M+1)];
        vI[j-1] = V[ID2(i,j,M+1)];
      }
      // Exit loop early if orbit has crashed
      if (!suborbital){
        //FIXME: Replace C_Req with primary radius that should be stored in the orbit file
        alt = sqrt(pow(xI[0],2)+pow(xI[1],2)+pow(xI[2],2))-C_Req;
        if (alt<0){
          suborbital=true;
        }
      }
      // Convert from ECI to ECEF
      eci2ecef(times[i-1],xI,vI,xECEF,vECEF);
      // Compute Variable Fidelity Gravity
      lunar_perturbed_gravity(times[i-1],xECEF,err,i,M,deg,hot,aECEF,tol,&itr,Feval,ITRs,del_G);
      //Calculate acceleration from drag
      Perturbed_Drag(xECEF, vECEF, orbit, drag_aECEF);

      //sum pertubed gravity and drag accelerations
      for(int k=0;k<3;k++){
        aECEF[k] = aECEF[k]+drag_aECEF[k];
      }

      // Convert from ECEF to ECI
      ecef2eci(times[i-1],aECEF,aECI);
      //calculate SRP and Third Body
      Perturbed_SRP(times[i-1], xI, orbit, ephem, SRP_aECI);
      Perturbed_three_body(times[i-1], xI, orbit, ephem, third_body_aECI);
      //Add perturbations to acceleration.
      for(int k=0;k<3;k++){
        aECI[k] = aECI[k] + SRP_aECI[k] + third_body_aECI[k];
      }
      for (int j=1; j<=3; j++){
        G[ID2(i,j,M+1)]      = aECI[j-1];
        xECIp[ID2(i,j,M+1)] = xI[j-1];
      }
    }
    
    // Velocity
    std::vector<double> tmp1;
    std::vector<double> tmp2;
    tmp1 = matmul(A,G,N-1,M+1,3,N-1,M+1);
    tmp2 = matmul(P1,tmp1,N,N-1,3,N,N-1);
    for (int i=1; i<=N; i++){
      for (int j=1; j<=3; j++){
        beta[ID2(i,j,N)] = w2*tmp2[ID2(i,j,N)];
        if (i == 1){
          beta[ID2(i,j,N)] = beta[ID2(i,j,N)] + Vint[j-1];
        }
      }
    }
    Vorig = matmul(T1,beta,M+1,N,3,M+1,N);

    // Position
    std::vector<double> tmp3;
    tmp3 = matmul(P2,beta,N+1,N,3,N+1,N);
    for (int i=1; i<=N+1; i++){
      for (int j=1; j<=3; j++){
        alpha[ID2(i,j,N+1)] = w2*tmp3[ID2(i,j,N+1)];
        if (i == 1){
          alpha[ID2(i,j,N+1)] = alpha[ID2(i,j,N+1)] + Xint[j-1];
        }
      }
    }
    Xorig = matmul(T2,alpha,M+1,N+1,3,M+1,N+1);

    for (int i=1; i<=M+1; i++){
      for (int j=1; j<=3; j++){
        xI[j-1] = Xorig[ID2(i,j,M+1)];
        vI[j-1] = Vorig[ID2(i,j,M+1)];
      }
      // Linear Error Correction Position
      for (int j=1; j<=3; j++){
        del_X[j-1] = xI[j-1] - xECIp[ID2(i,j,M+1)];
      }
      // Convert from ECI to ECEF
      eci2ecef(times[i-1],xI,vI,xECEF,vECEF);
      // Linear Error Correction Acceleration
      picard_error_feedback(xECEF,del_X,del_aECEF);
      // Convert from ECEF to ECI
      ecef2eci(times[i-1],del_aECEF,del_aECI);
      for (int j=1; j<=3; j++){
        del_a[ID2(i,j,M+1)] = del_aECI[j-1];
      }
    }

    // Linear Error Correction Velocity Coefficients
    std::vector<double> tmp4;
    std::vector<double> tmp5;
    tmp4 = matmul(A,del_a,N-1,M+1,3,N-1,M+1);
    tmp5 = matmul(P1,tmp4,N,N-1,3,N,N-1);
    for (int i=1; i<=N; i++){
      for (int j=1; j<=3; j++){
        gamma[ID2(i,j,N)] = w2*tmp5[ID2(i,j,N)];
      }
    }

    // Corrected Velocity
    for (int i=1; i<=N; i++){
      for (int j=1; j<=3; j++){
        Beta[ID2(i,j,N)] = beta[ID2(i,j,N)];
        if (err < 1e-13){
          Beta[ID2(i,j,N)] = beta[ID2(i,j,N)] + gamma[ID2(i,j,N)];
        }
      }
    }
    Vnew = matmul(T1,Beta,M+1,N,3,M+1,N);

    // Corrected Position
    std::vector<double> tmp6;
    tmp6 = matmul(P2,gamma,N+1,N,3,N+1,N);
    for (int i=1; i<=N+1; i++){
      for (int j=1; j<=3; j++){
        kappa[ID2(i,j,N+1)] = w2*tmp6[ID2(i,j,N+1)];
        Alpha[ID2(i,j,N+1)] = alpha[ID2(i,j,N+1)];
        if (err < 1e-13){
          Alpha[ID2(i,j,N+1)] = alpha[ID2(i,j,N+1)] + kappa[ID2(i,j,N+1)];
        }
      }
    }
    Xnew = matmul(T2,Alpha,M+1,N+1,3,M+1,N+1);

    // Non-dimensional Error
    double tmp = 0.0;
    double curr_err = 0.0;
    for (int i=1; i<=M+1; i++){
      for (int j=1; j<=6; j++){
        if (j<=3){
          tmp = fabs(Xnew[ID2(i,j,M+1)] - X[ID2(i,j,M+1)])/DU;
        }
        if (j>3){
          tmp = fabs(Vnew[ID2(i,j-3,M+1)] - V[ID2(i,j-3,M+1)])/DU*TU;
        }
        if (tmp > curr_err){
          curr_err = tmp;
        }
      }
    }
    err = curr_err;

    // Update
    X = Xnew;
    V = Vnew;

    // Iteration Counter
    if (itr == MaxIt){
      itr = itr - 1;
      break;
    }
  }
  
  

  

  if(suborbital){
    // Set suborbital flag to stop iterations early.
    orbit.SetSubOrbital();
    } 
}
