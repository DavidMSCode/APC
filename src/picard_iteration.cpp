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
*    T1      -- First Chebyshev matrix
*    A       -- Least squares operator
*    Feval   -- Function evaluation counter
*
* OUTPUTS:
*    X       -- Position solution for current segment (km)
*    V       -- Velocity solution for current segment (km/s)
*    Alpha   -- Position coefficients for current segment
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
#include "perturbed_gravity.h"
#include "picard_error_feedback.h"
#include "perturbations.h"
#include "Orbit.h"
#include "EGM2008.h"
#include "matrix_loader.h"
#include "Ephemeris.hpp"
#include "mee_dot.h"
#include "mee2rv.h"

void picard_iteration(std::vector<double> &X, std::vector<double> &V, std::vector<double> &MEE, mee0, std::vector<double> &times, int N, int M, double deg, int hot, double tol, std::vector<double> &P1, std::vector<double> &T1, std::vector<double> &A, double* Feval, std::vector<double> &Alpha, Orbit &orb, EphemerisManager ephem){
  
  // Initialization
  bool suborbital = false;
  double alt = 0.0;
  double xI[3]        = {0.0};
  double vI[3]        = {0.0};
  double xECEF[3]     = {0.0};
  double vECEF[3]     = {0.0};
  double aECEF[3]     = {0.0};
  double aECI[3]      = {0.0};
  double drag_aECEF[3] = {0.0};
  double SRP_aECI[3] = {0.0};
  double third_body_aECI[3] = {0.0};

  std::vector<double> MEEnew((M+1)*6,0.0);
  std::vector<double> MEE_dot((M+1)*6,0.0);
  std::vector<double> a_lvlh((M+1)*3,0.0);

  double TB[3] = {0.0};
  double ad_eci[3] = {0.0}; 
  double xrdl[3] = {0.0};
  double yrdl[3] = {0.0};
  double zrdl[3] = {0.0};
  double r_eci[3] = {0.0};
  double v_eci[3] = {0.0};
  double mee[6] = {0.0};
    
  //Perturbed Gravity iteration storage
  IterCounters ITRs;

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
        alt = sqrt(pow(xI[0],2)+pow(xI[1],2)+pow(xI[2],2))-C_Req;
        if (alt<0){
          suborbital=true;
        }
      }
      // Convert from ECI to ECEF
      eci2ecef(times[i-1],xI,vI,xECEF,vECEF);
      // Compute Variable Fidelity Gravity
      perturbed_gravity(times[i-1],xECEF,err,i,M,deg,hot,aECEF,tol,&itr,Feval,ITRs);
      //Calculate acceleration from drag
      Perturbed_Drag(xECEF, vECEF, orb, drag_aECEF);

      //sum pertubed gravity and drag accelerations
      for(int k=0;k<3;k++){
        aECEF[k] = aECEF[k]+drag_aECEF[k];
      }

      // Convert from ECEF to ECI
      ecef2eci(times[i-1],aECEF,aECI);
        
      // Compute Two-Body Acceleration
      for (int j=0; j<=2; j++){
          TB[j] = -C_MU*xI[j]/pow(sqrt(pow(xI[0],2)+pow(xI[1],2)+pow(xI[2],2)),3);
      }
      // Compute Perturbed Acceleration
      for (int j=0; j<=2; j++){
          ad_eci[j] = aECI[j] - TB[j];
      }
        
      //calculate SRP and Third Body
      Perturbed_SRP(times[i-1], xI, orb, ephem, SRP_aECI);
      Perturbed_three_body(times[i-1], xI, orb, ephem, third_body_aECI);
      //Add perturbations to acceleration.
      for(int k=0;k<3;k++){
        ad_eci[k] = ad_eci[k] + SRP_aECI[k] + third_body_aECI[k];
      }
        
      // Inertial to Radial
      inertial2radial(xI,vI,xrdl,yrdl,zrdl);
        
      // Acceleration Components
      a_lvlh[ID2(i,1,M+1)] = ad_eci[0]*xrdl[0]+ad_eci[1]*xrdl[1]+ad_eci[2]*xrdl[2];
      a_lvlh[ID2(i,2,M+1)] = ad_eci[0]*yrdl[0]+ad_eci[1]*yrdl[1]+ad_eci[2]*yrdl[2];
      a_lvlh[ID2(i,3,M+1)] = ad_eci[0]*zrdl[0]+ad_eci[1]*zrdl[1]+ad_eci[2]*zrdl[2];
        
    }
      
    // MEE Dynamics
    mee_dot(MEE,a_lvlh,MEE_dot,M);
    
    // Integrate
    std::vector<double> tmp1;
    tmp1 = matmul(A,MEE_dot,N+1,N,6,N+1,N);   // LSQ Coefficients
    Alpha = matmul(P1,tmp1,N+1,N,6,N+1,N);    // Integration Operator
    for (int i=1; i<=M+1; i++){
        for (int j=1; j<=6; j++){
          Alpha[ID2(i,j,M+1)] = w2*Alpha[ID2(i,j,M+1)];
            if (i==1){
                Alpha[ID2(i,j,M+1)] = Alpha[ID2(i,j,M+1)] + mee0[j-1]; 
            }
        }
    }

    MEEnew = matmul(T1,Alpha,M+1,N+1,6,M+1,N+1);

    // Non-dimensional Error
    double tmp = 0.0;
    double curr_err = 0.0;
    for (int i=1; i<=M+1; i++){
      for (int j=1; j<=6; j++){
        if (j==1){
          // Nondimensionalize p
          tmp = fabs(MEEnew[ID2(i,j,M+1)] - MEE[ID2(i,j,M+1)])/DU;
        }
        if (j>1){
          tmp = fabs(MEEnew[ID2(i,j,M+1)] - MEE[ID2(i,j,M+1)]);
        }
        if (tmp > curr_err){
          curr_err = tmp;
        }
      }
    }
    err = curr_err;

    // Update
    MEE = MEEnew;
      
    // MEE to RV
   for (int i=1; i<=M+1; i++){
       for (int j=1; j<=6; j++){
           mee[j-1] = MEE[ID2(i,j,M+1)];
       }
       mee2rv(mee,r_eci,v_eci);
       for (int j=0; j<=2; j++){
           X[ID2(i,j+1,M+1)] = r_eci[j];
           V[ID2(i,j+1,M+1)] = v_eci[j];
       }
   }

    // Iteration Counter
    if (itr == MaxIt){
      itr = itr - 1;
      break;
    }
  }

  if(suborbital){
    // Set suborbital flag to stop iterations early.
    orb.SetSubOrbital();
    } 
}
