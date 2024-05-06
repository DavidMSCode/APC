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
#include <iostream>
#include <omp.h>
#include "picard_iteration.h"
#include "const.h"
#include "c_functions.h"
#include "FandG.h"
#include "eci2ecef.h"
#include "ecef2eci.h"
#include "lunar_perturbed_gravity.h"
#include "perturbed_gravity.h"
#include "picard_error_feedback.h"
#include "perturbations.h"
#include "Orbit.h"
#include "GRGM1200b.h"
#include "matrix_loader.h"
#include "Ephemeris.hpp"
// #include "EphemerisRotation.h"
#include "flags.h"
using namespace std;
void picard_iteration(double *Feval, Orbit &orbit, EphemerisManager &ephem)
{
  // Load constants
  int N = orbit.N;                         // Degree of Chebyshev polynomial
  int M = orbit.M;                         // Number of sample points
  int hot = orbit.prep_HS;                 // Hot start on/off switch condition
  double tol = orbit.tol;                  // Tolerance
  int deg = orbit.deg;                     // Gravity degree
  vector<double> &times = orbit.times_seg; // Time array for current segment
  vector<double> &P1 = orbit.P1;           // First integration operator (Acceleration to Velocity)
  vector<double> &P2 = orbit.P2;           // Second integration operator (Velocity to Position)
  vector<double> &T1 = orbit.T1;           // First Chebyshev matrix
  vector<double> &T2 = orbit.T2;           // Second Chebyshev matrix
  vector<double> &A = orbit.A;             // Least squares operator
  double *Xint = orbit.r0_seg;             // initial position of segment
  double *Vint = orbit.v0_seg;             // initial velocity of segment

  // Output solution variables
  vector<double> &Beta = orbit.Beta_seg;   // Velocity coefficients for current segment
  vector<double> &Alpha = orbit.Alpha_seg; // Position coefficients for current segment
  vector<double> &X = orbit.X_seg;         // Position warm start for current segment
  vector<double> &V = orbit.V_seg;         // Velocity warm start for current segment

  // Initialization
  double alt = 0.0;
  double xI[3] = {0.0};
  double vI[3] = {0.0};
  double xPrimaryFixed[3] = {0.0};
  double vPrimaryFixed[3] = {0.0};
  double aPrimaryFixed[3] = {0.0};
  double aI[3] = {0.0};
  double del_X[3] = {0.0};
  double del_aECI[3] = {0.0};
  double drag_aECEF[3] = {0.0};
  double SRP_aI[3] = {0.0};
  double third_body_aI[3] = {0.0};

  std::vector<double> G((M + 1) * 3, 0.0);
  std::vector<double> beta(N * 3, 0.0);
  std::vector<double> gamma(N * 3, 0.0);
  std::vector<double> alpha((N + 1) * 3, 0.0);
  std::vector<double> kappa((N + 1) * 3, 0.0);
  std::vector<double> Xorig;
  std::vector<double> Vorig;
  std::vector<double> Xnew; // temporary per loop storage for X
  std::vector<double> Vnew; // temporary per loop storage for V
  std::vector<double> xECEFp((M + 1) * 3, 0.0);
  std::vector<double> xECIp((M + 1) * 3, 0.0);
  std::vector<double> del_a((M + 1) * 3, 0.0);
  // Perturbed Gravity iteration storage
  IterCounters ITRs;
  // Initialization
  if (hot == 0)
  {
    ITRs.ITR1 = 0;
    ITRs.ITR2 = 0;
    ITRs.ITR3 = 0;
    ITRs.ITR4 = 0;
    ITRs.MODEL = 0;
  }

  // Initialization with hot start
  if (hot == 1)
  {
    ITRs.ITR1 = -1;
    ITRs.ITR2 = -1;
    ITRs.ITR3 = -1;
    ITRs.ITR4 = -1;
    ITRs.MODEL = 0;
  }

  vector<double> del_G(3 * (Nmax + 1), 0.0);
  int itr, MaxIt;
  double err, w2;
  itr = 0;
  MaxIt = 300;
  err = 10.0;
  w2 = (times[M] - times[0]) / 2.0;

  if (hot == 1)
  {
    err = 1e-2; // Prevents low fidelity J2 to J6 computations
  }

  while (err > tol) // Iterate over same segment until max error at any node (diff between current iteration and previous iteration) meets the tolerance
  {

    // vector<double> G_copy = G;
    // vector<double> del_G_copy = del_G;
    // double Feval_copy[2];
    // Feval_copy[0] = Feval[0];
    // Feval_copy[1] = Feval[1];
    // IterCounters ITRs_copy = ITRs;
    // int itr_copy = itr;
    // prefetch the current iteration's gravities using
    picardSegmentGravity(times, X, G, del_G, err, M, deg, orbit.lowDeg, hot, tol, itr, Feval, ITRs, orbit);

    for (int i = 1; i <= M + 1; i++) // Loop over each node and get forces in inertial frame
    {
      for (int j = 1; j <= 3; j++) // Get current nodes position and velocity
      {
        xI[j - 1] = X[ID2(i, j, M + 1)];
        vI[j - 1] = V[ID2(i, j, M + 1)];
      }
      for (int j = 1; j <= 3; j++)
      {
        // G[ID2(i, j, M + 1)] = aI[j - 1];
        xECIp[ID2(i, j, M + 1)] = xI[j - 1];
      }
    } // Forces found for each node


    // // Compare the output of G and G_copy
    // // If they are the same parallel gravity works
    // for (int i = 0; i <= 3 * (M + 1) - 1; i++)
    // {
    //   if (G[i] != G_copy[i])
    //   {
    //     string error_string = "Parallel gravity differs at node: " + to_string(i) + "\n";
    //     error_string += "iteration: " + to_string(itr) + "\n";
    //     error_string += " G: " + to_string(G[i]) + "\nG_copy: " + to_string(G_copy[i])+"\n";
    //     cout << error_string;
    //   }
    // }
    // //check if Del_G and Del_G_copy are the same
    // for (int i = 0; i <= 3 * (M + 1) - 1; i++)
    // {
    //   if (del_G[i] != del_G_copy[i])
    //   {
    //     string error_string = "Parallel del_G differs at node: " + to_string(i) + "\n";
    //     error_string += "iteration: " + to_string(itr) + "\n";
    //     error_string += " del_G: " + to_string(del_G[i]) + "\n del_G_copy: " + to_string(del_G_copy[i])+"\n";
    //     cout << error_string;
    //   }
    // }


    // // cehck if Feval[0] and [1] are the same as they're copy
    // if (Feval[0] != Feval_copy[0] || Feval[1] != Feval_copy[1])
    // {
    //   string error_string = "Parallel Feval differs\n";
    //   error_string += "iteration: " + to_string(itr) + "\n";
    //   error_string += " Feval[0]: " + to_string(Feval[0]) + " Feval_copy[0]: " + to_string(Feval_copy[0])+"\n";
    //   error_string += " Feval[1]: " + to_string(Feval[1]) + " Feval_copy[1]: " + to_string(Feval_copy[1])+"\n";
    //   cout << error_string;
    // }

  // //Check if the IterCounters are the same
  //   if(ITRs.ITR1!=ITRs_copy.ITR1 || ITRs.ITR2 != ITRs_copy.ITR2 || ITRs.ITR3 != ITRs_copy.ITR3 || ITRs.ITR4 != ITRs_copy.ITR4){
  //     string error_string = "Parallel ITRs differs\n";
  //     error_string += "iteration: " + to_string(itr) + "\n";
  //     error_string += " ITR1: " + to_string(ITRs.ITR1) + " ITRs_copy.ITR1: " + to_string(ITRs_copy.ITR1)+"\n";
  //     error_string += " ITR2: " + to_string(ITRs.ITR2) + " ITRs_copy.ITR2: " + to_string(ITRs_copy.ITR2)+"\n";
  //     error_string += " ITR3: " + to_string(ITRs.ITR3) + " ITRs_copy.ITR3: " + to_string(ITRs_copy.ITR3)+"\n";
  //     error_string += " ITR4: " + to_string(ITRs.ITR4) + " ITRs_copy.ITR4: " + to_string(ITRs_copy.ITR4)+"\n";
  //     cout << error_string;
  //   }

    // Perform quadrature for velocity then position in inertial frame
    // Velocity
    std::vector<double> tmp1;
    std::vector<double> tmp2;
    tmp1 = matmul(orbit.A, G, N - 1, M + 1, 3, N - 1, M + 1);
    tmp2 = matmul(orbit.P1, tmp1, N, N - 1, 3, N, N - 1);
    for (int i = 1; i <= N; i++)
    {
      for (int j = 1; j <= 3; j++)
      {
        beta[ID2(i, j, N)] = w2 * tmp2[ID2(i, j, N)];
        if (i == 1)
        {
          beta[ID2(i, j, N)] = beta[ID2(i, j, N)] + Vint[j - 1];
        }
      }
    }
    Vorig = matmul(orbit.T1, beta, M + 1, N, 3, M + 1, N);

    // Position
    std::vector<double> tmp3;
    tmp3 = matmul(orbit.P2, beta, N + 1, N, 3, N + 1, N);
    for (int i = 1; i <= N + 1; i++)
    {
      for (int j = 1; j <= 3; j++)
      {
        alpha[ID2(i, j, N + 1)] = w2 * tmp3[ID2(i, j, N + 1)];
        if (i == 1)
        {
          alpha[ID2(i, j, N + 1)] = alpha[ID2(i, j, N + 1)] + Xint[j - 1];
        }
      }
    }
    Xorig = matmul(orbit.T2, alpha, M + 1, N + 1, 3, M + 1, N + 1);

    for (int i = 1; i <= M + 1; i++)
    {
      for (int j = 1; j <= 3; j++)
      {
        xI[j - 1] = Xorig[ID2(i, j, M + 1)];
        vI[j - 1] = Vorig[ID2(i, j, M + 1)];
      }
      // Linear Error Correction Position
      for (int j = 1; j <= 3; j++)
      {
        del_X[j - 1] = xI[j - 1] - xECIp[ID2(i, j, M + 1)];
      }
      // Convert from inertial to body fixed frame
      // InertialToBodyFixed(xI,vI,xPrimaryFixed,vPrimaryFixed,times[i-1],orbit);
      eci2ecef(orbit.et(times[i - 1]), xI, vI, xPrimaryFixed, vPrimaryFixed);
      // Linear Error Correction Acceleration
      double del_aECEF[3] = {0.0};
      picard_error_feedback_GRGM1200b(xPrimaryFixed, del_X, del_aECEF);
      // Convert from ECEF to ECI
      ecef2eci(orbit.et(times[i - 1]), del_aECEF, del_aECI);
      // BodyFixedAccelerationToInertial(del_aECEF,del_aECI,times[i-1],orbit);

      for (int j = 1; j <= 3; j++)
      {
        del_a[ID2(i, j, M + 1)] = del_aECI[j - 1];
      }
    }

    // Linear Error Correction Velocity Coefficients
    std::vector<double> tmp4;
    std::vector<double> tmp5;
    tmp4 = matmul(orbit.A, del_a, N - 1, M + 1, 3, N - 1, M + 1);
    tmp5 = matmul(orbit.P1, tmp4, N, N - 1, 3, N, N - 1);
    for (int i = 1; i <= N; i++)
    {
      for (int j = 1; j <= 3; j++)
      {
        gamma[ID2(i, j, N)] = w2 * tmp5[ID2(i, j, N)];
      }
    }

    // Corrected Velocity
    for (int i = 1; i <= N; i++)
    {
      for (int j = 1; j <= 3; j++)
      {
        Beta[ID2(i, j, N)] = beta[ID2(i, j, N)];
        if (err < 1e-13)
        {
          Beta[ID2(i, j, N)] = beta[ID2(i, j, N)] + gamma[ID2(i, j, N)];
        }
      }
    }
    Vnew = matmul(orbit.T1, Beta, M + 1, N, 3, M + 1, N);

    // Corrected Position
    std::vector<double> tmp6;
    tmp6 = matmul(orbit.P2, gamma, N + 1, N, 3, N + 1, N);
    for (int i = 1; i <= N + 1; i++)
    {
      for (int j = 1; j <= 3; j++)
      {
        kappa[ID2(i, j, N + 1)] = w2 * tmp6[ID2(i, j, N + 1)];
        Alpha[ID2(i, j, N + 1)] = alpha[ID2(i, j, N + 1)];
        if (err < 1e-13)
        {
          Alpha[ID2(i, j, N + 1)] = alpha[ID2(i, j, N + 1)] + kappa[ID2(i, j, N + 1)];
        }
      }
    }
    Xnew = matmul(orbit.T2, Alpha, M + 1, N + 1, 3, M + 1, N + 1);

    // Non-dimensional Error
    double tmp = 0.0;
    double curr_err = 0.0;
    for (int i = 1; i <= M + 1; i++)
    {
      for (int j = 1; j <= 6; j++)
      {
        if (j <= 3)
        {
          tmp = fabs(Xnew[ID2(i, j, M + 1)] - X[ID2(i, j, M + 1)]) / DU;
        }
        if (j > 3)
        {
          tmp = fabs(Vnew[ID2(i, j - 3, M + 1)] - V[ID2(i, j - 3, M + 1)]) / DU * TU;
        }
        if (tmp > curr_err)
        {
          curr_err = tmp;
        }
      }
    }
    err = curr_err;

    // Update
    X = Xnew;
    V = Vnew;

    // Iteration Counter
    if (itr == MaxIt)
    {
      break;
      if (g_DEBUG_PICARD)
      {
        cout << "Warning: Max iterations reached for current segment" << endl;
      }
    }
    itr++;
  }

  if (g_DEBUG_PICARD)
  {
    cout << "Segment " << orbit.DebugData.segments.size() << " converged in " << itr << " iterations." << endl;
  }

  // DEBUG: Store segment data
  if (g_DEBUG_SEGMENTS)
  {
    struct segment Segment;

    // iterate through the X array
    for (int i = 1; i <= M + 1; i++)
    {
      //
      double state[6] = {X[ID2(i, 1, M + 1)], X[ID2(i, 2, M + 1)], X[ID2(i, 3, M + 1)], V[ID2(i, 1, M + 1)], V[ID2(i, 2, M + 1)], V[ID2(i, 3, M + 1)]};
      // Calculate hamiltonian
      double H;
      double et = orbit.et(times[i - 1]);
      double r[3] = {X[ID2(i, 1, M + 1)], X[ID2(i, 2, M + 1)], X[ID2(i, 3, M + 1)]};
      double v[3] = {V[ID2(i, 1, M + 1)], V[ID2(i, 2, M + 1)], V[ID2(i, 3, M + 1)]};
      jacobiIntegral_GRGM1200b(et, state, &H, orbit.deg, orbit);
      if (orbit.DebugData.segments.empty() && i == 1)
      {
        orbit.DebugData.H0 = H;
      }
      // Store Data
      Segment.x.push_back(state[0]);
      Segment.y.push_back(state[1]);
      Segment.z.push_back(state[2]);
      Segment.vx.push_back(state[3]);
      Segment.vy.push_back(state[4]);
      Segment.vz.push_back(state[5]);
      Segment.t.push_back(times[i - 1]);
      Segment.et.push_back(et);
      Segment.dH.push_back((H - orbit.DebugData.H0) / orbit.DebugData.H0);
    }
    orbit.DebugData.segments.push_back(Segment);
  }

  // calculate altcitude at each node
  for (int i = 1; i <= M + 1; i++)
  {
    double r[3] = {X[ID2(i, 1, M + 1)], X[ID2(i, 2, M + 1)], X[ID2(i, 3, M + 1)]};
    double R = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    double alt = R - orbit.GetPrimaryRadius();
    if (alt < 0.0)
    {
      orbit.suborbital = true;
      break;
    }
  }
  return;
}

void picardSegmentGravity(vector<double> times, vector<double> X, vector<double> &accel, vector<double> &del_G, double err, int M, int deg, int lowdeg, int hot, double tol, int itr, double *Feval, IterCounters &ITRs, Orbit &orbit)
{
  // determine which gravity model to use
  bool fullGrav = fullgravswitch(err, tol, hot, itr, ITRs);


  // Calculate the acceleration at each node of the segment in parallel. F
  #pragma omp parallel for reduction(+:Feval[0] , Feval[1])
  for (int i = 1; i <= M + 1; i++)
  {
    double xI[3] = {0.0};
    double vI[3] = {0.0};
    double xECEF[3] = {0.0};
    double vECEF[3] = {0.0};
    double aI[3] = {0.0};
    double aECEF[3] = {0.0};
    for (int j = 1; j <= 3; j++) // Get current nodes position and velocity
    {
      xI[j - 1] = X[ID2(i, j, M + 1)];
    }
    eci2ecef(orbit.et(times[i - 1]), xI, vI, xECEF, vECEF);

    variableGrav(times[i - 1], xECEF, aECEF, &del_G[0], tol, i, deg, Feval, fullGrav, orbit.GetPrimaryBody());
    // lunar_perturbed_gravity(times[i - 1], xECEF, err, i, M, deg, hot, aECEF, tol, &itr, Feval, ITRs, &del_G[0], lowdeg);

    ecef2eci(orbit.et(times[i - 1]), aECEF, aI);

    for (int j = 1; j <= 3; j++)
    {
      accel[ID2(i, j, M + 1)] = aI[j - 1];
    }
  }
  return;
}