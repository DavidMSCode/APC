/*
 *  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
 *  DATE WRITTEN:     May 2017
 *  LAST MODIFIED:    May 2017
 *  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
 *  DESCRIPTION:      Computes the number of segments per orbit and the degree of the polynomial required to fit acceleration to the user specified accuracy
 *
 * INPUT:
 *    r0     -- Initial position vector (km)
 *    v0     -- Initial velocity vector (km/s)
 *    deg    -- Gravity Degree (max 100)
 *    tol    -- Tolerance
 *    Feval  -- Function evaluation counter
 *
 * OUTPUTS:
 *    seg    -- Segments per orbit
 *    degree -- Degree of Chebyshev Polynomial
 *    tp     -- Approximate time of Keplerian perigee passage (s)
 *    Period -- Keplerian orbit period (s)
 *
 * REFERENCES:
 * 1. Junkins, J.L., and Woollands, R., "Nonlinear Differential Equations Solvers via Adaptive Picard-Chebyshev Iteration: Applications in Astrodynamics",
 *    AAS/AIAA Astrodynamics Specialist Conference, Stevenson, WA, 2017.
 * 2. Junkins, J.L., and Woollands, R., "Adaptive-Picard-Chebyshev for Propagating Perturbed Two-Body Orbits",
 *    JGCD, submitted 2017.
 *
 * COMMENTS:
 *
 */

#include <math.h>
#include <vector>
#include <cstring>

#include "polydegree_segments.h"
#include "const.h"
#include "rv2elm.h"
#include "FandG.h"
#include "c_functions.h"
#include "perigee_approx.h"
#include "radial_gravity.h"
#include "Gravity.h"
#include "lsq_chebyshev_fit.h"
#include "matrix_loader.h"
#include "Orbit.h"
using namespace std;
// FIXME: remove redundant input parameters, pass and store orbital properties with orbit class
void polydegree_segments(Orbit &orbit, double *Feval)
{
  double *r0 = &orbit._In_r0[0]; // get input states as pointers to start of vector (acts like array)
  double *v0 = &orbit._In_v0[0];
  double &Period = orbit.Period;
  int &deg = orbit.deg;
  double tol = orbit.tol;
  int &seg = orbit.seg;
  int &degree = orbit.N;
  double &tp = orbit.tp;
  // Get gravitational parameter for primary body
  const double mu = orbit.GetPrimaryGravitationalParameter();
  // Tolerances
  double coeff_tol = tol / 1000.0;
  double fit_tol = tol / 100.0;

  // Compute Keplerian Orbit Period
  double a, e;
  double elm[10] = {0.0};
  rv2elm(r0, v0, tol, mu, elm);
  a = elm[1];
  e = elm[2];
  Period = 2.0 * C_PI * sqrt(pow(a, 3) / mu);

  // Compute F&G Analytical Solution for 1 Orbit
  int n = 100;
  double w1, w2, z0[6], z[6];
  double tau[300];
  memset(tau, 0.0, (300 * sizeof(double)));
  double times[300];
  memset(times, 0.0, (300 * sizeof(double)));
  double X[300 * 6];
  memset(X, 0.0, ((300 * 6) * sizeof(double)));
  double normX[300];
  memset(normX, 0.0, (300 * sizeof(double)));
  w1 = Period / 2.0;
  w2 = Period / 2.0;
  z0[0] = r0[0];
  z0[1] = r0[1];
  z0[2] = r0[2];
  z0[3] = v0[0];
  z0[4] = v0[1];
  z0[5] = v0[2];
  for (int cnt = 0; cnt <= n; cnt++)
  {
    tau[cnt] = -cos(cnt * C_PI / n);
    times[cnt] = tau[cnt] * w2 + w1;
    FandG(z0, z, times[cnt], mu);
    for (int i = 0; i <= 5; i++)
    {
      X[ID2(cnt + 1, i + 1, n + 1)] = z[i];
    }
    normX[cnt] = sqrt(pow(z[0], 2) + pow(z[1], 2) + pow(z[2], 2));
  }

  // Determine Approximate State and Time at Keplerian Perigee
  double rp[3], vp[3];
  int ind;
  perigee_approx(X, normX, times, n, rp, vp, &tp, &ind);

  /* If user specified ICs were not at perigee (i.e. ind!=0) then we need a
    finer grid to find perigee. I put 100 sample points in the vacinity of
    the perigee point and search again for a more accurate value. */
  if (ind != 0)
  {
    double diff;
    double tnew[300];
    memset(tnew, 0.0, (300 * sizeof(double)));
    diff = (times[ind + 1] - times[ind - 1]) / 99.0;
    tnew[0] = times[ind - 1];
    for (int cnt = 1; cnt <= n; cnt++)
    {
      tnew[cnt] = tnew[cnt - 1] + diff;
      FandG(z0, z, tnew[cnt], mu);
      for (int i = 0; i <= 5; i++)
      {
        X[ID2(cnt + 1, i + 1, n + 1)] = z[i];
      }
      normX[cnt] = sqrt(pow(z[0], 2) + pow(z[1], 2) + pow(z[2], 2));
    }
    perigee_approx(X, normX, tnew, n, rp, vp, &tp, &ind);
  }

  // Prepare for loop
  int fit_check, coeff, N, Nprev, jmax;
  fit_check = 0; // Loop check condition
  seg = 3;       // Minimum number of segments per orbit
  coeff = 3;     // Value of last 3 coefficients must be below tolerance
  jmax = 2;      // Loop Maximum (corresponds to Nmax = 80)

  // Perigee initial conditions
  z0[0] = rp[0];
  z0[1] = rp[1];
  z0[2] = rp[2];
  z0[3] = vp[0];
  z0[4] = vp[1];
  z0[5] = vp[2];

  double Nvec[4] = {0.0};
  Nvec[0] = 10;
  for (int i = 1; i <= jmax; i++)
  {
    Nvec[i] = 2 * Nvec[i - 1];
  }

  double f, E, Mf, MM, t0, tf1;

  // Loop through different combinations of segments & polynomial degree
  while (fit_check <= coeff)
  {

    f = 2.0 * C_PI / seg;                                     // True Anomaly
    E = 2.0 * atan2(tan(0.5 * f) * sqrt(1 - e), sqrt(1 + e)); // Eccentric anomaly
    if (E < 0.0)
    {
      E = 2.0 * C_PI + E;
    }
    Mf = E - e * sin(E);      // Mean anomaly
    MM = 2.0 * C_PI / Period; // Mean motion
    t0 = 0.0;                 // Initial Time
    tf1 = Mf / MM;            // Final Time (1 segment)
    // Scaling parameters
    w1 = (tf1 + t0) / 2.0;
    w2 = (tf1 - t0) / 2.0;

    // Loop through different values for N
    std::vector<double> G(300 * 3, 0.0);
    std::vector<double> Gprev(300 * 3, 0.0);
    for (int j = 0; j <= jmax; j++)
    {
      N = Nvec[j];

      // Generate discrete cosine nodes
      if (j == 0)
      {
        for (int cnt = 0; cnt <= N; cnt++)
        {
          tau[cnt] = -cos(cnt * C_PI / N);
          times[cnt] = tau[cnt] * w2 + w1;
          FandG(z0, z, times[cnt], mu);
          // FIXME: Use radial gravity for moon once radial gravity supports it
          // use radial gravity only if around the Earth for now.
          double grav = 0.0;
          if (orbit.InSet(orbit.GetPrimaryBody(), {"Earth"}))
          { // Check if primary body is Earth
            radial_gravity(z, tol, deg, &grav);
          }
          else
          {
            grav = deg;
          }
          double dstate[6] = {0.0};
          Gravity(z, &dstate[3], grav, orbit);
          Feval[0] = Feval[0] + pow(grav, 2) / pow(deg, 2);
          Gprev[ID2(cnt + 1, 1, N + 1)] = dstate[3];
          Gprev[ID2(cnt + 1, 2, N + 1)] = dstate[4];
          Gprev[ID2(cnt + 1, 3, N + 1)] = dstate[5];
        }
      }
      if (j > 0)
      {
        for (int cnt = 0; cnt <= N; cnt++)
        {
          if (cnt % 2 == 0)
          {
            G[ID2(cnt + 1, 1, N + 1)] = Gprev[ID2(cnt / 2 + 1, 1, Nprev + 1)];
            G[ID2(cnt + 1, 2, N + 1)] = Gprev[ID2(cnt / 2 + 1, 2, Nprev + 1)];
            G[ID2(cnt + 1, 3, N + 1)] = Gprev[ID2(cnt / 2 + 1, 3, Nprev + 1)];
          }
          if (cnt % 2 == 1)
          {
            tau[(cnt - 1) / 2] = -cos(cnt * C_PI / N);
            times[(cnt - 1) / 2] = tau[(cnt - 1) / 2] * w2 + w1;
            FandG(z0, z, times[(cnt - 1) / 2], mu);
            // FIXME: Use radial gravity for moon once radial gravity supports it
            // use radial gravity only if around the Earth for now.
            double grav = 0.0;
            if (orbit.InSet(orbit.GetPrimaryBody(), {"Earth"}))
            { // Check if primary body is Earth
              radial_gravity(z, tol, deg, &grav);
            }
            else
            {
              grav = deg;
            }
            double dstate[6] = {0.0};
            Gravity(z, &dstate[3], grav, orbit);
            Feval[0] = Feval[0] + pow(grav, 2) / pow(deg, 2);
            G[ID2(cnt + 1, 1, N + 1)] = dstate[3];
            G[ID2(cnt + 1, 2, N + 1)] = dstate[4];
            G[ID2(cnt + 1, 3, N + 1)] = dstate[5];
          }
        }
        std::fill(Gprev.begin(), Gprev.end(), 0.0);
        for (int cnt = 0; cnt <= N; cnt++)
        {
          Gprev[ID2(cnt + 1, 1, N + 1)] = G[ID2(cnt + 1, 1, N + 1)];
          Gprev[ID2(cnt + 1, 2, N + 1)] = G[ID2(cnt + 1, 2, N + 1)];
          Gprev[ID2(cnt + 1, 3, N + 1)] = G[ID2(cnt + 1, 3, N + 1)];
        }
      }

      // Least Squares Acceleration Approximation
      int M;
      M = N;
      std::vector<double> T((M + 1) * N, 0.0);
      std::vector<double> A(N * (M + 1), 0.0);
      std::vector<double> gamma;
      std::vector<double> Gapprox;
      lsq_chebyshev_fit(1.0, N - 1, M, T, A);
      gamma = matmul(A, Gprev, N, M + 1, 3, N, M + 1);
      Gapprox = matmul(T, gamma, M + 1, N, 3, M + 1, N);

      // Check magnitude of last few coefficients
      double max_gamma;
      for (int i = 1; i <= N; i++)
      {
        degree = i;
        max_gamma = 0.0;
        for (int k = 1; k <= 3; k++)
        {
          if (fabs(gamma[ID2(i, k, N)]) > max_gamma)
          {
            max_gamma = fabs(gamma[ID2(i, k, N)]);
          }
        }
        if (max_gamma < coeff_tol)
        {
          fit_check = fit_check + 1;
          if (fit_check == coeff)
          {
            break; // Break if last 3 coeffs are less than the tolerance
          }
        }
      }

      // Reinitialize
      Nprev = N;
      memset(tau, 0.0, ((300) * sizeof(double)));
      memset(times, 0.0, ((300) * sizeof(double)));
      memset(X, 0.0, ((300 * 6) * sizeof(double)));
    }
    seg = seg + 2;
    if (fit_check == coeff)
    {
      break;
    }
  }
  seg = seg - 2;
  return;
}
