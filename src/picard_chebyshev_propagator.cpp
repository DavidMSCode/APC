/*
 *  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
 *  DATE WRITTEN:     May 2017
 *  LAST MODIFIED:    May 2017
 *  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
 *  DESCRIPTION:      Propagates the initial conditions
 *
 * INPUT:
 *    r0            -- Initial position vector (km)
 *    v0            -- Initial velocity vector (km/s)
 *    t0            -- Initial time (s)
 *    t_final       -- Final time (s)
 *    deg           -- Gravity Degree (max 100)
 *    tol           -- Tolerance
 *    Period        -- Period (s)
 *    tvec          -- Segment start and end times (s)
 *    t_orig        -- Segment start and end times for first segment (s)
 *    seg           -- Segments per orbit
 *    N             -- Polynomial degree
 *    M             -- Sample points
 *    prep_HS       -- Hot start switch function
 *    coeff_size    -- Length of coefficient array
 *    soln_size     -- Size of solution array
 *    P1            -- First integration operator
 *    P2            -- Second integration operator
 *    T1            -- Chebyshev position matrix
 *    T2            -- Chebyshev velocity matrix
 *    A             -- Least squares operator
 *    Ta            -- Chebyshev acceleration matrix
 *    W1            -- Time scale factor 1
 *    W2            -- Time scale factor 2
 *    Feval         -- Function evaluation counter
 *
 * OUTPUTS:
 *    total_seg     -- Total segments
 *    ALPHA         -- Position coefficients
 *    BETA          -- Velocity coefficients
 *    segment_times -- Array of segment start and end times
 */

#include <math.h>

#include <string>
#include <vector>
#include <iostream>

#include "picard_chebyshev_propagator.h"
#include "const.h"
#include "prepare_propagator.h"
#include "picard_iteration.h"
#include "FandG.h"
#include "reosc_perigee.h"
#include "c_functions.h"
#include "Orbit.h"
#include "Ephemeris.hpp"
#include "flags.h"
#include "SpiceUsr.h"
using namespace std;

void picard_chebyshev_propagator(int coeff_size, double *Feval, Orbit &orbit, EphemerisManager &ephem)
{
  if (g_DEBUG_PICARD)
  {
    cout << "Entering picard_chebyshev_propagator" << std::endl;
  }

  // Initialize variables
  int loop = 0;                    // Break loop condition
  int M = orbit.M;                 // Number of Chebyshev nodes
  int N = orbit.N;                 // Chebyshev polynomial degree
  int &k = orbit.k;                // Counter: segments per orbit
  int &hot = orbit.prep_HS;        // Hot start switch
  int &seg_cnt = orbit.total_segs; // Counter: total segments
  int seg = orbit.seg;             // Segments per orbit
  double mu = orbit.GetPrimaryGravitationalParameter();
  double w1, w2, tf;
  vector<double> &ALPHA = orbit.CC.A;
  vector<double> &BETA = orbit.CC.B;
  vector<double> &segment_times = orbit.segment_end_times;
  // Temporary segment variables
  vector<double> &Alpha = orbit.Alpha_seg;
  vector<double> &Beta = orbit.Beta_seg;
  vector<double> &X = orbit.X_seg;
  vector<double> &V = orbit.V_seg;

  // Store initial position and velocity in temporary arrays for the integrator
  for (int i = 0; i < 3; i++)
  {
    orbit.r0_seg[i] = orbit._In_r0[i];
    orbit.v0_seg[i] = orbit._In_v0[i];
  }

  // PROPAGATION
  // #pragma omp critical(PI)
  // {
  while (loop == 0)
  {
    if (g_DEBUG_PICARD)
    {
      cout << "===" << std::endl;
      cout << "Starting Segment " << seg_cnt + 1 << "." << std::endl;
      double debug_altitude = 0.0;
      double debug_velocity = 0.0;
      for (int j = 0; j < 3; j++)
      {
        debug_altitude += orbit.r0_seg[j] * orbit.r0_seg[j];
        debug_velocity += orbit.v0_seg[j] * orbit.v0_seg[j];
      }
      debug_altitude = sqrt(debug_altitude);
      debug_velocity = sqrt(debug_velocity);
      cout << "debug_altitude: " << debug_altitude << std::endl;
      cout << "debug_velocity: " << debug_velocity << std::endl;
    }
    // A. CALCULATE CHEBYSHEV NODES
    //  Compute cosine time vector for the current segment
    double t0 = orbit.tvec[k];
    tf = orbit.tvec[k + 1];
    while (tf == 0.0)
    {
      k = k + 1;
      t0 = orbit.tvec[k];
      tf = orbit.tvec[k + 1];
    }
    if (tf > orbit._Integrator_tf)
    {
      tf = orbit._Integrator_tf;
    }
    w1 = (tf + t0) / 2.0;
    w2 = (tf - t0) / 2.0;
    orbit.W1[seg_cnt] = w1;
    orbit.W2[seg_cnt] = w2;

    double z[6] = {0.0};

    orbit.tau.resize(M + 1, 0.0);
    orbit.times_seg.resize(M + 1, 0.0);
    orbit.X_seg.resize((M + 1) * 3, 0.0);
    orbit.V_seg.resize((M + 1) * 3, 0.0);
    orbit.Beta_seg.resize(N * 3, 0.0);
    orbit.Alpha_seg.resize((N + 1) * 3, 0.0);

    // B. KEPLERIAN WARM START
    double et;                               // Ephermeris time
    string primary = orbit.GetPrimaryBody(); // The body of primary orbit
    string InertFrame = orbit.GetInertFrame();
    string InertCenter = orbit.GetInertCenter();
    double p_state[6] = {0.0};
    double lt;
    // Initial state
    for (int j = 0; j < 3; j++)
    {
      orbit.z0_seg[j] = orbit.r0_seg[j];
      orbit.z0_seg[j + 3] = orbit.v0_seg[j];
    }
    for (int cnt = 0; cnt <= M; cnt++)
    {
      // Warm start
      // Use F and G equations to get two body guess for inertial position and velocity  relative to the primary body
      orbit.tau[cnt] = -cos(cnt * C_PI / M);
      orbit.times_seg[cnt] = orbit.tau[cnt] * w2 + w1;
      et = orbit.et(orbit.times_seg[cnt]);
      FandG(orbit.z0_seg, z, orbit.times_seg[cnt] - t0, mu);
      // GET STATE OF PRIMARY BODY RELATIVE TO INERTIAL FRAME CENTER
      // spkezr_c(primary.c_str(),et,InertFrame.c_str(),"LT+S",InertCenter.c_str(),p_state,&lt);
      // Add two body to primary body state to get the inertial state for the satellite relative to the inertial frame center
      X[ID2(cnt + 1, 1, M + 1)] = z[0] + p_state[0];
      X[ID2(cnt + 1, 2, M + 1)] = z[1] + p_state[1];
      X[ID2(cnt + 1, 3, M + 1)] = z[2] + p_state[2];
      V[ID2(cnt + 1, 1, M + 1)] = z[3] + p_state[3];
      V[ID2(cnt + 1, 2, M + 1)] = z[4] + p_state[4];
      V[ID2(cnt + 1, 3, M + 1)] = z[5] + p_state[5];
    }

    // PICARD ITERATION
    picard_iteration(Feval, orbit, ephem);
    for (int cnt = 0; cnt <= M; cnt++)
    {
      // GET STATE OF PRIMARY BODY RELATIVE TO INERTIAL FRAME CENTER
      orbit.times_seg[cnt] = orbit.tau[cnt] * w2 + w1;
      et = orbit.et(orbit.times_seg[cnt]);
      // spkezr_c(primary.c_str(),et,InertFrame.c_str(),"LT+S",InertCenter.c_str(),p_state,&lt);
      // Add two body to primary body state to get the inertial state for the satellite relative to the primary body
      X[ID2(cnt + 1, 1, M + 1)] -= p_state[0];
      X[ID2(cnt + 1, 2, M + 1)] -= p_state[1];
      X[ID2(cnt + 1, 3, M + 1)] -= p_state[2];
      V[ID2(cnt + 1, 1, M + 1)] -= p_state[3];
      V[ID2(cnt + 1, 2, M + 1)] -= p_state[4];
      V[ID2(cnt + 1, 3, M + 1)] -= p_state[5];
    }
    // Loop exit condition (if tf is within 1e-12 of the final time)
    if (fabs(tf - orbit._Integrator_tf) / tf < 1e-12)
    {
      loop = 1;
    }
    if (orbit.suborbital)
    {
      loop = 1;
    }
    // ASSIGN NEXT SEGMENT INITIAL CONDITIONS
    // Assign new initial conditions for next segment. (The final conditions of the current segment)
    for (int j = 1; j <= 3; j++)
    {
      orbit.r0_seg[j - 1] = X[ID2(M + 1, j, M + 1)];
      orbit.v0_seg[j - 1] = V[ID2(M + 1, j, M + 1)];
    }

    /* REOSCULATE PERIGEE
    Reosculate Keplerian perigee after each orbit propagation. If this is not
    done then the precomputed segment times no longer align with true segment
    breaks (due to perturbations) and the precomputed segment scheme fails to
    produce a solution that satisfies the required tolerance. This effect
    increases with increasing eccentricity. */

    double &orb_end = orbit.orb_end;
    orb_end = 0.0;
    reosc_perigee(tf, orbit);
    // Segments per orbit counter
    k = k + 1;

    // STORE TRAJECTORY COEFFICIENTS
    for (int i = 1; i <= N; i++)
    {
      BETA[ID2(i + (seg_cnt * N), 1, coeff_size)] = Beta[ID2(i, 1, N)];
      BETA[ID2(i + (seg_cnt * N), 2, coeff_size)] = Beta[ID2(i, 2, N)];
      BETA[ID2(i + (seg_cnt * N), 3, coeff_size)] = Beta[ID2(i, 3, N)];
    }
    for (int i = 1; i <= N + 1; i++)
    {
      // Store X and V points
      ALPHA[ID2(i + seg_cnt * (N + 1), 1, coeff_size)] = Alpha[ID2(i, 1, N + 1)];
      ALPHA[ID2(i + seg_cnt * (N + 1), 2, coeff_size)] = Alpha[ID2(i, 2, N + 1)];
      ALPHA[ID2(i + seg_cnt * (N + 1), 3, coeff_size)] = Alpha[ID2(i, 3, N + 1)];
    }

    segment_times[seg_cnt + 1] = tf;
    if (orb_end != 0.0)
    {
      segment_times[seg_cnt + 1] = orb_end;
    }

    // debug print segment start and end time
    if (g_DEBUG_PICARD)
    {
      cout << "Segment end time: " << segment_times[seg_cnt + 1] << std::endl;
    }

    // Total segments counter
    seg_cnt = seg_cnt + 1;
  }
  // debug print end of integration reached
  if (g_DEBUG_PICARD)
  {
    cout << "===" << std::endl;
    cout << "Last segment reached." << std::endl;
  }
  return;
}
