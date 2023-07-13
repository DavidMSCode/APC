/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     May 2017
*  LAST MODIFIED:    May 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Adaptive Picard-Chebyshev Numerical Integration
*
* INPUT:
*    r0        -- Initial position vector (km)
*    v0        -- Initial velocity vector (km/s)
*    t0        -- Initial time (s)
*    tf        -- Final time (sec)
*    dt        -- Solution Output Time Interval (s)
*    deg       -- Gravity Degree (max 100)
*    tol       -- Tolerance
*    soln_size -- Length of solution array
*    Feval     -- Function evaluation counter
*
* OUTPUTS:
*    Soln      -- [Position (km), Velocity (km/s)] solution at user specified times
*
* REFERENCES:
* 1. Junkins, J.L., and Woollands, R., "Nonlinear Differential Equations Solvers via Adaptive Picard-Chebyshev Iteration: Applications in Astrodynamics",
*    AAS/AIAA Astrodynamics Specialist Conference, Stevenson, WA, 2017.
* 2. Junkins, J.L., and Woollands, R., "Nonlinear Differential Equation Solvers via Adaptve Picard-Chebyshev Iteration: Applications in Astrodynamics",
*    JGCD, submitted 2017.
*
* COMMENTS:
*
*/


#include <math.h>
#include <vector>

#include "adaptive_picard_chebyshev.h"
#include "const.h"
#include "polydegree_segments.h"
#include "prepare_propagator.h"
#include "picard_chebyshev_propagator.h"
#include "interpolate.h"
#include "c_functions.h"
#include "Orbit.h"
#include "Ephemeris.hpp"

std::vector<std::vector<double> > adaptive_picard_chebyshev(double* r0,double* v0, double t0, double tf, double dt, double deg, double tol, int soln_size, double* Feval, std::vector<double> &Soln, Orbit &orbit, EphemerisManager ephem){

  /* 1. DETERMINE DEGREE/SEGMENTATION SCHEME
  Compute the polynomial degree and number of segments per orbit that will
  result in a solution that satisfies the user specified tolerance. */
  polydegree_segments(orbit,deg,tol,Feval);

  // Array size for coefficients and solution
  orbit.coeff_size = int((tf/orbit.Period + 1.0)*(orbit.seg+2.0)*(orbit.N+1));
  /* 2. PREPARE PROPAGATOR
  Compute and store the begin and end times for each segment (based on true
  anomaly segmentation) and load the constant matrices corresponding to N. */
  int prep_HS = -1;         // Hot start switch condition
  prepare_propagator(tol,&prep_HS,orbit);

  /* 3. PICARD-CHEBYSHEV PROPAGATOR
  Propagate from t0 to tf, iterating on each segment (Picard Iteration), until
  completion. */
  orbit.CC.A.resize((orbit.coeff_size*3),0.0);
  orbit.CC.B.resize((orbit.coeff_size*3),0.0);
  orbit.CC.total_segs = 0;
  int sz = int(ceil(1.0*tf/orbit.Period*orbit.seg))+2;
  std::vector<double> segment_times(sz,0.0);
  std::vector<double> W1(sz,0.0);
  std::vector<double> W2(sz,0.0);
  std::vector<std::vector<double> > states;

  states = picard_chebyshev_propagator(r0,v0,t0,tf,deg,tol,orbit.Period,tvec,t_orig,orbit.seg,orbit.N,M,&prep_HS,orbit.coeff_size,soln_size,&total_seg,
    P1,P2,T1,T2,A,Ta,W1,W2,Feval,ALPHA,BETA,segment_times, orbit, ephem);

  //store chebyshev in orbit object
  orbit.SetCC(ALPHA,BETA,W1,W2,N,orbit.coeff_size,segment_times,tf,t0,total_seg);
  // /* 4. INTERPOLATE SOLUTION
  // The Chebyshev coefficients from each of the orbit segments are used to compute
  // the solution (position & velocity) at the user specified times. */
  if(orbit.USER_TIME)
  {
    //Interpoalte with user defined time vec
    Soln = interpolate(ALPHA,BETA,orbit.coeff_size,orbit.N,segment_times,W1,W2,total_seg,orbit.T);
  }
  else
  {
    //Interpolate with default dt spacing
    Soln = interpolate(ALPHA,BETA,soln_size,orbit.coeff_size,orbit.N,segment_times,W1,W2,t0,tf,dt,total_seg);
    int len;
    len = int(ceil(tf/dt));
    std::vector<double> time_out(len+1,0.0);
    time_out[0] = t0;
    for (int ii=1; ii<=len; ii++){
      double time = time_out[ii-1] + dt;
      if(time>tf)
      {
        time = tf;
      }
      time_out[ii] = time;
    }
    orbit.SetTimeVec(time_out);
    }
  return states;
}
