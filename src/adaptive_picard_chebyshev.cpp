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
#include <iostream>

#include "adaptive_picard_chebyshev.h"
#include "const.h"
#include "polydegree_segments.h"
#include "prepare_propagator.h"
#include "picard_chebyshev_propagator.h"
#include "interpolate.h"
#include "c_functions.h"
#include "Orbit.h"
#include "Ephemeris.hpp"

void adaptive_picard_chebyshev( double* Feval, std::vector<double> &Soln, Orbit &orbit, EphemerisManager ephem){

  /* 1. DETERMINE DEGREE/SEGMENTATION SCHEME
  Compute the polynomial degree and number of segments per orbit that will
  result in a solution that satisfies the user specified tolerance. */
  polydegree_segments(orbit,Feval);
  std::cout << "Nodes per segment: " << orbit.N << std::endl;
  std::cout << "Segments per orbit: " << orbit.seg << std::endl; 
  // Array size for coefficients and solution
  double &tf = orbit._Integrator_tf;
  orbit.coeff_size = int((tf/orbit.Period + 1.0)*(orbit.seg+2.0)*(orbit.N+1));
  /* 2. PREPARE PROPAGATOR
  Compute and store the begin and end times for each segment (based on true
  anomaly segmentation) and load the constant matrices corresponding to N. */
  int prep_HS = -1;         // Hot start switch condition
  orbit.prep_HS = prep_HS;
  prepare_propagator(orbit);

  /* 3. PICARD-CHEBYSHEV PROPAGATOR
  Propagate from t0 to tf, iterating on each segment (Picard Iteration), until
  completion. Stores solution as chebyshev nodes for each segment*/
  orbit.M = orbit.N;                                          // Number of Chebyshev nodes
  orbit.CC.A.resize((orbit.coeff_size*3),0.0);
  orbit.CC.B.resize((orbit.coeff_size*3),0.0);
  orbit.CC.total_segs = 0;                                    //running count of orbital path segments
  //reserve space in vectors for interpolation and solution
  int sz = int(ceil(tf/orbit.Period)*(orbit.seg+1));              //ensure sufficient space by overestimating the number of segments
  orbit.segment_end_times.resize(sz,0.0);
  orbit.W1.resize(sz,0.0);
  orbit.W2.resize(sz,0.0);

 picard_chebyshev_propagator(orbit.coeff_size,Feval, orbit, ephem);

  // /* 4. INTERPOLATE SOLUTION
  // The Chebyshev coefficients from each of the orbit segments are used to compute
  // the solution (position & velocity) at the user specified times. */
  if(orbit.USER_TIME)
  {
    //Interpoalte with user defined time vec
    Soln = interpolate(orbit);
  }
  else
  {
    //Interpolate with default dt spacing
    Soln = interpolateDefault(orbit);
    }
  return;
}
