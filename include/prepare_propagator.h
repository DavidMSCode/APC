/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Feb 2017
*  LAST MODIFIED:    Feb 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Header file
*/

#ifndef __PP__
#define __PP__

#include <vector>
#include "Orbit.h"
/**
 * @brief Prepares arrays for storing iteration computations
 * 
 * @param r0 Initial position (km)
 * @param v0 Initial velocity (km/s)
 * @param t0 Initial time (s)
 * @param t_final Final time (s)
 * @param dt Time interval (s)
 * @param tp Approximate time of perigee passage (s)
 * @param tol Tolerance
 * @param N Polynomial degree
 * @param M Sample points
 * @param seg Segments per orbit
 * @param prep_HS Hot start switch condition
 * @param t_orig Segment start and end times for first segment (s)
 * @param tvec Segment start and end times (s)
 * @param P1 First integration operator
 * @param P2 Second integrationi operator 
 * @param T1 Chebyshev velocity matrix
 * @param T2 Chebyshev position matrix
 * @param A Least squares operator
 * @param Ta Chebyshev acceleration matrix
 */
void prepare_propagator(double* r0, double* v0, double t0, double t_final, double dt, double tp, double tol,
  int N, int M, int seg, int* prep_HS, std::vector<double> &t_orig, std::vector<double> &tvec,
  std::vector<double> &P1, std::vector<double> &P2, std::vector<double> &T1, std::vector<double> &T2, std::vector<double> &A, std::vector<double> &Ta, Orbit &orbit);

#endif
