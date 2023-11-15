/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Feb 2017
*  LAST MODIFIED:    Feb 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Header file
*/

#ifndef __PROP__
#define __PROP__

#include <vector>

#include "Orbit.h"
#include "Ephemeris.hpp"

/**
 * @brief This function iterates through each segment of the predicted orbit and solves the picard integral to the desired degree of precision for each chebyshev node
 * 
 * @param r0 Initial position vector (km)
 * @param v0 Initial velocity vector (km/s)
 * @param t0 Iniitla time (s)
 * @param t_final Final time (s)
 * @param deg Gravity degree
 * @param tol Tolerance
 * @param Period Orbital period (s)
 * @param tvec Start and end times for all segments (s)
 * @param t_orig First segment start and end times
 * @param seg Number of segments per orbit
 * @param N Degree of Chebyshev polynomials
 * @param M Number of Chebyshev nodes
 * @param prep_HS Hot start switch flag
 * @param coeff_size Length of Chebyshev coefficient array
 * @param soln_size Size of solution vectors
 * @param total_seg Total number of segments
 * @param P1 First integrator operator
 * @param P2 Second integrator operator
 * @param T1 Chebyshev position matrix
 * @param T2 Chebyshev velocity matrix
 * @param A Least squares operator
 * @param Ta Chebyshev acceleration matrix
 * @param W1 Segment time scale factor 1
 * @param W2 Segment time scale factor 2
 * @param Feval Function evaluation counter
 * @param ALPHA Output: Position coefficients
 * @param BETA Output: Velocity coefficicents
 * @param segment_times Array of segment start and end times
 * @param orb Orbit object that stores properties of orbit
 * @param ephem Ephemeris manager that provides position of Sun and Moon
 * @return std::vector<std::vector<double> > 
 */
// std::vector<std::vector<double> >  picard_chebyshev_propagator(double* r0, double* v0, double t0, double t_final,double deg, double tol, double Period,
//    std::vector<double> &tvec, std::vector<double> &t_orig, int seg, int N, int M, int* prep_HS, int coeff_size, int soln_size, int* total_seg,
//    std::vector<double> &P1, std::vector<double> &P2, std::vector<double> &T1, std::vector<double> &T2, std::vector<double> &A, std::vector<double> &Ta, std::vector<double> &W1, std::vector<double> &W2, double* Feval,
//    std::vector<double> &ALPHA, std::vector<double> &BETA, std::vector<double> &segment_times, Orbit &orbit, EphemerisManager ephem);
void picard_chebyshev_propagator(int coeff_size, double* Feval, Orbit &orbit, EphemerisManager ephem);

#endif
