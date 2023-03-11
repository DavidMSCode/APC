/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Feb 2017
*  LAST MODIFIED:    Feb 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Header file
*/

#ifndef __PI__
#define __PI__


#include "Orbit.h"
#include "Ephemeris.hpp"

/**
 * @brief Performs a single iteration over an orbital segment
 * 
 * @param Xint Initial position (km)
 * @param Vint Initial velocity (km/s)
 * @param X Position warm start for current segment (km)
 * @param V Velocity warm start for current segment (km/s)
 * @param times Time array for current segment (s)
 * @param N Degree of Chebyshev polynomial
 * @param M Number of sample points
 * @param deg Degree of Gravity model
 * @param hot Hot start on/off switch condition
 * @param tol Tolerance
 * @param P1 First integration operator (Acceleration to Velocity)
 * @param P2 Second integration operator (Velocity to Position)
 * @param T1 First Chebyshev matrix
 * @param T2 Second Chebyshev matrix
 * @param A Least squares operator
 * @param Feval Function evaluation counter
 * @param Alpha Position coefficients for current segment
 * @param Beta Velocity coefficients for current segment
 * @param orb Orbit object that stores properties of orbit
 * @param ephem Ephemeris manager that provides position of Sun and Moon
 */
void picard_iteration(double* Xint, double* Vint, std::vector<double> &X, std::vector<double> &V, std::vector<double> &times, int N, int M, double deg, int hot, double tol,
  std::vector<double> &P1, std::vector<double> &P2, std::vector<double> &T1, std::vector<double> &T2, std::vector<double> &A, double* Feval, std::vector<double> &Alpha, std::vector<double> &Beta, Orbit &orbit, EphemerisManager &ephem);

#endif
