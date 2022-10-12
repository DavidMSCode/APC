/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Feb 2017
*  LAST MODIFIED:    Feb 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Header file
*/

#ifndef __PDS__
#define __PDS__

/**
 * @brief Computes the number of segments per orbit and the degree of the polynomial required to fit acceleration to the user specified accuracy
 * 
 * @param r0 Initial position vector (km)
 * @param v0 Initial velocity vector (km/s)
 * @param deg Gravity Degree (max 100)
 * @param tol Tolerance
 * @param Feval Function evaluation counter
 * @param seg Segments per orbit
 * @param degree Degree of Chebyshev Polynomial
 * @param tp Approximate time of Keplerian perigee passage (s)
 * @param Period Keplerian orbit period (s)
 */
void polydegree_segments( double* r0, double* v0, double deg, double tol, double* Feval, int* seg, int* degree, double* tp, double* Period );

#endif
