/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Feb 2017
*  LAST MODIFIED:    Feb 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Header file
*/

#ifndef __APC__
#define __APC__
#include "Orbit.h"
#include "Ephemeris.hpp"

/**
 * @brief This function controls the flow of the APC propagator from orbit segmentation to interpolating the final output vectors.
 * 
 * @param r0 Initial position in ECEF coordinates (km)
 * @param v0 Initial velocity in ECEF coordinates (km/s)
 * @param t0 Initial time from epoch
 * @param tf Final time from epoch
 * @param dt Default timestep
 * @param deg degree of the gravity model
 * @param tol tolerance for the solution
 * @param soln_size size of interpolated output
 * @param Feval Gravity function evaluations counter
 * @param Soln The solution array
 * @param orb The orbit object that stores various properties of the orbit
 * @param ephem the cached ephemeris object for interpolatinf moon and sun positions.
 * @return std::vector<std::vector<double> > Returns a vector of interpolated solution values
 */
std::vector<std::vector<double> > adaptive_picard_chebyshev(double* r0,double* v0, double t0, double tf, double dt, double deg, double tol, int soln_size, double* Feval, std::vector<double> &Soln, Orbit &orb, EphemerisManager ephem);

#endif
