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

 * @param Feval Function evaluation counter
 * @param orbit Object that stores properties of the computed orbit
 * @param ephem Ephemeris manager that provides position of Sun and Moon
 */
void picard_iteration(double* Feval, Orbit &orbit, EphemerisManager &ephem);

void picardSegmentGravity(vector<double> times, vector<double> X, vector<double> &accel, vector<double> &del_G, double err, int M, int deg, int lowdeg, int hot, double tol, int itr, double *Feval, IterCounters &ITRs, Orbit &orbit);
/**
 * @brief Computes the gravity acceleration for the whole segment in parallel
*/

#endif
