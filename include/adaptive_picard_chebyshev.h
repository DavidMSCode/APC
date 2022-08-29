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

std::vector<std::vector<double> > adaptive_picard_chebyshev(double* r0,double* v0, double t0, double tf, double dt, double deg,
double tol, int soln_size, double* Feval, std::vector<double> &Soln, Orbit &orb, EphemerisManager ephem, double t_start, double t_end, int back_prop);

#endif
