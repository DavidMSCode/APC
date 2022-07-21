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

void picard_iteration(std::vector<double> &X, std::vector<double> &V, std::vector<double> &times, int N, int M, double deg, int hot, double tol,
  std::vector<double> &P1, std::vector<double> &T1, std::vector<double> &A, double* Feval, std::vector<double> &Alpha, Orbit &orb, EphemerisManager ephem);

#endif
