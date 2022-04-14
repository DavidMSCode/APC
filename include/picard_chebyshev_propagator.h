/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Feb 2017
*  LAST MODIFIED:    Feb 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Header file
*/

#ifndef __PROP__
#define __PROP__

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include "const.h"
#include <Orbit.h>
#include <vector>
#include "Ephemeris.hpp"

std::vector<std::vector<double> >  picard_chebyshev_propagator(double* r0, double* v0, double t0, double t_final,double deg, double tol, double Period,
   std::vector<double> &tvec, std::vector<double> &t_orig, int seg, int N, int M, int* prep_HS, int coeff_size, int soln_size, int* total_seg,
   std::vector<double> &P1, std::vector<double> &P2, std::vector<double> &T1, std::vector<double> &T2, std::vector<double> &A, std::vector<double> &Ta, std::vector<double> &W1, std::vector<double> &W2, double* Feval,
   std::vector<double> &ALPHA, std::vector<double> &BETA, std::vector<double> &segment_times, Orbit &orb, EphemerisManager ephem);

#endif
