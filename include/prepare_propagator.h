/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Feb 2017
*  LAST MODIFIED:    Feb 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Header file
*/

#ifndef __PP__
#define __PP__

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include "const.h"
#include <vector>

void prepare_propagator(double* r0, double* v0, double t0, double t_final, double dt, double tp, double tol,
  int N, int M, int seg, int* prep_HS, std::vector<double> &t_orig, std::vector<double> &tvec,
  std::vector<double> &P1, std::vector<double> &P2, std::vector<double> &T1, std::vector<double> &T2, std::vector<double> &A, std::vector<double> &Ta);

#endif
