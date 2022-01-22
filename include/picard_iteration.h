/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Feb 2017
*  LAST MODIFIED:    Feb 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Header file
*/

#ifndef __PI__
#define __PI__

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include "const.h"
#include <vector>

void picard_iteration(double* Xint, double* Vint, std::vector<double> X, std::vector<double> V, std::vector<double>times, int N, int M, double deg, int hot, double tol,
  std::vector<double> P1, std::vector<double> P2, std::vector<double> T1, std::vector<double> T2, std::vector<double> A, double* Feval, std::vector<double> Alpha, std::vector<double> Beta);

#endif
