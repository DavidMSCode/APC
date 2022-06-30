/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Feb 2017
*  LAST MODIFIED:    Feb 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Header file
*/

#ifndef __CC__
#define __CC__

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>

#include "const.h"

void clenshaw_curtis_ivpII( int N, int M, std::vector<double> &T2, std::vector<double> &P2, std::vector<double> &T1, std::vector<double> &P1, std::vector<double> &Ta, std::vector<double> &A );

#endif
