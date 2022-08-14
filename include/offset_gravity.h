/*
*  AUTHORS:          David Stanley (DavidMS4@illinois.edu) based on code by Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Aug 2022
*  LAST MODIFIED:    Aug 2022
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Header file
*/

#ifndef __OFFSET__
#define __OFFSET__

#include "perturbed_gravity.h"

void offset_gravity(double t, double* Xo, double err, int ii, int N, double deg, int hot, double* G, double tol, int* itr, double* Feval, IterCounters& ITRs, double* del_G);

#endif
