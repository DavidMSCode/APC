/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Feb 2017
*  LAST MODIFIED:    Feb 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Header file
*/

#ifndef __INT__
#define __INT__

#include <vector>

std::vector<double> interpolate(std::vector<double>  ALPHA, std::vector<double> BETA, int soln_size, int coeff_size, int N, std::vector<double> seg_times,
  std::vector<double> W1, std::vector<double> W2, double t0, double tf, double dt, int total_segs, int back_prop);

#endif
