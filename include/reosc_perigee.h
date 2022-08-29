/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     May 2017
*  LAST MODIFIED:    May 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Header file
*/

#ifndef __PER__
#define __PER__

#include <vector>

void reosc_perigee(std::vector<double> &X, std::vector<double> &V, std::vector<double> &times, std::vector<double> &Alpha, std::vector<double> &Beta,
  double tf, double t_final, std::vector<double> &t_orig, int N, int M, int* k, int seg, int* prep_HS,
  double tol, double* orb_end, std::vector<double> &tvec, double* r0, double* v0, double t_end, int back_prop);

#endif
