/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Feb 2017
*  LAST MODIFIED:    Feb 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Header file
*/

#ifndef __LSQ__
#define __LSQ__


void lsq_chebyshev_fit(double s, int N, int M, std::vector<double> &T, std::vector<double> &A);

#endif
