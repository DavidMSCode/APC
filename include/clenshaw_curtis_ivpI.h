/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Jul 2022
*  LAST MODIFIED:    Jul 2022
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Header file
*/

#ifndef __CC__
#define __CC__

#include <vector>

void clenshaw_curtis_ivpI( int N, int M, std::vector<double> &T1, std::vector<double> &P1, std::vector<double> &Ta, std::vector<double> &A );

#endif
