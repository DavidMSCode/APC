/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Jun 2022
*  LAST MODIFIED:    Jun 2022
*  AFFILIATION:      Department of Aerospace Engineering, University of Illinois, Urbana-Champaign
*  DESCRIPTION:      Header file
*/

#ifndef __MEEDOT__
#define __MEEDOT__

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include "const.h"
#include <vector>

void mee_dot(std::vector<double> &MEE, std::vector<double> &a_lvlh, std::vector<double> &MEE_dot, int M);

#endif