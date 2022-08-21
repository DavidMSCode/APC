/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Jun 2022
*  LAST MODIFIED:    Jun 2022
*  AFFILIATION:      Department of Aerospace Engineering, University of Illinois, Urbana-Champaign
*  DESCRIPTION:      Header file
*/

#ifndef __INT2RAD__
#define __INT2RAD__

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include "const.h"

void inertial2radial(double* r_eci,double* v_eci,double* xrdl,double* yrdl,double* zrdl);
    
#endif