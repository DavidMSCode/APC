/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Feb 2017
*  LAST MODIFIED:    Feb 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Header file
*/

#ifndef __PP__
#define __PP__

#include <vector>
#include "Orbit.h"

/**
 * @brief 
 * 
 * @param tol 
 * @param prep_HS 
 * @param orbit 
 */
void prepare_propagator(double tol, int* prep_HS, Orbit &orbit);

#endif
