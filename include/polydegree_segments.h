/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Feb 2017
*  LAST MODIFIED:    Feb 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Header file
*/

#ifndef __PDS__
#define __PDS__

#include "Orbit.h"

/**
 * @brief Computes the number of segments per orbit and the degree of the polynomial required to fit acceleration to the user specified accuracy
 * 
 * @param orbit Orbit class object
 * @param deg Gravity Degree (max 100)
 * @param Feval Function evaluation counter
 */
void polydegree_segments(Orbit &orbit, double* Feval);

#endif
