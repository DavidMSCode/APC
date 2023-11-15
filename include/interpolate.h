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
#include "Orbit.h"
/**
 * @brief Interpolates solution using a default time step
 * 
 * @param orbit the orbit object that stores the chebyshev coefficients
 */
std::vector<double> interpolateDefault(Orbit &orbit);

/**
 * @brief Overload function that Interpolates solution using a user provided timestep
 * 
 * @param orbit the orbit objects that stores the chebyshev coefficients
 */
std::vector<double> interpolate(Orbit &orbit);

#endif
