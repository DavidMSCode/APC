/**
 * @file EphemerisRotation.cpp
 * @author David Stanley (davidms4@illinois.edu)
 * @brief 
 * @version 0.1
 * @date 2023-03-26
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef __EPHEMROT__
#define __EPHEMROT__

#include "Orbit.h"

using namespace std;

void BodyFixedToInertial(vector<double> xF, vector<double> vF, vector<double> &xI, vector<double> &vI, double t, Orbit &orbit);
void BodyFixedToInertial(double* xF, double* vF, double* xI, double *vI, double t, Orbit &orbit);

void BodyFixedAccelerationToInertial(double* aF, double* aI, double t, Orbit &orbit);

void InertialToBodyFixed(vector<double> xI, vector<double> vI, vector<double> &xF, vector<double> &vF, double t, Orbit &orbit);
void InertialToBodyFixed(double* xI, double* vI, double* xF, double* vF, double t, Orbit &orbit);

#endif