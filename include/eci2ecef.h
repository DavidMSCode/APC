/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Jun 2016
*  LAST MODIFIED:    Jun 2016
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Header file
*/

#ifndef __ECI__
#define __ECI__
#include <vector>
void eci2ecef(double t, double* X, double* V, double* xB, double* vB);


void eci2ecef(double t, std::vector<double> X, std::vector<double> V, std::vector<double> A, std::vector<double> &xB, std::vector<double> &vB, std::vector<double> &aB);

/**
 * @brief Converts ECI (Earth-Centered Inertial) coordinates to ECEF (Earth-Centered Earth-Fixed) coordinates, but do not include transport terms for velocity or acceleration.
 */
void eci2ecef_notransport(double t, std::vector<double> X, std::vector<double> V, std::vector<double> A, std::vector<double> &xB, std::vector<double> &vB, std::vector<double> &aB);

#endif
