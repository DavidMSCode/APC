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
/**
 * @brief Interpolates solution using a default time step
 * 
 * @param ALPHA Position Coefficients
 * @param BETA Velocity Coefficients
 * @param soln_size Length of solution
 * @param coeff_size Length of Coefficients
 * @param N Degree of polynomial
 * @param seg_times t0 and tf for each segment
 * @param W1 Time scale factor 1
 * @param W2 Time scale factor 2
 * @param t0 Initial time (s)
 * @param tf Final time (s)
 * @param dt User desired output time interval
 * @param total_segs Total number of segments
 * @return std::vector<double> The array of interpolated solutions 
 */
std::vector<double> interpolate(std::vector<double>  ALPHA, std::vector<double> BETA, int soln_size, int coeff_size, int N, std::vector<double> seg_times,
  std::vector<double> W1, std::vector<double> W2, double t0, double tf, double dt, int total_segs);

/**
 * @brief Overload function that Interpolates solution using a user provided timestep
 * 
 * @param ALPHA Position Coefficients
 * @param BETA Velocity Coefficients
 * @param coeff_size Length of Coefficients
 * @param N Degree of polynomial
 * @param seg_times t0 and tf for each segment
 * @param W1 Time scale factor 1
 * @param W2 Time scale factor 2
 * @param total_segs Total number of segments
 * @param time_out User defined time output vector
 * @return std::vector<double> The array of interpolated solutions
 */
std::vector<double> interpolate(std::vector<double>  ALPHA, std::vector<double> BETA, int coeff_size, int N, std::vector<double> seg_times,
  std::vector<double> W1, std::vector<double> W2, int total_segs, std::vector<double > time_out);

#endif
