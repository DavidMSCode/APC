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
 * @param ALPHA 
 * @param BETA 
 * @param soln_size 
 * @param coeff_size 
 * @param N 
 * @param seg_times 
 * @param W1 
 * @param W2 
 * @param t0 
 * @param tf 
 * @param dt 
 * @param total_segs 
 * @return std::vector<double> 
 */
std::vector<double> interpolate(std::vector<double>  ALPHA, std::vector<double> BETA, int soln_size, int coeff_size, int N, std::vector<double> seg_times,
  std::vector<double> W1, std::vector<double> W2, double t0, double tf, double dt, int total_segs);

/**
 * @brief Interpolates solution using a user provided timestep
 * 
 * @param ALPHA 
 * @param BETA 
 * @param coeff_size 
 * @param N 
 * @param seg_times 
 * @param W1 
 * @param W2 
 * @param total_segs 
 * @param time_out 
 * @return std::vector<double> 
 */
std::vector<double> interpolate(std::vector<double>  ALPHA, std::vector<double> BETA, int coeff_size, int N, std::vector<double> seg_times,
  std::vector<double> W1, std::vector<double> W2, int total_segs, std::vector<double > time_out);

#endif
