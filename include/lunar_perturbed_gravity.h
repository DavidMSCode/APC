/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Feb 2017
*  LAST MODIFIED:    Feb 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Header file
*/

#ifndef __PERT__
#define __PERT__
 
 /**
  * @brief Iteration struct for controlling full or approximate calculations in the pertubed_gravity() function
  * 
  */
struct IterCounters{
    int ITR1 = 0;
    int ITR2 = 0;
    int ITR3 = 0;
    int ITR4 = 0;
    int MODEL = 0;
};

/**
 * @brief  Returns the gravitational acceleration at a position to a desired degree of gravity. Uses the full fidelity model or an approximation depending on the current tolerance and iteration count.
 * 
 * @param t time from epoch (s)
 * @param Xo position in ECEF coordinates (km)
 * @param err error
 * @param i current index along array 
 * @param N number of sampled points
 * @param deg Gravity model degree
 * @param hot Hot start switch condition
 * @param G Gravitational acceleration for output (km/s^2)
 * @param tol Solution tolerance
 * @param itr  Picard iteration counter
 * @param Feval Function evaluations counter
 * @param ITRs Iteration struct for controlling full or approx calculations
 * @param del_G Difference between full and approximate gravity evaluation (km/s^2)
 */
void lunar_perturbed_gravity(double t, double* Xo, double err, int i, int N, double deg, int hot, double* G, double tol, int* itr, double* Feval, IterCounters& ITRs, double* del_G);

/**
 * @brief Returns the approximate gravity acceleration through the J6 zonal harmonic.
 * 
 * @param t time from epoch (s)
 * @param X position in ECEF coordinates (km)
 * @param dX Approximate gravitational acceleration (km/s^2)
 * @param Feval Function evaluations counter
 */
void lunar_Grav_Approx(double t, double* X, double* dX, double* Feval);

/**
 * @brief Returns the full gravity acceleration through the input degree
 * 
 * @param t times from epoch (s)
 * @param Xo position in ECEF coordinates (km)
 * @param acc full gravitational acceleration (km/s^2)
 * @param tol output tolerance
 * @param deg degree of gravity model to use
 * @param Feval Function evaluations counter
 */
void lunar_Grav_Full(double t, double* Xo, double* acc, double tol, double deg, double* Feval);

#endif
