/*
 *  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
 *  DATE WRITTEN:     Feb 2017
 *  LAST MODIFIED:    Feb 2017
 *  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
 *  DESCRIPTION:      Header file
 */

#ifndef __PERT__
#define __PERT__
#include <string>
/**
 * @brief Iteration struct for controlling full or approximate calculations in the pertubed_gravity() function
 *
 */
struct IterCounters
{
  int ITR1;
  int ITR2;
  int ITR3;
  int ITR4;
  int MODEL;
  bool NEXT_FULL;
  bool DID_FULL;
  double TOL_SCALE;
  double SF;
  IterCounters(): ITR1(0), ITR2(0), ITR3(0), ITR4(0), MODEL(0), NEXT_FULL(false), DID_FULL(false), TOL_SCALE(1e-2), SF(1){};
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
void perturbed_gravity(double t, double *Xo, double err, int i, int N, double deg, int hot, double *G, double tol, int *itr, double *Feval, IterCounters &ITRs, double *del_G);

/**
 * @brief Returns the gravitational acceleration using the J2-J6 zonal harmonics of the Earth given a position in ECEF coordinates.
 *
 * @param X ECEF coordinate position
 * @param dX Gravitational acceleration
 */
void Grav_EarthJ2J6(double *X, double *dX);

/**
 * @brief Returns the approximate gravity acceleration through the J6 zonal harmonic.
 *
 * @param t time from epoch (s)
 * @param X position in ECEF coordinates (km)
 * @param dX Approximate gravitational acceleration (km/s^2)
 * @param Feval Function evaluations counter
 */
void Grav_Approx(double t, double *X, double *dX, double *Feval);

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
void Grav_Full(double t, double *Xo, double *acc, double tol, double deg, double *Feval);

bool fullgravswitch(double err, double tol, int hot, int itr, IterCounters &ITRs);

void variableGrav(double t, double *Xo, double *G, double *del_G, double tol, int i, double deg, double *Feval, bool fullgravswitch, std::string primary);


#endif
