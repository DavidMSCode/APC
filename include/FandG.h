
#ifndef _F_and_G_
#define _F_and_G_

#include "const.h"
/*! \brief Calculates the analytic solution to the unperturbed orbit problem using F&G solution
 *   Used as a "warm-start" for MCPI.
 *  \param[in] z0  Initial state vector [r0,v0] in units of [e.r] and [e.r/ctu]
 *  \param[out] z  Output state vector [rf,vf] in units of [e.r] and [e.r/ctu]
 *  \param[in] dt  Desired time of output in units of [ctu]
*/
void FandG( const double* z0, double* zf, const double dt, const double mu);
/*! \brief Newton solver to solve Kepler's equation in the F&G solution.
 *  \param[in] a     Semimajor axis [e.r]
 *  \param[in] dt    Desired time of solution [ctu]
 *  \param[in] rMag  Magnitude of the radius vector [e.r.]
 *  \param[in] sig0  Measure of orthogonality between instantaneous positiona and velocity vector
 *  \param[in] tol   Convergence threshold
 *  \param[out] Ehat Change in eccentric anomaly
*/
double newtonFandG( const double a, const double dt, const double rMag,
                    const double sig0, const double tol, const double mu);

#endif
