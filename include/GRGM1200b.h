#ifndef __GRGM1200B_H__
#define __GRGM1200B_H__

#include "Orbit.h"

const int Max_Degree_GRGM1200b = 251;
const double Re_GRGM1200b = 1.7380000000000000e+03;                     //Reference radius in km
const double GM_GRGM1200b = 4.9028001224453001e+03;                     //Grvitational Parameter km^3/s^2

/* Macro for looking up indices using "1" based addressing ( a la fortran and matlab )
 and column-major array format [Note ld = Num of rows]
 */
#define IDX2F(i,j,ld) ((((j)-1)*(ld))+((i)-1))
#define IDX3F(i,j,k,ld1,ld2) ((((k)-1)*(ld1*ld2))+(((j)-1)*(ld1))+((i)-1))

/*!
* \brief Gravity Evaluation
* This is the function that evaluates the spherical harmonic
* series to provide acceleration
*
* \param[in] p 3 element position in ECEF
* \param[in] Deg Degree and order of the serries to be used
* \param[out] Gxyz Gravitational Acceloration Output
*/
void GRGM1200b( double* p, double* Gxyz, int DEG);

/*!
* \brief Gravity Potential Evaluationc++ co
* This is the function that evaluates the spherical harmonic
* serries to provide gravitational potential
*
* \param[in] p 3 element position in ECEF
* \param[in] Deg Degree and order of the serries to be used
* \param[out] Pot Gravitational Potential Output
*/
void GRGM1200bPot(double *p, double *Pot, int DEG);
/*!
* \brief Legendre Polyniomial Evaluation
* This is the function that computes the normalized associated
* legendre polynomials based on geocentric latitude
*
* \param[in] phi Geocentric latitude
* \param[in] Deg Degree and order of the serries to be used
* \param[out] P associated Legendre polynomial matrix
* \param[out] scaleFactor_GRGM1200b Legendre scale factor
*/
void loc_gravLegendre_GRGM1200b(double phi, double *scaleFactor_GRGM1200b, double *P, int DEG);
/*!
* \brief Internal Gravitational Acceloration Evaluation
* This is the function that computes the gravitational acceloration based on
* the associated Legendre polynomials and the state
*
* \param[in] p Position vector in ECEF
* \param[in] P associated Legendre polynomial matrix
* \param[in] scaleFactor_GRGM1200b Legendre scale factor
* \param[in] Deg Degree and order of the serries to be used
* \param[in] r Position vector in ECEF
* \param[in] smlambda Trigonometric function of longitude
* \param[in] smlambda Trigonometric function of longitude
* \param[out] Gxyz Gravitational Acceloration Output
*/
void loc_gravityPCPF_GRGM1200b(double *p, double *P, int DEG, double *smlambda, double *cmlambda, double r, double *scaleFactor_GRGM1200b, double *Gxyz);
/*!
* \brief Internal Gravity Potential Evaluation
* This is the function that computes the gravitational acceloration based on
* the associated Legendre polynomials and the state
*
* \param[in] p Position vector in ECEF
* \param[in] P associated Legendre polynomial matrix
* \param[in] scaleFactor_GRGM1200b Legendre scale factor
* \param[in] Deg Degree and order of the serries to be used
* \param[in] r Position vector in ECEF
* \param[in] smlambda Trigonometric function of longitude
* \param[in] smlambda Trigonometric function of longitude
* \param[out] Pot Gravitational Potential Output
*/
void loc_gravityPot_GRGM1200b(double *p, double *P, int DEG, double *smlambda, double *cmlambda, double r, double *scaleFactor_GRGM1200b, double *Pot);

/*!
* \brief Jacobi Integral
* This is the function computes the Jacobi Integral based on position and
* state vector. Needs to be modififed for lunar fixed frame
*
* This is a usefule tool for checking the accuracy of conservitive
* orbit propigation
*
* \param[in] solN State (position and velocity) vector in ECEF
* \param[in] Deg Degree and order of the serries to be used
* \param[out] H Jacobi Integral Output
*/

void jacobiIntegral_GRGM1200b(double t, double *solN, double *H, int Deg,Orbit &orbit);

#endif