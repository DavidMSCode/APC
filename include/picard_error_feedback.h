/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     Feb 2017
*  LAST MODIFIED:    Feb 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Header file
*/

#ifndef __PEF__
#define __PEF__

/**
 * @brief Generates a linear correction term for acceleration based on the position error relative to the previous iteration
 * 
 * @param X The position (km)
 * @param del_X  The position error (km)
 * @param del_a  The acceleration correction
 */
void picard_error_feedback(double* X, double* del_X, double* del_a);
void picard_error_feedback_GRGM1200b(double* X, double* del_X, double* del_a);
#endif
