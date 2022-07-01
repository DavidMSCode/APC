/*
*  AUTHORS:          David Stanley (davidms4@illinois.edu)
*  DATE WRITTEN:     Feb 2022
*  LAST MODIFIED:    Feb 2022
*  AFFILIATION:      Department of Aerospace Engineering, University of Illinois Urbana-Champaign
*  DESCRIPTION:      Methods for drag, solar radiation pressure and third body perturbations
   REFERENCE:       1. Woollands, R., and Junkins, J., "Nonlinear Differential Equation Solvers
                    via Adaptive Picard-Chebyshev Iteration: Applications in Astrodynamics", JGCD, 2016
                    2. Vallado, David A, and Wayne D McClain. Fundamentals of Astrodynamics and 
                    Applications. 3rd ed., Springer, 2007.
*/
#ifndef __PERTURBATIONS__
#define __PERTURBATIONS__

#include <Orbit.h>
#include "Ephemeris.hpp"

void Perturbed_SRP(double time, double* X, Orbit orb, EphemerisManager ephem, double* SRP_aECI);

void Perturbed_Drag(double* X, double* V, Orbit orb, double* drag_aECEF);

void Perturbed_three_body(double time, double* X, Orbit orb, EphemerisManager ephem, double* third_body_aECI);

int LastFirstSearch(double *p, int length_t, double key);

double atmospheric_density(double alt);

#endif