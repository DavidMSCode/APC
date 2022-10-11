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

#include "Orbit.h"
#include "Ephemeris.hpp"

/**
 * @brief Returns the acceleration perturbation due to solar radiation pressure
 * 
 * @param time time since epoch (s)
 * @param X position in ECI coordinates (km)
 * @param orb orbital object which stores parameters
 * @param ephem ephemeris caching manager that returns the Sun's position
 * @param SRP_aECI The SRP acceleration in ECI (km/s^2)
 */
void Perturbed_SRP(double time, double* X, Orbit orb, EphemerisManager ephem, double* SRP_aECI);

/**
 * @brief Returns the acceleration perturbation due to atmospheric drag using Vallado's exponential atmospheric model
 * 
 * @param X position in ECEF coordinates (km)
 * @param V velocity in ECEF coordinates (km/s)
 * @param orb orbital object which contains necessary parameters
 * @param drag_aECEF The drag acceleration in ECEF (km/s^2)
 */
void Perturbed_Drag(double* X, double* V, Orbit orb, double* drag_aECEF);

/**
 * @brief Returns the perturbed third body acceleration due to the Moon and the Sun
 * 
 * @param time time since epoch (s)
 * @param X position in ECI coordinates (km)
 * @param orb orbital object which contains necessary parameters
 * @param ephem ephemeris caching manager that returns the Sun's and the Moon's positions
 * @param third_body_aECI the third body acceleration in ECI (km/s^2)
 */
void Perturbed_three_body(double time, double* X, Orbit orb, EphemerisManager ephem, double* third_body_aECI);

/**
 * @brief This function attempts to find the index for the value given by the key by searching from largest
     to smallest. If the key is not present the function returns the index of the first item smaller 
     than the key value assuming a sorted list. Used to search the atmospheric density tables.
 * 
 * @param p Table pointer
 * @param length_t Length of table
 * @param key Value to find in the table
 * @return int 
 */
int LastFirstSearch(double *p, int length_t, double key);

/**
 * @brief  Returns the qtmospheric density utilizing Vallado's exponential atmospheric model given an altitude above sea-level
 * 
 * @param alt Altitude of orbiting object (km)
 * @return The atmospheric density (kg/m^3)
 */
double atmospheric_density(double alt);

#endif