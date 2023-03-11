/**
 * @file Gravity.h
 * @author David Stanley (davidms4@illinois.edu)
 * @brief Functions for choosing gravity models
 * @version 0.1
 * @date 2023-03-10
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef __GRAVITY__
#define __GRAVITY__

#include "EGM2008.h"
#include "GRGM1200b.h"
#include "Orbit.h"
using namespace std;

void Gravity(double* p, double* Gxyz, int deg, Orbit orbit);



#endif