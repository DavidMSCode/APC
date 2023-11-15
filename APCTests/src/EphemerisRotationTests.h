/**
 * @file EphemerisRotationTests.h
 * @author David Stanley (davidms4@illinois.edu)
 * @brief 
 * @version 0.1
 * @date 2023-03-26
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef __ROTATIONTEST__
#define __ROTATIONTEST__

int EphemerisRotationTests(int argc, char** argv);

int EphemerisRotationTest(vector<double> r0, vector<double> v0, double rtol=1e-15);

#endif