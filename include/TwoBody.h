/**
 * @file TwoBody.h
 * @author David Stanley (davidms4@illinois.edu)
 * @brief 
 * @version 0.1
 * @date 2023-12-06
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef TWOBODY_H
#define TWOBODY_H

double calculateE(double M, double e, double tol);

std::vector<std::vector<double>> elms2rv(double a, double e, double inc, double Om, double w, double M, double mu);

#endif // !TWOBODY_H