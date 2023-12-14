/**
 * @file TwoBody.cpp
 * @author David Stanley (davidms4@illinois.edu)
 * @brief Functions for calculating astrodynamic quantities for a two-body orbit.
 * @version 0.1
 * @date 2023-12-06
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include <cmath>
#include <vector>
#include "TwoBody.h"


// Function to calculate E using Newton's method
double calculateE(double M, double e, double tol = 1e-16)
{
    double E = M;
    double g = 1;
    int iters = 0;
    int max_iter = 100000;

    while (std::abs(g) > tol && iters < max_iter)
    {
        g = E - e * std::sin(E) - M;
        double dgde = 1 - e * std::cos(E);
        E = E - g / dgde;
        iters++;
    }

    return E;
}

// Function to convert orbital elements to position and velocity vectors
std::vector<std::vector<double>> elms2rv(double a, double e, double inc, double Om, double w, double M, double mu)
{
    double tol = 1e-10;

    // Calculate E
    double E = calculateE(M, e, tol);

    // Calculate true anomaly
    double f = 2 * std::atan(std::sqrt((1 + e) / (1 - e)) * std::tan(E / 2));

    // Calculate radius and speed
    double h = std::sqrt(mu * a * (1 - e * e));
    double r = a * (1 - e * std::cos(E));
    double T = w + f;

    // Calculate position vector in ECI coordinates
    std::vector<double> r_vec = {
        r * (std::cos(Om) * std::cos(T) - std::sin(Om) * std::sin(T) * std::cos(inc)),
        r * (std::sin(Om) * std::cos(T) + std::cos(Om) * std::sin(T) * std::cos(inc)),
        r * (std::sin(T) * std::sin(inc))};

    // Calculate velocity vector
    std::vector<double> v_vec = {
        mu / h * (-(std::cos(Om) * (std::sin(T) + e * std::sin(w))) - std::sin(Om) * (std::cos(T) + e * std::cos(w)) * std::cos(inc)),
        mu / h * (-(std::sin(Om) * (std::sin(T) + e * std::sin(w))) + std::cos(Om) * (std::cos(T) + e * std::cos(w)) * std::cos(inc)),
        mu / h * ((std::cos(T) + e * std::cos(w)) * std::sin(inc))};

    return {r_vec, v_vec};
}
