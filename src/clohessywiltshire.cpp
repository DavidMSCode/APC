 // File: clohessywiltshire.cpp
// Author: David Stanley
// Date: September 15
// Description: Implements analytic relative motion for spacecraft in a circular orbit around a central body

#include "clohessywiltshire.h"
#include "c_functions.h"
#include <cmath>
#include <vector>

using namespace std;

/// @brief 
/// @param t    Time since initial state
/// @param ds0  Initial state vector [x, y, z, vx, vy, vz]
/// @param n    Mean motion of the central body
/// @return     New state vector after time t
vector<double> ClohessyWiltshire(double t, vector<double> ds0, double n)
{
    vector<double> ds(6);
    ds[0] = 4 * ds0[0] + 2 / n * ds0[4] + ds0[3] / n * sin(n * t) - (3 * ds0[0] + 2 / n * ds0[4]) * cos(n * t);
    ds[1] = ds0[1] - 2 / n * ds0[3] - 3 * (2 * n * ds0[0] + ds0[4]) * t + 2 * (3 * ds0[0] + 2 / n * ds0[0] + 2 / n * ds0[4]) * sin(n * t) + 2 / n * ds0[3] * cos(n * t);
    ds[2] = 1 / n * ds0[2] * sin(n * t) + ds0[2] * cos(n * t);
    ds[3] = 3 * n * sin(n * t) * ds0[0] + cos(n * t) * ds0[3] + 2 * sin(n * t) * ds0[4];
    ds[4] = 6 * n * (cos(n * t) - 1) * ds0[0] - 2 * sin(n * t) * ds0[3] + (4 * cos(n * t) - 3) * ds0[4];
    ds[5] = -n * sin(n * t) * ds0[2] + cos(n * t) * ds0[5];
    return ds;
} // end ClohessyWiltshire

pair<vector<double>, vector<double>> orb2lvlh(vector<double> r,vector<double> v, vector<double> dr, vector<double> dv){
//convert dr dv in the eci frame to the lvlh frame defined by r and h

    vector<double> h = {r[1]*v[2] - r[2]*v[1], r[2]*v[0] - r[0]*v[2], r[0]*v[1] - r[1]*v[0]}; // angular momentum
    double h_norm = sqrt(h[0]*h[0] + h[1]*h[1] + h[2]*h[2]);
    vector<double> h_hat = {h[0]/h_norm, h[1]/h_norm, h[2]/h_norm}; // unit vector in the direction of h

    vector<double> r_hat = {r[0]/sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]), r[1]/sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]), r[2]/sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2])}; // unit vector in the direction of r

    //v_hat = h_hatxr_hat
    vector<double> v_hat = {h_hat[1]*r_hat[2] - h_hat[2]*r_hat[1], h_hat[2]*r_hat[0] - h_hat[0]*r_hat[2], h_hat[0]*r_hat[1] - h_hat[1]*r_hat[0]}; // unit vector in the direction of v

    vector<double> dr_lvlh = {Cdot(r_hat, dr), Cdot(v_hat, dr), Cdot(h_hat, dr)}; // dr in lvlh frame
    vector<double> dv_lvlh = {Cdot(r_hat, dv), Cdot(v_hat, dv), Cdot(h_hat, dv)}; // dv in lvlh frame
    return {dr_lvlh, dv_lvlh};
}

pair<vector<double>, vector<double>> lvlh2orb(vector<double> r,vector<double> v, vector<double> dr, vector<double> dv){
//convert dr dv in the elvlh frame defined by r and h to the eci frame 

    vector<double> h = {r[1]*v[2] - r[2]*v[1], r[2]*v[0] - r[0]*v[2], r[0]*v[1] - r[1]*v[0]}; // angular momentum
    double h_norm = sqrt(h[0]*h[0] + h[1]*h[1] + h[2]*h[2]);
    vector<double> h_hat = {h[0]/h_norm, h[1]/h_norm, h[2]/h_norm}; // unit vector in the direction of h

    vector<double> r_hat = {r[0]/sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]), r[1]/sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]), r[2]/sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2])}; // unit vector in the direction of r

    //v_hat = h_hatxr_hat
    vector<double> v_hat = {h_hat[1]*r_hat[2] - h_hat[2]*r_hat[1], h_hat[2]*r_hat[0] - h_hat[0]*r_hat[2], h_hat[0]*r_hat[1] - h_hat[1]*r_hat[0]}; // unit vector in the direction of v

    vector<double> row1 = {r_hat[0], v_hat[0], h_hat[0]};
    vector<double> row2 = {r_hat[1], v_hat[1], h_hat[1]};
    vector<double> row3 = {r_hat[2], v_hat[2], h_hat[2]};
    vector<double> dr_eci = {Cdot(row1, dr), Cdot(row2, dr), Cdot(row3, dr)}; // dr in lvlh frame
    vector<double> dv_eci = {Cdot(row1, dv), Cdot(row2, dv), Cdot(row3, dv)}; // dv in lvlh frame
    return {dr_eci, dv_eci};

}