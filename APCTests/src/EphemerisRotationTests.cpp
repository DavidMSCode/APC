/**
 * @file EphemerisRotationTests.cpp
 * @author David Stanley (davidms4@illinois.edu)
 * @brief
 * @version 0.1
 * @date 2023-03-26
 *
 * @copyright Copyright (c) 2023
 *
 */

#include <vector>
#include <iostream>
#include "Orbit.h"
#include "EphemerisRotationTests.h"
#include "EphemerisRotation.h"
#include "SpiceUsr.h"
using namespace std;

int EphemerisRotationTests(int argc, char **argv)
{
    vector<double> r0 = {1838., 0., 0.};              // Initial Position (km)
    vector<double> v0 = {0., 0.74844211, 1.49688755}; // Initial Velocity (km/s)
    int failuresRot = EphemerisRotationTest(r0,v0);

    //test using numbers from Mark H
    r0 = {1.368858356478835e6, 1.226051809239666e6, 3.560563099011000e6};              // Initial Position (km)
    v0 = {0., 0.74844211, 1.49688755}; // Initial Velocity (km/s)
    failuresRot += EphemerisRotationTest(r0,v0);
    return failuresRot;
}

int EphemerisRotationTest(vector<double> r0, vector<double> v0, double rtol/*=1e-15*/)
{
       // Load spice kernels
    string spkFile = "de440.bsp";
    string lskFile = "naif0012.tls";
    furnsh_c(spkFile.c_str());
    furnsh_c(lskFile.c_str());
    furnsh_c("moon_de440_220930.tf");
    furnsh_c("moon_pa_de440_200625.bpc");
    furnsh_c("pck00010.tpc");
    furnsh_c("naif0012.tls");
    string date = "1969-07-22T00:00:00.000000"; // Date of interest
    double et;
    str2et_c(date.c_str(), &et);

    int failures = 0;
 
    Orbit orbit("MOON", "MOON_PA", "J2000");


    double T = 7690.61; // Orbital period (s)
    double t0 = 0;      // initial time (s)
    double tf = T;

    orbit.SetPosition0(r0);
    orbit.SetVelocity0(v0);
    orbit.SetIntegrationTime(3000, tf);
    orbit.SetComputeThirdBody();
    orbit.SetComputeHamiltonian();
    vector<double> xI;
    vector<double> vI;
    BodyFixedToInertial(r0, v0, xI, vI, 3000, orbit);
    vector<double> xF;
    vector<double> vF;
    InertialToBodyFixed(xI, vI, xF, vF, 3000, orbit);
    vector<double> r0_error(3);
    vector<double> v0_error(3);
    for (int i = 0; i < 3; i++)
    {
        r0_error[i]=r0[i] - xF[i];
        v0_error[i]=v0[i] - vF[i];
        // if (abs(r0[i])>tol){
        //     failures+=1;
        // }
        // if (abs(v0[i])>tol){
        //     failures+=1;
        // }
    }
    //calculate the error magnitude
    double abs_r0mag = sqrt(r0_error[0] * r0_error[0] + r0_error[1] * r0_error[1] + r0_error[2] * r0_error[2]);
    double abs_v0mag = sqrt(v0_error[0] * v0_error[0] + v0_error[1] * v0_error[1] + v0_error[2] * v0_error[2]);
    //calculate the relative error
    double r0mag = abs_r0mag / sqrt(r0[0] * r0[0] + r0[1] * r0[1] + r0[2] * r0[2]);
    double v0mag = abs_v0mag / sqrt(v0[0] * v0[0] + v0[1] * v0[1] + v0[2] * v0[2]);
    // check that the magnitude of each difference vector is less than the tolerance
    if (r0mag > rtol)
    {
        failures += 1;
    }
    if (v0mag > rtol)
    {
        failures += 1;
    }
    // print the two error vectors their absolute magnitude and relative error
    cout << "r0_error: " << r0_error[0] << " " << r0_error[1] << " " << r0_error[2] << endl;
    cout << "v0_error: " << v0_error[0] << " " << v0_error[1] << " " << v0_error[2] << endl;
    cout << "abs_r0mag: " << abs_r0mag << endl;
    cout << "abs_v0mag: " << abs_v0mag << endl;
    cout << "relative error r: " << r0mag << endl;
    cout << "relative error v: " << v0mag << endl;

    // Unload spice kernels

    unload_c(spkFile.c_str());
    unload_c(lskFile.c_str());
    unload_c("moon_de440_220930.tf");
    unload_c("moon_pa_de440_200625.bpc");
    unload_c("pck00010.tpc");
    return failures;
}