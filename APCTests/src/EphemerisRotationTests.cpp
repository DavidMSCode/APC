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
#include "Orbit.h"
#include "EphemerisRotationTests.h"
#include "EphemerisRotation.h"
#include "SpiceUsr.h"
using namespace std;

int main()
{
    int failuresRot = EphemerisRotationTest();
    return failuresRot;
}

int EphemerisRotationTest()
{
    double tol = 1e-12;
    int failures=0;
    // Load spice kernels
    string spkFile = "de440.bsp";
    string lskFile = "naif0012.tls";
    furnsh_c(spkFile.c_str());
    furnsh_c(lskFile.c_str());
    furnsh_c("moon_de440_220930.tf");
    furnsh_c("moon_pa_de440_200625.bpc");
    furnsh_c("pck00010.tpc");
    Orbit orbit("MOON", "MOON_PA", "J2000");
    vector<double> r0 = {1838., 0., 0.};              // Initial Position (km)
    vector<double> v0 = {0., 0.74844211, 1.49688755}; // Initial Velocity (km/s)

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
    InertialToBodyFixed(xI,vI,xF,vF,3000,orbit);
    for(int i=0;i<3;i++){
        r0[i]-=xF[i];
        v0[i]-=vF[i];
        if (abs(r0[i])>tol){
            failures+=1;
        }
        if (abs(v0[i])>tol){
            failures+=1;
        }
    }
    unload_c(spkFile.c_str());
    unload_c(lskFile.c_str());
    unload_c("moon_de440_220930.tf");
    unload_c("moon_pa_de440_200625.bpc");
    unload_c("pck00010.tpc");
    return failures;
}