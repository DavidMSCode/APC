/**
 * @file FandG.cpp
 * @author Robyn Woollands (rmw@illinois.edu)
 * @brief 
 * @version 0.1
 * @date 2022-06-30
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <math.h>

#include "FandG.h"
#include "const.h"



//Define F and G magic numbers.
const int newtonMaxIt = 300;
const double newtonTol = 1E-13;

void FandG( const double* z0, double* zf, const double dt, const double mu = C_MU)
{
    int kk;

    double rMag = sqrt( z0[0]*z0[0] + z0[1]*z0[1] + z0[2]*z0[2] );
    double vSq  = z0[3]*z0[3] + z0[4]*z0[4] + z0[5]*z0[5];
    double a    = 1.0 / ( 2.0/rMag - vSq/mu );
    double sig0 = (z0[0]*z0[3] + z0[1]*z0[4] + z0[2]*z0[5])/sqrt(mu);
    double EHat = newtonFandG( a, dt, rMag, sig0, newtonTol, mu);
    double r    = a + (rMag-a)*cos(EHat) + sqrt(a)*sig0*sin(EHat);
    double F    = 1.0 - a/rMag*( 1 - cos(EHat) );
    double G    = dt + sqrt( a*a*a/mu )*( sin(EHat) - EHat )  ;
    double Fdot = -sqrt(a*mu)*sin(EHat)/(r*rMag);
    double Gdot = 1.0 - a/r*(1.0-cos(EHat));

    // Current state is a linear combination of the initial condition vectors
    int ii;
    for ( ii=0; ii<3; ii++ ) {
        zf[ii]   = F*z0[ii]    + G*z0[ii+3];
        zf[ii+3] = Fdot*z0[ii] + Gdot*z0[ii+3];
    }
}

double newtonFandG( const double a, const double dt, const double rMag,
                    const double sig0, const double tol, const double mu)
{

    // Initial guess
    //double EHat = C_PI;
    double EHat = C_PI;
    int ctr = 0;
    double fx, dfx;
    double dE = 1.0;

    while ( ( fabs(dE/EHat) > tol ) && ( ctr <= newtonMaxIt ) ) {

        // Calculate dE correction
        fx = sqrt(mu)/(a*sqrt(a))*dt  - EHat + (1.0-rMag/a)*sin(EHat) +  \
            sig0/sqrt(a)*(cos(EHat)-1.0);
        dfx = -1.0 + ((1.0-rMag/a)*cos(EHat) - sig0/sqrt(a)*sin(EHat));
        dE  = -1.0*fx/dfx;
        if ( ctr == newtonMaxIt ) {
        }

        // Add correction and increment ctr
        EHat += dE;
        ctr++;

    }

    return EHat;
}
