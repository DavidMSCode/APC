/* const.h
*  AUTHOR:           Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     May 2016
*  LAST MODIFIED:    Feb 2021
*  AFFILIATION:      Department of Aerospace Engineering, Texas A & M University, College Station, TX
*  DESCRIPTION:      Constants
*/

#include <string>
#include <vector>
using namespace std;

#ifndef _CONSTANTS_
#define _CONSTANTS_
#define C_PI 3.1415926535897932      // Pi
#define C_MU_EARTH 3.9860043543609598E+05          // Gravitational Constant Earth [km^3/s^2]
#define C_MU C_MU_EARTH          // Gravitational Constant of the Earth [km^3/s^2]
#define C_MU_MOON 4.9028000661637961E+03
#define C_MU_SUN 1.3271244004193938E+11
#define C_MUSun 1.32712440018e11     // Gravitational Constant of the Sun [km^3/s^2]
#define C_MUMoon 4.9048695e3         //Gravitational Constant of the Moon [km^3/s^2]
#define C_MU_EARTH_CAN 1                    // Gravitational Constant Canonical Units
#define C_MU_MOON_CAN C_MU_MOON/C_MU_EARTH
#define C_MU_SUN_CAN C_MU_SUN/C_MU_EARTH
// #define C_omega 7292115.0e-011       // Angular Speed of Earth [rad/s]
// #define C_omega 2.661699624635926e-06   //angular speed of Moon [rad/s]
#define C_omega 0.0
#define C_Req 6378.137               // Equatorial Radius of Earth [km]
#define C_Rmoon 1738.0               // Reference Radius of the Moon [km] from GRAIL model
#define C_ckm 299792.458             // Speed of Light (km/s)
#define C_Gsc 1362                   // Solar constant (W/m^2 or kg/s^3)
#define g0 9.8065                    // Earth's Surface Gravity (m/s^2)
#define C_J2 0.00108263              // Earth's Second Zonal Harmonic
#define RS 696000.0                  // Solar Radius [km]
#define DU C_Req
#define TU sqrt(pow(DU,3.0)/C_MU_EARTH)
// #define DU C_Rmoon
// #define TU sqrt(pow(DU,3.0)/C_MU_MOON)
#define VU DU/TU
#define si2can pow(TU,2.0)/(DU*1000.0) // Canonical unit conversion factor
#define JD2000 2451545.0             // Julian days of J2000
const vector<string> C_VALID_PRIMARIES = {"earth","moon"};  //List of valid primary body names
#endif