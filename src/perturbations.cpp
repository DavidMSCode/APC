/*
*  AUTHORS:          David Stanley (davidms4@illinois.edu)
*  DATE WRITTEN:     Feb 2022
*  LAST MODIFIED:    Feb 2022
*  AFFILIATION:      Department of Aerospace Engineering, University of Illinois Urbana-Champaign
*  DESCRIPTION:      Methods for drag, solar radiation pressure and third body perturbations
*  REFERENCE:        Woollands, R., and Junkins, J., "Nonlinear Differential Equation Solvers
*                    via Adaptive Picard-Chebyshev Iteration: Applications in Astrodynamics", JGCD, 2016.
*/

#include "SpiceUsr.h"
#include <c_functions.h>
#include <math.h>
#include <string>
#include <vector>
#include <Orbit.h>
#include <const.h>
#include <omp.h>
#include "Ephemeris.hpp"


int LastFirstSearch(double *p, int length_t, double key){
    /*This function attempts to find the index for the value given by the key by searching from largest
     to smallest. If the key is not present the function returns the index of the first item smaller 
     than the key value assuming a sorted list. A binary search would be more efficient, 
     but the arrays are very short.
    */
    //Set flags and index for iterating through vector
    bool indexFound=false;  
    int i = length_t-1;
    int index_sol=-1;

    while(!indexFound && i>=0){
        if (*p <= key){
            //printf("%f<=%f\n",*p,key);
            //Smaller or equal value found. End loop.
            index_sol=i;
            indexFound=true;
        }
        //increment search index and array pointer
        i--;
        p--;
    }
    //Return found index. If value is -1 then an index was not found.
    return index_sol;
}

double atmospheric_density(double alt){
    // Returns Atmospheric density  utilizing Vallado's exponential atmospheric model given an altitude above sea-level
    double density_vector[28] = {
    1.225,
    3.899e-2,
    1.774e-2,
    3.972e-3,
    1.057e-3,
    3.206e-4,
    8.770e-5,
    1.905e-5,
    3.396e-6,
    5.297e-7,
    9.661e-8,
    2.438e-8,
    8.484e-9,
    3.845e-9,
    2.070e-9,
    5.464e-10,
    2.789e-10,
    7.248e-11,
    2.418e-11,
    9.518e-12,
    3.725e-12,
    1.585e-12,
    6.967e-13,
    1.454e-13,
    3.614e-14,
    1.170e-14,
    5.245e-15,
    3.019e-15};

    double base_alt_vector[28] = {
    0,
    25,
    30,
    40,
    50,
    60,
    70,
    80,
    90,
    100,
    110,
    120,
    130,
    140,
    150,
    180,
    200,
    250,
    300,
    350,
    400,
    450,
    500,
    600,
    700,
    800,
    900,
    1000};
  
    double scale_height_vector[28] = {
    7.249,
    6.349,
    6.682,
    7.554,
    8.382,
    7.714,
    6.549,
    5.799,
    5.382,
    5.877,
    7.263,
    9.473,
    12.636,
    16.149,
    22.523,
    29.740,
    37.105,
    45.546,
    53.628,
    53.298,
    58.515,
    60.828,
    63.822,
    71.835,
    88.667,
    124.64,
    181.05,
    268.00};

    //Find the index for the base altitude and get the corresponding base altitude, scaling height and atmospheric density
    int length_t = sizeof(base_alt_vector)/sizeof(base_alt_vector[0]);                   //Get length of array
    double *p = (double *)(&base_alt_vector + 1) - 1;                                    //Get pointer for last element in array
    int index = LastFirstSearch(p, length_t, alt);                                       //Find base altitude index
    double density = 0.0;
    if (index!=-1){
        //Get all values at index
        double base_density = density_vector[index];
        double base_alt = base_alt_vector[index];
        double scale_height = scale_height_vector[index];
        //calculate atmospheric density at given altitude
        density = base_density*exp(-(alt-base_alt)/scale_height);
    }   
    return density;
}

void Perturbed_SRP(double time, double* X, Orbit orb, EphemerisManager ephem, double* SRP_aECI){
    double satvec[3];
    double sunvec[3];
    double satsununitvec[3];
    double satsunvec[3];
    double norm_sunpos;
    double norm_satpos;
    double norm_satsunpos;
    std::vector<double> sunstate;

    for(int i=0;i<3;i++){
        SRP_aECI[i]=0.0;
        satvec[i] = X[i];
    }
    if (orb.Compute_SRP){
        double Mass = orb.GetMass();                            //Sat mass (kg)
        double Area_m = orb.GetArea();                          //Sat cannonball area (m^2)
        double Area_km = Area_m/pow(1000,2);                    //Sat cannonball area (km^2)
        double Cr = orb.GetReflectance();                        //Coeffecicient of reflectance
        double r_eq = C_Req;                                    //Equatorial radius (km)
        double G_sc = C_Gsc;                                    //Solar constant (kg/s^3))
        double C = C_ckm;                                       //Speed of light (km/s)

        //Get Earth to Sun vector
        // ConstSpiceChar target[4] = "Sun";
        // ConstSpiceChar observer[6] = "Earth";
        // ConstSpiceChar iframe[6] = "J2000";
        double epoch = time;
        sunstate = ephem.getState("SUN",epoch);
        for(int i=0;i<3;i++){
            sunvec[i]=sunstate[i];
            satsunvec[i] = sunvec[i]-satvec[i];
        }
        //calculate vector lengths and sun invit vector
        Cnorm(sunvec,norm_sunpos);
        Cnorm(satvec,norm_satpos);
        Cnorm(satsunvec,norm_satsunpos);
        for(int i=0;i<3;i++){
            satsununitvec[i]=satsunvec[i]/norm_satsunpos;
        }
        //Use vector angles to find occlusion state
        double satsunangle;
        double angle1;
        double angle2;
        double sumangle;
        double p;
        satsunangle = acos(Cdot(sunvec,satvec)/(norm_satpos*norm_sunpos));
        angle1 = acos(r_eq/norm_satpos);
        angle2 = acos(r_eq/norm_sunpos);
        sumangle = angle1+angle2;
        if (sumangle>satsunangle){
            //calcualte srp acceleration in km/s^2
            p = G_sc/C*Cr*Area_km/Mass;
            for(int i=0;i<3;i++){
                SRP_aECI[i]=-p*satsununitvec[i];
            }
        }
    }
};

void Perturbed_Drag(double* X_ECEF, double* V_ECEF, Orbit orbit, double* drag_aECEF){
    /* Returns the atmospheric drag acceleration on the sattelite in a given orbit around Earth. 
    */
    for(int i=0;i<3;i++){
        drag_aECEF[i]=0.0;
    }
    if (orbit.Compute_Drag){
        double r = sqrt(pow(X_ECEF[0],2)+pow(X_ECEF[1],2)+pow(X_ECEF[2],2));   //distance of sattelite from center of the earth in km
        double s = sqrt(pow(V_ECEF[0],2)+pow(V_ECEF[1],2)+pow(V_ECEF[2],2));   //speed of sattelite in km/s

        double Mass = orbit.GetMass();                            //Sat mass (kg)
        double Area_m = orbit.GetArea();                          //Sat cannonball area (m^2)
        double Area_km = Area_m/pow(1000,2);                    //Sat cannonball area (km^2)
        double Cd = orbit.GetDragCoefficient();                   //Sat drag coefficient 
        double r_eq = C_Req;                                    //Equatorial radius (km)
        double alt = r-r_eq;                                    //Altitude (km)

        //Get the atmospheric density kg/km^3
        double rho = atmospheric_density(alt)*pow(1000,3);
        //Calculate drag acceleration km/s^2
        double drag_a = 0.5*rho*pow(s,2)*Area_km*Cd/Mass;

        for (int i=0;i<3;i++){
            //calc drag vector in km/s^2
            drag_aECEF[i]= -drag_a*V_ECEF[i]/s;
        }
    }
};

void Perturbed_three_body(double time, double* X, Orbit orb, EphemerisManager ephem, double* third_body_aECI){
    double satvec[3];
    double sunvec[3];
    double moonvec[3];
    double satsununitvec[3];
    double satmoonunitvec[3];
    double satsunvec[3];
    double satmoonvec[3];
    double norm_sunpos;
    //double norm_satpos;
    //double norm_moonpos
    double norm_satsunpos;
    double norm_satmoonpos;
    std::vector<double> sunstate;
    std::vector<double> moonstate;

    for(int i=0;i<3;i++){
        third_body_aECI[i]=0.0;
        satvec[i] = X[i];
    }
    if (orb.Compute_Third_Body){
        //Get Earth to Sun vector and Earth to moon vector
        double epoch = time;
        //Read ephem data at current time
        sunstate = ephem.getState("SUN",epoch);
        moonstate = ephem.getState("MOON",epoch);
        //Get vectors from satellite to third bodies
        for(int i=0;i<3;i++){
            sunvec[i]=sunstate[i];
            moonvec[i] = moonstate[i];
            satsunvec[i] = sunvec[i]-satvec[i];
            satmoonvec[i] = moonvec[i]-satvec[i];
        }
        //calculate vector lengths and sun/moon unit vector from sat
        Cnorm(sunvec,norm_sunpos);
        //Cnorm(satvec,norm_satpos);
        Cnorm(satsunvec,norm_satsunpos);
        Cnorm(satmoonvec,norm_satmoonpos);
        for(int i=0;i<3;i++){
            satsununitvec[i]=satsunvec[i]/norm_satsunpos;
            satmoonunitvec[i]=satmoonvec[i]/norm_satmoonpos;
        }
        //Use vector angles to find occlusion state
        double gSun = C_MUSun/pow(norm_satsunpos,2);
        double gMoon = C_MUMoon/pow(norm_satmoonpos,2);
        double q;
        double Fq;
        double tempvec[3];
        for(int i=0;i<3;i++){
            tempvec[i] = (2*sunvec[i]-satvec[i]);
        }
        q = Cdot(satvec,tempvec)/pow(norm_sunpos,2);
        Fq = q*(pow(q,2)-3*q+3)/(1+pow(1-q,1.5));
            //calcualte srp acceleration in km/s^2
        for(int i=0;i<3;i++){
            third_body_aECI[i]=gSun/norm_satsunpos*(Fq*sunvec[i]-satvec[i]) + gMoon*satmoonunitvec[i];
        }
    }
};





