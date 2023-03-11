/**
 * @file Gravity.cpp
 * @author David Stanley (davidms4@illinois.edu)
 * @brief Chooses gravity model based on primary body
 * @version 0.1
 * @date 2023-03-10
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include <string>
#include "Gravity.h"
using namespace std;

void Gravity(double *p, double *Gxyz, int deg, Orbit orbit)
{
    //Get the name of the primary body
    string primary = orbit.GetPrimaryBody();
    //convert primary body string to lowercase for consistency
    transform(primary.begin(),primary.end(),primary.begin(),[](unsigned char c){ return std::tolower(c);});
    //Find the right gravity model
    if(primary.compare("earth")==0){
        EGM2008(p,Gxyz,deg);
    }
    else if(primary.compare("moon")==0){
        GRGM1200b(p,Gxyz,deg);
    }
    else{
        //No model found return 0 gravitational acceleration
        for(int i=0;i<3;i++){
            Gxyz[i] = 0.0;
        }
    }
}