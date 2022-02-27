/*
*  AUTHORS:          David Stanley (davidms4@illinois.edu)
*  DATE WRITTEN:     Feb 2022
*  LAST MODIFIED:    Feb 2022
*  AFFILIATION:      Department of Aerospace Engineering, University of Illinois Urbana-Champaign, IL
*  DESCRIPTION:      small methods to test linking against cspice library
*  INPUT:
*  OUTPUT:
*  COMMENTS:
*/

#include "SpiceUsr.h"
#include <iostream>
#include <linktest.h>

void Linktest(){
    furnsh_c("/Users/davidstanley/Documents/Github/adaptive_picard_chebyshev_C/adaptive_picard_chebyshev_v4/bin/de430.bsp");
    nestLinkMethod();
    kclear_c();
    std::cout << "finished linktest\n";
}

void nestLinkMethod(){
    SpiceDouble     et;
    SpiceDouble     lt;
    SpiceDouble     state [6];
    spkezr_c( "SUN", 0, "J2000", "LT+S", "EARTH", state, &lt);
    std::cout << "The sun is at position (" << state[0]<<", "<<state[1]<<", "<<state[2]<<")\n";
}