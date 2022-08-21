/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     June 2022
*  LAST MODIFIED:    June 2022
*  AFFILIATION:      Department of Aerospace Engineering, University of Illinois, Urbana-Champaign
*  DESCRIPTION:      Compute Modified Equinocitial Element Dynamics
*
* INPUT:
*    MEE     -- Modified Equinoctial Elements
*    a_lvlh  -- Perturbing acceleration in the LVLH frame (km/s^2)
*
* OUTPUTS:
*    MEE_dot -- Modified Equinoctial Element dynamics
*
*/

#include "const.h"
#include "mee_dot.h"
#include "c_functions.h"
#include <math.h>
#include <vector>

void mee_dot(std::vector<double> &MEE, std::vector<double> &a_lvlh, std::vector<double> &MEE_dot, int M){

    double cosL = 0.0;
    double sinL = 0.0;
    double sqpm = 0.0;
    double q = 0.0;
    double s_sq = 0.0;
    
    for (int i=1; i<=M+1; i++){
    
        // MEE Equations
        cosL = cos(MEE[ID2(i,6,M+1)]);
        sinL = sin(MEE[ID2(i,6,M+1)]);
        sqpm = sqrt(MEE[ID2(i,1,M+1)]/C_MU);
        q    = 1.0 + MEE[ID2(i,2,M+1)]*cosL + MEE[ID2(i,3,M+1)]*sinL;
        s_sq = 1.0 + pow(MEE[ID2(i,4,M+1)],2) + pow(MEE[ID2(i,5,M+1)],2);
        MEE_dot[ID2(i,1,M+1)] = sqpm*2.0*MEE[ID2(i,1,M+1)]*a_lvlh[ID2(i,2,M+1)]/q;
        MEE_dot[ID2(i,2,M+1)] = sqpm*(a_lvlh[ID2(i,1,M+1)]*sinL + (((q+1.0)*cosL + MEE[ID2(i,2,M+1)])*a_lvlh[ID2(i,2,M+1)])/q - (MEE[ID2(i,3,M+1)]*(MEE[ID2(i,4,M+1)]*sinL - MEE[ID2(i,5,M+1)]*cosL)*a_lvlh[ID2(i,3,M+1)])/q);
        MEE_dot[ID2(i,3,M+1)] = sqpm*(-a_lvlh[ID2(i,1,M+1)]*cosL + (((q+1.0)*sinL + MEE[ID2(i,3,M+1)])*a_lvlh[ID2(i,2,M+1)])/q + (MEE[ID2(i,2,M+1)]*(MEE[ID2(i,4,M+1)]*sinL - MEE[ID2(i,5,M+1)]*cosL)*a_lvlh[ID2(i,3,M+1)])/q);
        MEE_dot[ID2(i,4,M+1)] = sqpm*s_sq/2.0/q*cosL*a_lvlh[ID2(i,3,M+1)];
        MEE_dot[ID2(i,5,M+1)] = sqpm*s_sq/2.0/q*sinL*a_lvlh[ID2(i,3,M+1)];
        MEE_dot[ID2(i,6,M+1)] = (sqrt(C_MU*MEE[ID2(i,1,M+1)])*pow(q/MEE[ID2(i,1,M+1)],2) + sqpm*(MEE[ID2(i,4,M+1)]*sinL - MEE[ID2(i,5,M+1)]*cosL)/q*a_lvlh[ID2(i,3,M+1)]);
        
    };
    
}