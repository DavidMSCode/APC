/*
*  AUTHORS:          David Stanley (davidms4@illinois.edu)
*  DATE WRITTEN:     Feb 2022
*  LAST MODIFIED:    Feb 2022
*  AFFILIATION:      Department of Aerospace Engineering, University of Illinois Urbana-Champaign
*  DESCRIPTION:      Class that stores orbit solution and properties
*  REFERENCE:        Woollands, R., and Junkins, J., "Nonlinear Differential Equation Solvers
*                    via Adaptive Picard-Chebyshev Iteration: Applications in Astrodynamics", JGCD, 2016.
*
*/

#include <string>
#include <vector>
#include <array>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include "Orbit.h"
#include "const.h"

using namespace std;
//Orbit Constructors
Orbit::Orbit(){
//Empty Orbit Constructor
}
Orbit::Orbit(string primary, string frame){
    //Initialize orbit with primary body and input data frame
    if(InSet(primary,C_VALID_PRIMARIES)){
        _primary = primary;
        _validPrimary = true;
        SetMu(_primary);
    }
}
Orbit::Orbit(std::vector<std::vector<double > > Solution){
    SetSolution(Solution);
}
Orbit::Orbit(double area, double reflectivity, double mass, double Cd, bool compute_drag, bool compute_SRP, bool compute_third_body,bool compute_hamiltonian, int id){
    //Set primary body for two body targets
    _primary = "Earth";     //defaults to orbiting around Earth
    _validPrimary = true;
    //Set orbital properties
    SetProperties(area, reflectivity, mass, Cd, compute_drag, compute_SRP, compute_third_body, compute_hamiltonian, id);
    suborbital = false; 
}

Orbit::Orbit(double area, double reflectivity, double mass, double Cd, bool compute_drag, bool compute_SRP, bool compute_third_body,bool compute_hamiltonian, int id, string primary = "Earth"){
    // Set primary body for two body orbits after checking if valid
    if(InSet(primary,C_VALID_PRIMARIES)){
        _primary = primary;
        _validPrimary = true;
        SetMu(_primary);
    }
    
    SetProperties(area, reflectivity, mass, Cd, compute_drag, compute_SRP, compute_third_body, compute_hamiltonian, id);
    suborbital = false;
}

//Getter Functions
std::vector<std::vector<double> > Orbit::getPosition(){
    std::vector<std::vector<double> >::const_iterator first = Soln.begin() + 1;
    std::vector<std::vector<double> >::const_iterator last = Soln.begin() + 4;
    std::vector<std::vector<double> > Position(first, last);

    return Position;
}
std::vector<std::vector<double> > Orbit::getVelocity(){
    std::vector<std::vector<double> >::const_iterator first = Soln.begin() + 4;
    std::vector<std::vector<double> >::const_iterator last = Soln.begin() + 7;
    std::vector<std::vector<double> > Velocity(first, last);

    return Velocity;
}

//Setter Functions
void Orbit::SetSolution(std::vector<std::vector<double > > Solution){
    //store the solution
    Soln = Solution;
}

void Orbit::SetTimeVec(std::vector<double> time_vec){
    //Store user time vector and raise flag indicating that a time vector was given
    T = time_vec;
    USER_TIME = true;
}

void Orbit::SetProperties(double area, double reflectance, double mass, double Cd, bool compute_drag, bool compute_SRP, bool compute_third_body, bool compute_hamiltonian, int id, string primary){
    //Set satellite properties and flags for perturbations
    satproperties._Area = area;
    satproperties._Mass = mass;
    satproperties._Reflectance = reflectance;
    satproperties._Cd = Cd;
    //Perturbation flags
    Compute_Drag = compute_drag;
    Compute_SRP = compute_SRP;
    Compute_Third_Body = compute_third_body;
    //Other flags
    Compute_Hamiltonian = compute_hamiltonian;
    //ID
    ID = id;
}

void Orbit::SetSubOrbital(){
    suborbital = true;
}

void Orbit::SetCC(std::vector<double> A, std::vector<double> B, std::vector<double> W1, std::vector<double> W2, int N, int coeff_size, std::vector<double> seg_times, double TF, double T0,int total_segs){
    //Set values for chebyshev coefficients and associated variables
    CC.A = A;
    CC.B = B;
    CC.W1 = W1;
    CC.W2 = W2;
    CC.N = N;
    CC.coeff_size = coeff_size;
    CC.seg_times = seg_times;
    CC.TF = TF;
    CC.T0 = T0;
    CC.total_segs = total_segs;
}

void Orbit::SetMu(string primary){
    //lower primary name case
    std::transform(primary.begin(),primary.end(),primary.begin(),
    [](unsigned char c){ return std::tolower(c);});
    //Assign mu based on name of primary body. Defined in const.h
    if(primary.compare("earth")==0){
        _mu = C_MU_EARTH;
    }
    else if (primary.compare("moon")==0)
    {
        _mu = C_MU_MOON;
    }
}

bool Orbit::InSet(string item, vector<string> validset){
    bool inSet = false;
    //convert item to lower case for consistency
    std::transform(item.begin(), item.end(), item.begin(),
    [](unsigned char c){ return std::tolower(c);});
    //Iterate through the valid set
    for(string validItem : validset){
        string lowerValid = validItem;
        //convert valid item to lowercase for consistency
        std::transform(lowerValid.begin(),lowerValid.end(),lowerValid.begin()
        ,[](unsigned char c){ return std::tolower(c);});
        //Check if item is the same as the validItems
        if(item.compare(lowerValid)==0){
            //Found in valid set, end the loop
            inSet = true;
            break;
        }
    }
    return inSet;
}

