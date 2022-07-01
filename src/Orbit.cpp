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


#include <vector>

#include "Orbit.h"

//Orbit Constructors
Orbit::Orbit(){

}
Orbit::Orbit(std::vector<std::vector<double > > Solution){
    SetSolution(Solution);
}
Orbit::Orbit(double area, double reflectivity, double mass, double Cd, bool compute_drag, bool compute_SRP, bool compute_third_body, int id){
    SetProperties(area, reflectivity, mass, Cd, compute_drag, compute_SRP, compute_third_body, id);
    suborbital = false;
}

//Copy constructor
Orbit::Orbit(const Orbit &O){
    //Sat properties
    satproperties = O.satproperties;
    //Perturbation flags
    Compute_Drag = O.Compute_Drag;
    Compute_SRP = O.Compute_SRP;
    Compute_Third_Body = O.Compute_Third_Body;
    //
    ID = O.ID;
    //Solution
    Soln = O.Soln;
    //suborbital flag
    suborbital = O.suborbital;
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
    Soln = Solution;
}

void Orbit::SetProperties(double area, double reflectance, double mass, double Cd, bool compute_drag, bool compute_SRP, bool compute_third_body, int id){
    //Set satellite properties and flags for perturbations
    satproperties.Area = area;
    satproperties.Mass = mass;
    satproperties.Reflectance = reflectance;
    satproperties.Cd = Cd;
    //Perturbation flags
    Compute_Drag = compute_drag;
    Compute_SRP = compute_SRP;
    Compute_Third_Body = compute_third_body;
    //ID
    ID = id;
}

void Orbit::SetSubOrbital(){
    suborbital = true;
}
