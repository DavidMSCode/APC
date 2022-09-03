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

// //Copy constructor
// Orbit::Orbit(const Orbit &O){
//     //Sat properties
//     satproperties = O.satproperties;
//     //Perturbation flags
//     Compute_Drag = O.Compute_Drag;
//     Compute_SRP = O.Compute_SRP;
//     Compute_Third_Body = O.Compute_Third_Body;
//     //
//     ID = O.ID;
//     //Solution
//     Soln = O.Soln;
//     //suborbital flag
//     suborbital = O.suborbital;
//     //User defineds time vector flag
//     USER_TIME = false;
// }
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
