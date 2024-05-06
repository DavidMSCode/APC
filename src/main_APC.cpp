/*
*  AUTHORS:          David Stanley (DavidMS4@Illinois.edu)
*  DATE WRITTEN:     Nov 2023
 * @ Modified by: Your name
 * @ Modified time: 2023-12-24 17:44:44
*  DESCRIPTION:      Set up an Adaptive-Picard-Chebyshev integration test case
*  REFERENCE:        Woollands, R., and Junkins, J., "Nonlinear Differential Equation Solvers
*                    via Adaptive Picard-Chebyshev Iteration: Applications in Astrodynamics", JGCD, 2016.
*/
#include <string>
#include <iostream>
#include <math.h>
#include <vector>
#include <iostream>
#include <utility>
#include <unistd.h>
#include <fstream>
#include <APC.h>
#include <adaptive_picard_chebyshev.h>
#include <c_functions.h>
#include <Orbit.h>
#include <EGM2008.h>
#include <Ephemeris.hpp>
#include "matrix_loader.h"
#include "flags.h"
#include "TwoBody.h"
#include "const.h"
#include "EphemerisRotation.h"
using namespace std;

int main(int argc, char** argv){
  //PRINT THE CURRENT WORKING DIRECTORY
  char cwd[1024];
  if (getcwd(cwd, sizeof(cwd)) != NULL)
  {
    printf("Current working dir: %s\n", cwd);
  }
  else
  {
    perror("getcwd() error");
    return 1;
  }
  
  // MATRICES_LOADED=false;
  //satellite properties
  double mass = 212;                               //sat mass (kg)
  double area = 10;                                 //sat wetted area (m^2)
  double reflectance = 1.5;                        //sat refelction absorption ratio
  double drag_C = 2.0;                              //sat coefficient of drag
  //Perturbation calc flags
  bool compute_drag = false;                         //atmostpheric drag toggle
  bool compute_SRP = false;                          //Solar radiation pressure toggle
  bool compute_third_body = false;                   //Third body gravity toggle
  bool compute_hamiltonian = true;                 //whether or not the hamiltonian should be compuited for the output
  //Ephemeris
  string spk = "de440.bsp";
  string lsk = "naif0012.tls";
  list<string> bodies = {"SUN","EARTH"};
  string center = "MOON";
  string frame = "J2000";

  double alt = 300 ; //km
  double a = C_Rmoon+alt;
  double e = 0.0;
  double i = 0.0;
  double raan = 0.0;
  double aop = 0.0;
  double ta = 0.0;
  vector<vector<double>> states = elms2rv(a,e,i,raan,aop,ta,C_MU_MOON);
  double followtime = 3;
  vector<double> r0 = states[0];                // Initial Position (km)
  vector<double> v0 = states[1];
  double T = 2*C_PI*sqrt(pow(a,3)/C_MU_MOON);                             //Orbital period (s)
  double t0 = 0;                                              //initial time (s)
  double tf = t0+T;                                         //final time (s)

  //Orbit orb = SinglePropagate(r0, v0, time_vec,  area,  reflectance,  mass,  drag_C,  compute_drag,  compute_SRP,  compute_third_body);
  Orbit orbit("MOON","MOON_PA","J2000");
  orbit.SetProperties(area,reflectance,mass,drag_C);
  orbit.SetPosition0(r0);
  orbit.SetVelocity0(v0);
  orbit.SetIntegrationTime(t0,tf);
  // orbit.SetComputeThirdBody();
  // orbit.SetComputeSRP();
  orbit.SetComputeHamiltonian();
  orbit.SetMaxDegree(200);
  orbit.SetTolerance(1e-15);
  orbit.SetPolyDegreeParams(3,3);
  //run propagation
  orbit.SinglePropagate();
  // orb = SinglePropagate(r0, v0, t0 , tf,  area,  reflectance,  mass,  drag_C,  compute_drag,  compute_SRP,  compute_third_body);
  // orb = SinglePropagate(r0, v0, t0 , tf,  area,  reflectance,  mass,  drag_C,  compute_drag,  compute_SRP,  compute_third_body);
  std::cout << "Single Propagation Test Complete" << std::endl << "====================================" << std::endl;
  
  return 0;
}
