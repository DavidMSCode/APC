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
#include <sstream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <iostream>
#include <utility>
#include <unistd.h>

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
#include "interpolate.h"
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
  string center = "EARTH";
  string frame = "J2000";

  double alt = 360 ; //km
  double a = 7000;
  double e = 0;
  double i = 0.0;
  double raan = 0.0;
  double aop = 0.0;
  double ta = 0.0;

  vector<vector<double>> states = elms2rv(a,e,i,raan,aop,ta,C_MU_EARTH);
  double followtime = 0;
  vector<double> r0 = states[0];                // Initial Position (km)
  vector<double> v0 = states[1];
  double T = 2*C_PI*sqrt(pow(a,3)/C_MU_EARTH);                             //Orbital period (s)
  double t0 = 0;                                              //initial time (s)
  double tf = t0+T;                                         //final time (s)

  //Orbit orb = SinglePropagate(r0, v0, time_vec,  area,  reflectance,  mass,  drag_C,  compute_drag,  compute_SRP,  compute_third_body);
  Orbit orbit("EARTH","EARTH_IAU","J2000");
  orbit.SetProperties(area,reflectance,mass,drag_C);
  orbit.SetPosition0(r0);
  orbit.SetVelocity0(v0);
  orbit.SetIntegrationTime(t0,tf);
  // orbit.SetComputeThirdBody();
  // orbit.SetComputeSRP();
  orbit.SetComputeHamiltonian();
  orbit.SetMaxDegree(70);
  orbit.SetTolerance(1e-14);
  //run propagation
  InterpolatedOrbit InterpolatedOrbit(orbit, followtime);
  InterpolatedOrbit.InterpolatePropagate();
  vector<double> node_soln = nodes(InterpolatedOrbit);
  int M  = InterpolatedOrbit.M;
  int segs = InterpolatedOrbit.total_segs;
  int soln_size = (M+1)*segs;
  //seperate out the x y z components of the nodes from the soln vector
  vector<double> x_nodes(soln_size,0.0);
  vector<double> y_nodes(soln_size,0.0);
  vector<double> z_nodes(soln_size,0.0);
  vector<double> t_nodes(soln_size,0.0);
  for(int i=1;i<=soln_size;i++)
  {
    x_nodes[i-1] = node_soln[ID2(i,1,soln_size)];
    y_nodes[i-1] = node_soln[ID2(i,2,soln_size)];
    z_nodes[i-1] = node_soln[ID2(i,3,soln_size)];
    t_nodes[i-1] = node_soln[ID2(i,7,soln_size)];
  }
  
  stringstream ss;
  ss << std::fixed << std::setprecision(3);
  ss << "Nodes_a="<<a<<"_e="<<e<<".csv";
  string filename = ss.str();

  ofstream myfile;
  myfile.open(filename);
  myfile << fixed << setprecision(16);
  //number of nodes poly degree and the number of segments

  myfile << "Orbit had the following properties a="<<a<<" e="<<e<<" i="<<i<<" raan="<<raan<<" aop="<<aop<<" ta="<<ta;
  myfile << " Number of nodes per segment="<<M+1<<" Number of segments="<<segs <<" Total nodes="<<soln_size<<" Polynomial degree="<<M<<"\n";
  //write the orbital elements to the file header
  myfile << "x,y,z,t\n";
  for(int i=0;i<soln_size;i++)
  {
    myfile << x_nodes[i] << "," << y_nodes[i] << "," << z_nodes[i] << "," << t_nodes[i] << "\n";
  }
  myfile.close();


  std::cout << "Single Propagation Test Complete" << std::endl << "====================================" << std::endl;
  
  std::cout << "Bootstrap orbit test starting" << std::endl << "====================================" << std::endl;

  std::cout << "Bootstrap orbit test complete" << std::endl << "====================================" << std::endl;
  return 0;
  
}
