/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     May 2017
 * @ Modified by: Your name
 * @ Modified time: 2022-09-03 13:46:11
*  DESCRIPTION:      Set up an Adaptive-Picard-Chebyshev integration test case
*  REFERENCE:        Woollands, R., and Junkins, J., "Nonlinear Differential Equation Solvers
*                    via Adaptive Picard-Chebyshev Iteration: Applications in Astrodynamics", JGCD, 2016.
*/
#include <string>
#include <iostream>

#include <vector>
#include <iostream>
#include <utility>


#include <APC.h>
#include <adaptive_picard_chebyshev.h>
#include <c_functions.h>
#include <Orbit.h>
#include <EGM2008.h>
#include <Ephemeris.hpp>
#include "matrix_loader.h"
#include "flags.h"

using namespace std;

int main(){
  // MATRICES_LOADED=false;
  //satellite properties
  double mass = 1000;                               //sat mass (kg)
  double area = 10;                                 //sat wetted area (m^2)
  double reflectance = 1.5;                        //sat refelction absorption ratio
  double drag_C = 2.0;                              //sat coefficient of drag
  //Perturbation calc flags
  bool compute_drag = true;                         //atmostpheric drag toggle
  bool compute_SRP = true;                          //Solar radiation pressure toggle
  bool compute_third_body = true;                   //Third body gravity toggle
  //Ephemeris
  string spk = "de440.bsp";
  string lsk = "naif0012.tls";
  list<string> bodies = {"SUN","MOON"};
  string center = "Earth";
  string frame = "J2000";

  // Initialize Input Variables
  // LEO
  // vector<double> r0 = {7000, 0.0, 0.0};         // Initial Position (km)
  // vector<double> v0 = {0.0,  7.54604911, 0.};   // Initial Velocity (km/s)
  // double t0    = 0.0;                                // Initial Times (s)
  // double tf    = 10*5059.648765;                     // Final Time (s)

  // //SSO
  // vector<double> r0 = {6678,  0,  0};                  // Initial Position (km)
  // vector<double> v0 = {0,  0,  7.7451257}; // Initial Velocity (km/s)
  // double T = 5431.013011331035;                               //Orbital period (s)
  // double t0 = 0;                                              //initial time (s)
  // double tf = 10*T;                                              //final time (s)
       
  // MEO
  //std::vector<double> r0 = {9000.0, 0.0, 0.0};                                // Initial Position (km)
  //std::vector<double> v0 = {0.0, 6.7419845635570, 1.806509319188210};         // Initial Velocity (km/s)
  // double t0    = 0.0;                                               // Initial Times (s)
  // double tf    = 3.0*9.952014050491189e+03;                         // Final Time (s)
  // GEO
  // double r0[3] = {42000, 0.0, 0.0};                              // Initial Position (km)
  // double v0[3] = {0.0, 3.080663355435613, 0.0};                  // Initial Velocity (km/s)
  // double t0    = 0.0;                                            // Initial Times (s)
  // double tf    = 3.0*8.566135031791795e+04;                      // Final Time (s)
  // GTO
  // double r0[3] = {8064, 0.0, 0.0};                               // Initial Position (km)
  // double v0[3] = {0.0, 9.112725097814229, 0.0};                  // Initial Velocity (km/s)
  // double t0    = 0.0;                                            // Initial Times (s)
  // double tf    = 3.0*3.981179798339227e+04;                      // Final Time (s)
  // Molniya
  // double r0[3] = {7435.12, 0.0, 0.0};                            // Initial Position (km)
  // double v0[3] = {0.0, 4.299654205302486, 8.586211043023614};    // Initial Velocity (km/s)
  // double t0    = 0.0;                                            // Initial Times (s)
  // double tf    = 5.0*4.306316113361824e+04;                      // Final Time (s
  // Nan Orbit
  vector<double> r0 = {8000, 0, 0};                  // Initial Position (km)
  vector<double> v0 = {0.0, 7.2329973040227244, 0.0}; // Initial Velocity (km/s)
  double T = 5431.013011331035;                               //Orbital period (s)
  double t0 = 0;                                              //initial time (s)
  double tf = 7680;     
  double dt = 30;
  int steps = tf/dt+1;
  std::vector<double> time_vec;
  for(int jj=0;jj<steps;jj++)
  {
    double time = jj*dt;
    if(time>tf)
    {
      time=tf;
    }
    time_vec.push_back(time);
  }
  // EphemerisManager ephem(spk,lsk,t0,tf,bodies,center,frame);
  // MPGetTest(ephem, t0, tf);
  // std::cout << "Parallel Ephemeris Fetching Test Complete" << std::endl << "================================================" << std::endl;



  Orbit orb = SinglePropagate(r0, v0, time_vec,  area,  reflectance,  mass,  drag_C,  compute_drag,  compute_SRP,  compute_third_body);
  Orbit orb2 = SinglePropagate(r0, v0, t0, tf,  area,  reflectance,  mass,  drag_C,  compute_drag,  compute_SRP,  compute_third_body);
  // orb = SinglePropagate(r0, v0, t0 , tf,  area,  reflectance,  mass,  drag_C,  compute_drag,  compute_SRP,  compute_third_body);
  // orb = SinglePropagate(r0, v0, t0 , tf,  area,  reflectance,  mass,  drag_C,  compute_drag,  compute_SRP,  compute_third_body);
  std::cout << "Single Propagation Test Complete" << std::endl << "====================================" << std::endl;
  
  
  std::vector<SatState> sigma13 = GenSigma13(r0,v0,10,.1);
  int j = 0;
  std::vector<SatState> largelist;
  for (int i=0;i<100;i++){
    if (j>12){
      j=0;
    }
    largelist.push_back(sigma13[j]);
    j++;
  }
  std::vector<Orbit> orbits = ParallelPropagate(largelist, t0 , tf,  area,  reflectance,  mass,  drag_C,  compute_drag,  compute_SRP,  compute_third_body);
  std::cout << "Parallel Propagation Test Complete" << std::endl << "====================================" << std::endl;

  // std::pair<int,double> bench = Benchmark1000(8);
  // std::cout << "Benchmark with " << bench.first << " threads finished in " << bench.second << " seconds.\n";
}
