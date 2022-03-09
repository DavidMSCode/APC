/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     May 2017
*  LAST MODIFIED:    Aug 2021
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Set up an Adaptive-Picard-Chebyshev integration test case
*  REFERENCE:        Woollands, R., and Junkins, J., "Nonlinear Differential Equation Solvers
*                    via Adaptive Picard-Chebyshev Iteration: Applications in Astrodynamics", JGCD, 2016.
*/

#include <APC.h>
#include <adaptive_picard_chebyshev.h>
#include <c_functions.h>
#include <Orbit.h>
#include <EGM2008.h>
#include <time.h>  
#include <errno.h>
#include <vector>

int main(){
  //satelltie properties
  double mass = 1000;                               //sat mass (kg)
  double area = 10;                                 //sat wetted area (m^2)
  double reflectance = 1.5;                        //sat refelction absorption ratio
  double drag_C = 2.0;                              //sat coefficient of drag
  //Perturbation calc flags
  bool compute_drag = true;                         //atmostpheric drag toggle
  bool compute_SRP = true;                          //Solar radiation pressure toggle
  bool compute_third_body = true;                   //Third body gravity toggle

  // Initialize Input Variables
  // LEO
  std::vector<double> r0 = {0.0, -8000.0, 0.0};      // Initial Position (km)
  std::vector<double> v0 = {7.0586826,  0.0, 0.0};   // Initial Velocity (km/s)
  double t0    = 0.0;                                // Initial Times (s)
  double tf    = 10*5059.648765;                     // Final Time (s)
  // MEO
  // double r0[3] = {9000.0, 0.0, 0.0};                                // Initial Position (km)
  // double v0[3] = {0.0, 6.7419845635570, 1.806509319188210};         // Initial Velocity (km/s)
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
  // double tf    = 5.0*4.306316113361824e+04;                      // Final Time (s)
  Orbit orb = PropagateOrbit(r0, v0, t0 , tf,  area,  reflectance,  mass,  drag_C,  compute_drag,  compute_SRP,  compute_third_body);
  //free(Soln);

}
