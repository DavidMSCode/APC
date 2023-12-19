/**
 * @file APC.cpp
 * @author David Stanley (davidms4@illinois.edu)
 * @brief 
 * @version 0.1
 * @date 2022-04-13
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include <time.h> 
#include <errno.h>
#include <math.h>

#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iterator>
#include <utility>
#include <chrono>
#include <list>
#include <string>
#include <omp.h>
#include "SpiceUsr.h"

#include "APC.h"
#include "adaptive_picard_chebyshev.h"
#include "c_functions.h"
#include "EGM2008.h"
#include "GRGM1200b.h"
#include "Orbit.h"
#include "Ephemeris.hpp"
#include "flags.h"

using namespace std;

EphemerisManager cacheEphemeris(double t0, double tf, const string center, const list<string> bodies,string frame){
  //Ephemeris
  //FIXME: These strings should be user modifiable
  string spk = "de440.bsp";               //Ephemeris file for sun earth and moon
  string lsk = "naif0012.tls";            //leap second kernel
  EphemerisManager ephem(spk,lsk,t0,tf,bodies,center,frame);
  return ephem;
}

void printState(SatState state, int i){
  std::stringstream output1;
  std::string pos;
  std::copy(state.r.begin(), state.r.end(), std::ostream_iterator<double>(output1, ", "));
  pos = output1.str();
  pos.pop_back();
  pos.pop_back();
  std::cout << "Orbit " << i << " position: {" << pos << "}"<< std::endl;
  std::stringstream output2;
  std::string vel;
  std::copy(state.v.begin(), state.v.end(), std::ostream_iterator<double>(output2, ", "));
  vel = output2.str();
  vel.pop_back();
  vel.pop_back();
  std::cout << "Orbit " << i << " velocity: {" << vel <<"}"<< std::endl;
}

string vec2prettystring(std::vector<double> vector){
  //Takes a vector of doubles and returns a string of the form {x,y,z}
  std::stringstream output;
  std::string vecstring;
  std::copy(vector.begin(), vector.end(), std::ostream_iterator<double>(output, ", "));
  vecstring = output.str();
  vecstring.pop_back();
  vecstring.pop_back();
  return "{"+vecstring+"}";
}

void printStateList(std::vector<SatState> statelist){
  size_t N = statelist.size();
  for(auto i=0;i<N;i++){
      printState(statelist[i],i);
      if (i!=N-1){
        std::cout << std::endl;
      }
    }
}

void LoadKernels(){
  if(not g_KERNELS_LOADED){
    //Load spice kernels
    string spkFile = "de440.bsp";
    string lskFile = "naif0012.tls";
    furnsh_c(spkFile.c_str());
    furnsh_c(lskFile.c_str());
    furnsh_c("moon_de440_220930.tf");
    furnsh_c("moon_pa_de440_200625.bpc");
    furnsh_c("pck00010.tpc");
    g_KERNELS_LOADED = true;
  }
}

void UnloadKernels(){
  if(g_KERNELS_LOADED){
    string spkFile = "de440.bsp";
    string lskFile = "naif0012.tls";
    unload_c(spkFile.c_str());
    unload_c(lskFile.c_str());
    unload_c("moon_de440_220930.tf");
    unload_c("moon_pa_de440_200625.bpc");
    unload_c("pck00010.tpc");
    g_KERNELS_LOADED = false;
  }
}

std::vector<std::vector<double> > PropagateICs(std::vector<double> r, std::vector<double> v, double t0, double tf, Orbit &orbit, EphemerisManager &ephem){
  LoadKernels();
  //Convert vectors to array since pybind wants vectors but the functions are coded for arrays
  double* r0 = &r[0];
  double* v0 = &v[0];
  double dt    = 30.0;                              // Soution Output Time Interval (s)
  //FIXME: deg should be user modifiable
  double deg   = 250;                               // Gravity Degree (max 100)
  orbit.deg = deg;
  // Initialize Output Variables
  int soln_size = int(ceil((tf/dt)))+1;
  if (soln_size == 1){
    soln_size = 2;
  } 
  std::vector<double> Soln(soln_size*6,0.0);
  double Feval[2] = {0.0};
  std::vector<std::vector<double> > states;
  adaptive_picard_chebyshev(Feval,Soln,orbit,ephem);
  int total;
  total = int(ceil(Feval[0] + Feval[1]*pow(6.0,2)/pow(deg,2)));

  // Assemble solution vector from solution array
  std::vector<std::vector<double> > States(6);
  double state[6] = {0.0};
  std::vector<double> Hs;
  std::vector<double> Ts;
  double H    = 0.0;
  double H0   = 0.0;
  double Hmax = 0.0;
  std::vector<double> t_vec = orbit.T;
  double t_curr;
  soln_size = t_vec.size();
  std:string HmaxStr;
  for (int i=1; i<=soln_size; i++){
    t_curr = orbit.et(t_vec[i-1]);
    for (int j=1; j<=6; j++){
      States[j-1].push_back(Soln[ID2(i,j,soln_size)]);
      state[j-1] = Soln[ID2(i,j,soln_size)];
    }
    //Compute Hamiltonian
    if(orbit.Compute_Hamiltonian)
    {
      jacobiIntegral_GRGM1200b(t_curr,state,&H,deg,orbit);
      if (i == 1){
        H0 = H;
      }
      if (fabs((H-H0)/H0) > Hmax){
        Hmax = fabs((H-H0)/H0);
      }
    }
    //Store hamiltonian
    Hs.push_back(fabs((H-H0)/H0));
  std::ostringstream ss;
  ss << Hmax;
  HmaxStr = ss.str();
  }
  std::cout << to_string(orbit.ID)+":\tFunc Evals: " + to_string(total) + "  \t Hmax: " + HmaxStr + "\n";
  //Assemble solution vector
  std::vector<std::vector<double> > Solution;
  Solution.push_back(t_vec);
  for(int i=0; i<=5; i++){
    Solution.push_back(States[i]);
  }
  Solution.push_back(Hs);
  UnloadKernels();        //Unload the kernels from SPICE
  return Solution;
}

//FIXME: This should be an orbit class method
class Orbit PropagateOrbit(std::vector<double> r, std::vector<double> v, double t0, double tf, Orbit &orbit, EphemerisManager &ephem){
  //Generates 13 orbits with slight perturbations in state +/- on each coordinate
  std::vector<std::vector<double> > solution;
  solution = PropagateICs(r, v, t0 , tf, orbit, ephem);
  orbit.SetSolution(solution);
  return orbit;
}

std::vector<SatState> GenSigma13(std::vector<double> r, std::vector<double> v, double pos_error, double vel_error){
  std::vector<SatState> Sigma13(13);
  std::vector<double> plusR(3);
  std::vector<double> minR(3);
  std::vector<double> plusV(3);
  std::vector<double> minV(3);
  //IC with no error
  Sigma13[0].r = r;
  Sigma13[0].v = v;

  
  for(int i=0;i<3;i++){
    //Generate 2 ICs per DOF 
    plusR = r;
    minR = r;
    plusV = v;
    minV = v;
    plusR[i] += pos_error;
    minR[i] -= pos_error;
    plusV[i] += vel_error;
    minV[i] -= vel_error;
    
    //plus position
    Sigma13[1+4*i].r = plusR;
    Sigma13[1+4*i].v = v;
    //min position
    Sigma13[2+4*i].r = minR;
    Sigma13[2+4*i].v = v;
    //plus velocity
    Sigma13[3+4*i].r = r;
    Sigma13[3+4*i].v = plusV;
    //min velocity
    Sigma13[4+4*i].r = r;
    Sigma13[4+4*i].v = minV;
  }
  return Sigma13;
}

std::vector<SatState> GenSigma3(std::vector<double> r, std::vector<double> v, double pos_error, double vel_error){
  //Generates 3 orbits with slight perturbations in first state coordinate +/-
  std::vector<SatState> Sigma3(3);
  std::vector<double> plusR(3);
  std::vector<double> minR(3);
  std::vector<double> plusV(3);
  std::vector<double> minV(3);
  //IC with no error
  Sigma3[0].r = r;
  Sigma3[0].v = v;

  
  int i=0;
  //Generate 2 new ICs 
  plusR = r;
  minR = r;
  plusR[i] += pos_error;
  minR[i] -= pos_error;
  //plus position
  Sigma3[1+4*i].r = plusR;
  Sigma3[1+4*i].v = v;
  //min position
  Sigma3[2+4*i].r = minR;
  Sigma3[2+4*i].v = v;
  return Sigma3;
}

std::vector<Orbit> ParallelPropagate(std::vector<SatState> StateList, double t0, double tf, double area, double reflectance, double mass, double drag_C, bool compute_drag, bool compute_SRP, bool compute_third_body,bool compute_hamiltonian){
  size_t n = StateList.size();
  int threads;
  //std::cout<<"There are "+to_string(threads)+" available threads.\n";
  //cache ephemeris data for time range
  EphemerisManager ephem = cacheEphemeris(t0-1000,tf+1000);
  std::vector<Orbit> orbits;
  //Spawn threads and start parallel for loop
  #pragma omp parallel for shared(orbits, n)
    for (int i=0;i<n;i++){
      //Thread debug print
      //std::cout<< "Orbit "+to_string(i)+" assigned to thread "+to_string(omp_get_thread_num())+".\n";
      //Private varaibles
      double t0_priv = t0;
      double tf_priv = tf;
      double area_priv = area;
      double reflectance_priv = reflectance;
      double mass_priv = mass;
      double drag_C_priv = drag_C;
      bool compute_drag_priv = compute_drag;
      bool compute_SRP_priv = compute_SRP;
      double compute_third_body_priv = compute_third_body;
      //Private initial state
      std::vector<double> r0 = StateList[i].r;
      std::vector<double> v0 = StateList[i].v;
      //Private orbit object
      Orbit orbit(area_priv,reflectance_priv,mass_priv,drag_C_priv,compute_drag_priv,compute_SRP_priv,compute_third_body_priv,compute_hamiltonian,i);
      //thread debug

      //std::cout << "Running thread "+to_string(omp_get_thread_num())+".\n";
      orbit = PropagateOrbit(r0, v0,  t0_priv,  tf_priv,  orbit, ephem);
      #pragma omp critical(writeout)
      {
        //Store solution
        orbits.push_back(orbit);
        //std::cout << "Thread "+to_string(omp_get_thread_num())+" finished sucessfully.\n";
      }
    }
return orbits;
}


class Orbit SinglePropagate(std::vector<double> r, std::vector<double> v, double t0, double tf, double area, double reflectance, double mass, double drag_C, bool compute_drag, bool compute_SRP, bool compute_third_body, bool compute_hamiltonian, string center){
  EphemerisManager ephem = cacheEphemeris(t0,tf+3600);
  //Initialize orbit object
  Orbit orbit(area,reflectance,mass,drag_C,compute_drag,compute_SRP,compute_third_body,compute_hamiltonian,1,center);
  double dt  = 30;
  double len = int(ceil(tf/dt));
  // std::vector<double> t_vec(len+1,0.0);
  // t_vec[0] = t0;
  // for (int ii=1; ii<=len; ii++){
  //   double time = t_vec[ii-1] + dt;
  //   if(time>tf)
  //   {
  //     time = tf;
  //   }
  //   t_vec[ii] = time;

  // }
  // orbit.SetTimeVec(t_vec);
  orbit = PropagateOrbit(r, v,  t0,  tf,  orbit, ephem);
  return orbit;
}
//overload using a user defined time vector
class Orbit SinglePropagate(std::vector<double> r, std::vector<double> v, std::vector<double> time_vec, double area, double reflectance, double mass, double drag_C, bool compute_drag, bool compute_SRP, bool compute_third_body, bool compute_hamiltonian, string center){
  double t0 = time_vec[0];
  double tf = time_vec.back()+300;
  EphemerisManager ephem = cacheEphemeris(t0,tf+3600);
  Orbit orbit(area,reflectance,mass,drag_C,compute_drag,compute_SRP,compute_third_body, compute_hamiltonian,1);
  orbit.SetTimeVec(time_vec);
  Orbit orbit2 = PropagateOrbit(r, v,  t0,  tf,  orbit, ephem);
  return orbit2;
}



std::pair<int,double>  Benchmark1000(int max_threads){
  //Benchmarking function that solves 1000 orbits and reports the time for completion.
  omp_set_num_threads(max_threads);
  int procs = omp_get_num_procs();
  auto start = std::chrono::steady_clock::now();
  double t0 = 0;
  double tf = 5059.648765;
  double time = 0;
  //number of seperate orbits
  int n = 1000;
  EphemerisManager ephem = cacheEphemeris(t0-1000,tf+1000);
  std::vector<Orbit> orbits;
  int j = 0;
  //Spawn threads and start parallel for loop
  #pragma omp parallel for shared(orbits, n, j)
    for (int i=0;i<n;i++){
      //Thread debug print
      //std::cout<< "Orbit "+to_string(i)+" assigned to thread "+to_string(omp_get_thread_num())+".\n";
      //Private varaibles
      double t0_priv = t0;
      double tf_priv = tf;
      double area_priv = 10;
      double reflectance_priv = 1.5;
      double mass_priv = 1000;
      double drag_C_priv = 2.0;
      //Private initial state
      std::vector<double> r0 = {8000, 0, 0};
      std::vector<double> v0 = {0, 8, 0};
      //Private orbit object
      Orbit orbit(area_priv,reflectance_priv,mass_priv,drag_C_priv,false,false,false,false,i);
      //thread debug

      //std::cout << "Running thread "+to_string(omp_get_thread_num())+".\n";
      orbit = PropagateOrbit(r0, v0,  t0_priv,  tf_priv,  orbit, ephem);
      #pragma omp critical(writeout)
      {
        //Store solution even though it isn't returned as the writeout lock causes some overhead
        orbits.push_back(orbit);
        //std::cout << "Thread "+to_string(omp_get_thread_num())+" finished sucessfully.\n";
      }
    }
    //get time at end
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff = end-start;
    //return number of threads and the time to finish
    std::pair<int,double> out = {std::min(max_threads,procs), diff.count()};
    return out;
}

void printMatrixState(){
  if(g_MATRICES_LOADED){
    std::cout << "The matrices are loaded.\n";
  }
  else{
    std::cout << "The matrices are not loaded.\n";
  }
}

bool MatricesLoaded(){
  return g_MATRICES_LOADED;
}