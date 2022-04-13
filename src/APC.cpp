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

#include "SpiceUsr.h"
#include <adaptive_picard_chebyshev.h>
#include <c_functions.h>
#include <EGM2008.h>
#include <time.h> 
#include <errno.h>
#include <vector>
#include <Orbit.h>
#include <APC.h>
#include <iostream>
#include <sstream>
#include <iterator>
#include <Ephemeris.hpp>

class EphemerisManager cacheEphemeris(double t0, double tf){
  //Ephemeris
  string spk = "de440.bsp";               //Ephemeris file for sun earth and moon
  string lsk = "naif0012.tls";            //leap second kernel
  list<string> bodies = {"SUN","MOON"};   //Required bodies to store
  string center = "Earth";                //Observing body
  string frame = "J2000";                 //Frame
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
  std::stringstream output;
  std::string vecstring;
  std::copy(vector.begin(), vector.end(), std::ostream_iterator<double>(output, ", "));
  vecstring = output.str();
  vecstring.pop_back();
  vecstring.pop_back();
  return "{"+vecstring+"}";
}

void printStateList(std::vector<SatState> statelist){
  int N = statelist.size();
  for(int i=0;i<N;i++){
      printState(statelist[i],i);
      if (i!=N-1){
        std::cout << std::endl;
      }
    }
}

std::vector<std::vector<double> > PropagateICs(std::vector<double> r, std::vector<double> v, double t0, double tf, Orbit orb, EphemerisManager ephem){
  //Convert vectors to array since pybind wants vectors but the functions are coded for arrays
  double* r0 = &r[0];
  double* v0 = &v[0];
  double dt    = 30.0;                             // Soution Output Time Interval (s)
  double deg   = 70.0;                             // Gravity Degree (max 100)
  double tol   = 1.0e-15;                          // Tolerance
  // Initialize Output Variables
  int soln_size = int(1.1*(tf/dt));
  if (soln_size == 1){
    soln_size = 2;
  }
  std::vector<double> Soln(soln_size*6,0.0);
  //Soln = static_cast<double*>(calloc(soln_size*6,sizeof(double)));       // Position (km) & Velocity (km/s)

  double Feval[2] = {0.0};
  std::vector<std::vector<double> > states;
  states = adaptive_picard_chebyshev(r0,v0,t0,tf,dt,deg,tol,soln_size,Feval,Soln,orb,ephem);
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

  double t_curr = t0;
  for (int i=1; i<=soln_size; i++){
    Ts.push_back(t_curr);
    for (int j=1; j<=6; j++){
      States[j-1].push_back(Soln[ID2(i,j,soln_size)]);
      state[j-1] = Soln[ID2(i,j,soln_size)];
    }
    jacobiIntegral(t_curr,state,&H,deg);
    if (i == 1){
      H0 = H;
    }
    if (fabs((H-H0)/H0) > Hmax){
      Hmax = fabs((H-H0)/H0);
    }
    Hs.push_back(fabs((H-H0)/H0));
    t_curr = t_curr + dt;
    if (t_curr > tf){
      break;
    }
  }
  std::ostringstream ss;
  ss << Hmax;
  std::string HmaxStr(ss.str());

  std::cout << to_string(orb.ID)+":\tFunc Evals: " + to_string(total) + "  \t Hmax: " + HmaxStr + "\n";
  // printf("Func Evals: %i\t",total);
  // printf("Hmax %1.16E\n",Hmax);
  //Assemble solution vector
  std::vector<std::vector<double> > Solution;
  Solution.push_back(Ts);
  for(int i=0; i<=5; i++){
    Solution.push_back(States[i]);
  }
  Solution.push_back(Hs);
  //free(Soln);
  return Solution;
}


class Orbit PropagateOrbit(std::vector<double> r, std::vector<double> v, double t0, double tf, Orbit orbit, EphemerisManager ephem){
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

std::vector<Orbit> ParallelPropagate(std::vector<SatState> StateList, double t0, double tf, double area, double reflectance, double mass, double drag_C, bool compute_drag, bool compute_SRP, bool compute_third_body){
  int n = StateList.size();
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
      Orbit orbit(area_priv,reflectance_priv,mass_priv,drag_C_priv,compute_drag_priv,compute_SRP_priv,compute_third_body_priv,i);
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


class Orbit SinglePropagate(std::vector<double> r, std::vector<double> v, double t0, double tf, double area, double reflectance, double mass, double drag_C, bool compute_drag, bool compute_SRP, bool compute_third_body){
  EphemerisManager ephem = cacheEphemeris(t0,tf+3600);
  Orbit orbit(area,reflectance,mass,drag_C,compute_drag,compute_SRP,compute_third_body,1);
  Orbit orbit2 = PropagateOrbit(r, v,  t0,  tf,  orbit, ephem);
  return orbit2;
}

void MPGetTest(EphemerisManager ephem, double t0, double tf){
  int N = 1000;
  double dt = (tf-t0)/N;
  #pragma omp parallel for
    for(int i=0;i<1000;i++){
      double epoch = dt*i+t0;
      std::vector<double> sunstate = ephem.getState("SUN",epoch);
      std::vector<double> moonstate = ephem.getState("MOON",epoch);
      std::string sunvec = vec2prettystring(sunstate);
      std::string moonvec = vec2prettystring(moonstate);
      std::cout << to_string(epoch)+"s: Sun @ "+sunvec+" Moon @ " + moonvec + "\n";
    }
}