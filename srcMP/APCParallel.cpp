/**
 * @file APCParallel.cpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2022-03-23
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
#include <omp.h>
#include <APCParallel.h>

void ParallelPropagateTest(std::vector<SatState> StateList, double t0, double tf, double area, double reflectance, double mass, double drag_C, bool compute_drag, bool compute_SRP, bool compute_third_body){
//Spawn threads
  int n = StateList.size();
  furnsh_c("de440.bsp");
  #pragma omp parallel for 
  {
    for (int i=0;i<n;i++){
      #pragma omp critical
      {
        std::cout << "Using thread " << omp_get_thread_num()<<" for orbit #"<<i<<std::endl;
      }
      std::vector<double> r0 = StateList[i].r;
      std::vector<double> v0 = StateList[i].v;
      Orbit orbit = PropagateOrbit(r0, v0,  t0,  tf,  area,  reflectance,  mass,  drag_C,  compute_drag,  compute_SRP,  compute_third_body);
    }
  }
  kclear_c();

}