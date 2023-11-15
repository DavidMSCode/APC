/**
 * @file EphemerisManagerTests.cpp
 * @author David Stanley (davidms4@illinois.edu)
 * @ Modified by: Your name
 * @ Modified time: 2023-10-17 14:12:20
 * @date 2023-03-10
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include <iostream>
#include <stdexcept>
#include "Ephemeris.hpp"
#include "EphemerisManagerTests.h"
#include <string>
#include "APC.h"
#include "omp.h"
#include "SpiceUsr.h"
using namespace std;

int EphemerisManagerTests(int argc, char** argv){
    furnsh_c("naif0012.tls");
    SpiceDouble t0;
    str2et_c("Mon Sep 30 09:59:10 PDT 2023", &t0);
    double tf = t0 + 5*24*60*60;
    unload_c("naif0012.tls");
    EphemerisManager ephem("de440.bsp","naif0012.tls",t0,tf,{"SUN","MOON","EARTH"},"SOLAR SYSTEM BARYCENTER","J2000");
    int errs = 0;
    errs += MPGetTest(ephem, t0+800,tf-800);
    return errs;
}

int MPGetTest(EphemerisManager &ephem, double t0, double tf){
  int N = 1000;
  double dt = (tf-t0)/N;
  int errs = 0;
  list<string> bodies = ephem.getBodies();
//   list<string> bodies = {"SATURN"};
  // Track errors per thread then sum together at end of parallel section
  #pragma omp parallel for reduction (+ : errs)
    for(int i=0;i<1000;i++){
        double epoch = dt*i+t0;
        string fullstring = to_string(epoch)+"s: ";
        for(list<string>::iterator it = bodies.begin(); it != bodies.end(); it++){
            string target = *it;
            try {
                std::vector<double> targetstate = ephem.getState(target,epoch);
                std::string targetvec = vec2prettystring(targetstate);
                // fullstring = fullstring+ *it +" @ "+targetvec+" ";
            }
            catch(std::invalid_argument const& ex){
                // fullstring = fullstring + *it +" state could not be retrieved.";
                errs += 1;
                if(errs%100==0){
                    cout << "MPGetTest: Thread " + to_string(omp_get_thread_num()) + " has had " + to_string(errs) + " retrieval errors so far." << endl;
                }
            }
        }
        // std::cout << fullstring << endl;
    }
    if(errs>0){
        cout<<"Ephemeris Manager failed to retrieve target states " + to_string(errs) + " time(s)."<<endl; 
    }
    return errs;
}
 
