/*
 *  AUTHORS:          David Stanley (DavidMS4@Illinois.edu)
 *  DATE WRITTEN:     Nov 2023
 * @ Modified by: Your name
 * @ Modified time: 2023-12-24 17:44:44
 *  DESCRIPTION:      Calculates the average computation time for differnet segment lengths and node counts whil calculating the gravity per segment in parallel.
 *  REFERENCE:        Woollands, R., and Junkins, J., "Nonlinear Differential Equation Solvers
 *                    via Adaptive Picard-Chebyshev Iteration: Applications in Astrodynamics", JGCD, 2016.
 */
#include <string>
#include <iostream>
#include <sstream>
#include <time.h>
#include <numeric>
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
#include <omp.h>
using namespace std;

void print_progress(double progress)
{
    int barWidth = 70;

    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}

double cpu_time(void)
{
  /*Calculates cpu time in cpu seconds*/
  double value;
  value = (double)clock() / (double)CLOCKS_PER_SEC;
  return value;
}

int main(int argc, char **argv)
{

  // PRINT THE CURRENT WORKING DIRECTORY
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

  int n_trials = 20;
  std::vector<int> min_seg_params = {1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 51, 101}; // minimum number of segments per orbit
  std::vector<int> j_max_params = {5, 4, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1};              // max node doublings while finding desired polynomial degree. jmax={1,2,3,4,5} corrersponds to Nmax = {20,40,80,160,320}
  int n_params = min_seg_params.size();
  std::vector<int> Ns(n_params, 0);   // preallocate vector for storing the polynomial degree for each run
  std::vector<int> segs(n_params, 0); // preallocate vector for storing the number of segments for each run

  // preallocate vectors for storing 10 runs of real time
  std::vector<double> avg_time(n_params, 0.0);
  std::vector<std::vector<double>> real_times(n_params, std::vector<double>(n_trials, 0.0));
  // preallocate vectors for storing 10 runs of total cpu time
  std::vector<double> avg_tick_time(n_params, 0.0);
  std::vector<std::vector<double>> tick_times(n_params, std::vector<double>(n_trials, 0.0));

  // MATRICES_LOADED=false;
  // satellite properties
  double mass = 212;        // sat mass (kg)
  double area = 10;         // sat wetted area (m^2)
  double reflectance = 1.5; // sat refelction absorption ratio
  double drag_C = 2.0;      // sat coefficient of drag
  // Perturbation calc flags
  bool compute_drag = false;       // atmostpheric drag toggle
  bool compute_SRP = false;        // Solar radiation pressure toggle
  bool compute_third_body = false; // Third body gravity toggle
  bool compute_hamiltonian = true; // whether or not the hamiltonian should be compuited for the output
  // Ephemeris
  string spk = "de440.bsp";
  string lsk = "naif0012.tls";
  list<string> bodies = {"SUN", "EARTH"};
  string center = "MOON";
  string frame = "J2000";

  double alt = 500; // km
  double a = C_Rmoon + alt;
  double e = 0.0;
  double inc = 0.0;
  double raan = 0.0;
  double aop = 0.0;
  double ta = 0.0;
  vector<vector<double>> states = elms2rv(a, e, inc, raan, aop, ta, C_MU_MOON);
  vector<double> r0 = states[0]; // Initial Position (km)
  vector<double> v0 = states[1];
  double T = 2 * C_PI * sqrt(pow(a, 3) / C_MU_MOON); // Orbital period (s)
  double t0 = 0;                                     // initial time (s)
  double tf = t0 + 10 * T;                           // final time (s)

  // set random seed
  srand(time(NULL));
  int i;
  std::vector<int> js(n_params, 0);
  omp_set_num_threads(12);
  int thread_count = 0;
#pragma omp parallel
  {
    thread_count = omp_get_num_threads();
  }
  while (std::accumulate(js.begin(), js.end(), 0) < n_trials * n_params)
  {
    // choose a random int for i from 0 to n_params-1
    while (true)
    {
      i = rand() % n_params;
      // cout << "rolled i = " << i << endl;
      if (js[i] < n_trials)
      {
        js[i]++;
        break;
      }
    }
    int j_max = j_max_params[i];
    int min_seg = min_seg_params[i];
    // cout << "Running test with min_seg = " << min_seg << " and j_max = " << j_max << endl;
    // choose a random

    Orbit orbit("MOON", "MOON_PA", "J2000");
    orbit.SetProperties(area, reflectance, mass, drag_C);
    orbit.SetPosition0(r0);
    orbit.SetVelocity0(v0);
    orbit.SetIntegrationTime(t0, tf);
    orbit.SetComputeHamiltonian(false);
    orbit.SetMaxDegree(200);
    orbit.SetTolerance(1e-15);
    orbit.SetPolyDegreeParams(min_seg, j_max);

    // start real world timer
    double start = omp_get_wtime();
    // start cpu timer
    double start_cpu = cpu_time();
    orbit.SinglePropagate();
    // stop omp parallel timer
    double end = omp_get_wtime();
    // stop cpu timer
    double end_cpu = cpu_time();

    // store real time in secs
    real_times[i][js[i] - 1] = end - start;
    // store cpu time in secs
    tick_times[i][js[i] - 1] = end_cpu - start_cpu;
    // print progress
    print_progress((double)std::accumulate(js.begin(), js.end(), 0) / (n_trials * n_params));
  }
  cout<<endl;
  std::vector<double> max_Hs(n_params, 0.0);
  std::vector<int> Fevals(n_params, 0);
  // run one more set but don't record the times just get the max hamiltonian error and the number of function evals
  for (int i = 0; i < n_params; i++)
  {
    int j_max = j_max_params[i];
    int min_seg = min_seg_params[i];
    cout << "Running test with min_seg = " << min_seg << " and j_max = " << j_max << endl;
    Orbit orbit("MOON", "MOON_PA", "J2000");
    orbit.SetProperties(area, reflectance, mass, drag_C);
    orbit.SetPosition0(r0);
    orbit.SetVelocity0(v0);
    orbit.SetIntegrationTime(t0, tf);
    orbit.SetComputeHamiltonian(true);
    orbit.SetMaxDegree(200);
    orbit.SetTolerance(1e-15);
    orbit.SetPolyDegreeParams(min_seg, j_max);
    orbit.SinglePropagate();
    Ns[i] = orbit.N;
    segs[i] = orbit.seg;
    Fevals[i] = orbit.TotalFuncEvals;
    max_Hs[i] = orbit.Hmax;
  }

  // write the real times to a file with columns Ns, segs, Fevals, max_Hs, real_times trial 1, 2, 3, etc.
  ofstream file;
  // get date
  time_t now = time(0);
  tm *ltm = localtime(&now);
  string date = to_string(1 + ltm->tm_mon) + "-" + to_string(ltm->tm_mday) + "-" + to_string(1900 + ltm->tm_year);
  // open the file
  string filename_string = "Parallel segment times " + to_string(thread_count)+" threads_"+ date + ".csv";
  file.open(filename_string.c_str());
  // set precision
  file.precision(15);
  // write the data
  for (int i = 0; i < n_params; i++)
  {
    file << Ns[i] << ", " << segs[i] << ", " << Fevals[i] << ", " << max_Hs[i] << ", ";
    for (int j = 0; j < n_trials; j++)
    {
      file << real_times[i][j];
      file << ", ";
    }
    for (int j = 0; j < n_trials; j++)
    {
      file << tick_times[i][j];
      if (j < n_trials - 1)
      {
        file << ", ";
      }
    }
    file << endl;
  }
  file.close();

  string filename_string2 = "Parallel segment times " + to_string(thread_count)+" threads header_"+ date + ".txt";
  file.open(filename_string2.c_str());
  file << "The corresponding csv file for this header is " << filename_string << endl;
  file << "the file contains the real time and cpu time for each run of the APC integrator for different segment lengths and node counts." << endl;
  file << "The Real Time is measured in seconds elapsed for the entire run." << endl;
  file << "The CPU Time is measured in seconds elapsed for each thread during the entire run. e.g. two threads running for 1 real second corresonds to 2 CPU seconds." << endl;
  file << "The data was generated on " << date << endl;
  file << "The columns are as follows:" << endl;

  // write the correspponding header file
  file << "Nodes (#), Orbit Segments (#), Full Function Evaluations (#), Relative Hamiltonian Error (dE/E_0), ";
  for (int i = 0; i < n_trials; i++)
  {
    file << "Real Time Trial " << i + 1 << " (secs)";
    file << ", ";
  }
  for (int i = 0; i < n_trials; i++)
  {
    file << "CPU Time Trial " << i + 1 << " (secs)";
    if (i < n_trials - 1)
    {
      file << ", ";
    }
    else
    {
      file << endl;
    }
  }
  file.close();

  return 0;
}