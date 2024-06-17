/*
 *  AUTHORS:          David Stanley (DavidMS4@Illinois.edu)
 *  DATE WRITTEN:     Nov 2023
 * @ Modified by: Your name
 * @ Modified time: 2023-12-24 17:44:44
 * DESCRIPTION:      Set up an Adaptive-Picard-Chebyshev bootstrap trade study
 * REFERENCE:        Woollands, R., and Junkins, J., "Nonlinear Differential Equation Solvers
 *                   via Adaptive Picard-Chebyshev Iteration: Applications in Astrodynamics", JGCD, 2016.
 */

#include <string>
#include <iostream>

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

struct outputs
{
    double TotalFevals;
    double Hmax;
    double PrepFevals;
    double ForPIFevals;
    double AftPIFevals;
    double BootstrapFevals;
};

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}

vector<double> linspace(double a, double b, int n)
{
    vector<double> array;
    double step = (b - a) / (n - 1);

    while (a <= b)
    {
        array.push_back(a);
        a += step;
    }
    return array;
}

outputs run_trade_study(double alt, double d, bool hot_finish, bool DisableBootstrap)
{
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

    double a = C_Rmoon + alt;
    double e = 0.0;
    double i = 0.0;
    double raan = 0.0;
    double aop = 0.0;
    double ta = 0.0;
    vector<vector<double>> states = elms2rv(a, e, i, raan, aop, ta, C_MU_MOON);

    vector<double> r0 = states[0]; // Initial Position (km)
    vector<double> v0 = states[1];
    double speed = sqrt(pow(v0[0], 2) + pow(v0[1], 2) + pow(v0[2], 2));
    double followtime = d / speed;                     // Follow time (s)
    double T = 2 * C_PI * sqrt(pow(a, 3) / C_MU_MOON); // Orbital period (s)
    double t0 = 0;                                     // initial time (s)
    double tf = t0 + T;                                // final time (s)

    // Orbit orb = SinglePropagate(r0, v0, time_vec,  area,  reflectance,  mass,  drag_C,  compute_drag,  compute_SRP,  compute_third_body);
    Orbit orbit("MOON", "MOON_PA", "J2000");
    orbit.SetProperties(area, reflectance, mass, drag_C);
    orbit.SetPosition0(r0);
    orbit.SetVelocity0(v0);
    orbit.SetIntegrationTime(t0, tf);
    orbit.SetComputeHamiltonian();
    orbit.SetMaxDegree(200);
    orbit.SetTolerance(1.1e-15);
    BootstrapOrbit bootstrap(orbit, followtime);
    if(DisableBootstrap)
    {
        bootstrap.DisableBootstrap();
    }
    

    bootstrap.SetBootstrapHotFinish(hot_finish);
    bootstrap.BootstrapPropagate();

    //get fevals as catergories
    double partialFevalRatio = pow(bootstrap.lowDeg,2)/pow(bootstrap.deg,2);
    double PrepFevals = bootstrap.Feval.Prepare[0]+bootstrap.Feval.Prepare[1]*partialFevalRatio;
    double ForPIFevals = bootstrap.forOrbit.Feval.PicardIteration[0]+bootstrap.forOrbit.Feval.PicardIteration[1]*partialFevalRatio;
    double AftPIFevals = bootstrap.aftOrbit.Feval.PicardIteration[0]+bootstrap.aftOrbit.Feval.PicardIteration[1]*partialFevalRatio;
    double BootstrapFevals = bootstrap.Feval.Bootstrap[0]+bootstrap.Feval.Bootstrap[1]*partialFevalRatio+ bootstrap.Feval.PicardIteration[0]+bootstrap.Feval.PicardIteration[1]*partialFevalRatio;

    outputs output = {bootstrap.TotalFuncEvals, bootstrap.dHmax, PrepFevals, ForPIFevals, AftPIFevals, BootstrapFevals};
    return output;
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

    vector<double> alts = {30,120,480,960};
    // vector<double> alts = {960};
    vector<double> spacings = linspace(0, 0.5, 51);
    // vector<double> spacings = {0.26};
    map<pair<double, double>, outputs> DataBootstrap;
    map<pair<double, double>, outputs> DataFinish;
    map<pair<double, double>, outputs> DataDisabled;
    int counter = 0;
    printProgress(0.0);
    for (double alt : alts)
    {
        for (double d : spacings)
        {
            // cout << "Running alt: " << alt << " spacing: " << d << endl;
            outputs output = run_trade_study(alt, d, false, false);
            DataBootstrap.insert(make_pair(make_pair(alt, d), output));
            output = run_trade_study(alt, d, true, false);
            DataFinish.insert(make_pair(make_pair(alt, d), output));
            output = run_trade_study(alt, d, false, true);
            DataDisabled.insert(make_pair(make_pair(alt, d), output));
            counter++;
            printProgress((double)counter / (alts.size() * spacings.size()));
        }
    }
    // write Data to csv
    string headers = "Altitude,Spacing,Fevals,Hmax,Prepare Segmentation,Forward Satellite,Aftward Satellite,Bootstrapped Satellite";
    string filename = "bootstrap_tradestudy_short960.csv";
    ofstream myfile;
    myfile.open(filename);
    myfile << headers<<"\n";
    myfile << std::setprecision(15);
    for (auto const &x : DataBootstrap)
    {
        myfile << x.first.first << "," << x.first.second << "," << x.second.TotalFevals << "," << x.second.Hmax << "," << x.second.PrepFevals << "," << x.second.ForPIFevals << "," << x.second.AftPIFevals << ","  << x.second.BootstrapFevals << "\n";
    }
    myfile.close();
    filename = "bootstrap_finish_tradestudy_short960.csv";
    myfile.open(filename);
    myfile << headers<<"\n";
    myfile << std::setprecision(15);
    for (auto const &x : DataFinish)
    {
        myfile << x.first.first << "," << x.first.second << "," << x.second.TotalFevals << "," << x.second.Hmax << "," << x.second.PrepFevals << "," << x.second.ForPIFevals << "," << x.second.AftPIFevals << "," << x.second.BootstrapFevals << "\n";
    }
    myfile.close();
    filename = "bootstrap_disabled_tradestudy_short960.csv";
    myfile.open(filename);
    myfile << headers<<"\n";
    myfile << std::setprecision(15);
    for (auto const &x : DataDisabled)
    {
        myfile << x.first.first << "," << x.first.second << "," << x.second.TotalFevals << "," << x.second.Hmax << "," << x.second.PrepFevals << "," << x.second.ForPIFevals << "," << x.second.AftPIFevals << "," << x.second.BootstrapFevals << "\n";
    }
    myfile.close();

    return 0;
}