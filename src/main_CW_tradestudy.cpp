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
#include "clohessywiltshire.h"
using namespace std;
struct state_store
{
    vector<double> x;
    vector<double> y;
    vector<double> z;
    vector<double> vx;
    vector<double> vy;
    vector<double> vz;
    vector<double> t;
    vector<double> r0;
    vector<double> v0;
};

struct outputs
{
    double TotalFevals;
    double Hmax;
    double PrepFevals;
    double ForPIFevals;
    double AftPIFevals;
    double BootstrapFevals;
    state_store chase_state;
    state_store target_state;
    state_store bootstrap_state;
    state_store CW_state;
    state_store zonal_state;
    vector<double> bootstrap_error;
    vector<double> CW_error;
    vector<double> zonal_error;
    vector<double> bootstrap_Hamiltonian;
    vector<double> CW_Hamiltonian;
    vector<double> zonal_Hamiltonian;
};

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

void printProgress(double percentage)
{
    int val = (int)(percentage * 100);
    int lpad = (int)(percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}

vector<double> linspace(double a, double b, int n)
{
    vector<double> array;
    double step = (b - a) / (n - 1);
    for (int i = 0; i < n; i++)
    {
        array.push_back(a + (step * i));
    }

    return array;
}

outputs run_trade_study_bootstrap(double alt, double d)
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
    bool compute_hamiltonian = true; // whether or not the hamiltonian should be computed for the output
    // Ephemeris
    string spk = "de440.bsp";
    string lsk = "naif0012.tls";
    string center = "EARTH";
    string frame = "J2000";

    double a = C_Req + alt;
    double e = 0.0;
    double i = 0.0;
    double raan = 0.0;
    double aop = 0.0;
    double ta = 0.0;
    vector<vector<double>> states = elms2rv(a, e, i, raan, aop, ta, C_MU_EARTH);

    vector<double> r0 = states[0]; // Initial Position (km)
    vector<double> v0 = states[1];
    double speed = sqrt(pow(v0[0], 2) + pow(v0[1], 2) + pow(v0[2], 2));
    double followtime = d / speed;                      // Follow time (s)
    double T = 2 * C_PI * sqrt(pow(a, 3) / C_MU_EARTH); // Orbital period (s)
    double t0 = 0;                                      // initial time (s)
    double tf = t0 + T;                                 // final time (s)

    // Orbit orb = SinglePropagate(r0, v0, time_vec,  area,  reflectance,  mass,  drag_C,  compute_drag,  compute_SRP,  compute_third_body);
    Orbit orbit("EARTH", "EARTH_IAU", "J2000");
    orbit.SetProperties(area, reflectance, mass, drag_C);
    orbit.SetPosition0(r0);
    orbit.SetVelocity0(v0);
    orbit.SetIntegrationTime(t0, tf);
    orbit.SetComputeHamiltonian();
    orbit.SetMaxDegree(200);
    orbit.SetTolerance(1.1e-15);
    BootstrapOrbit chaser(orbit, followtime);
    chaser.DisableBootstrap();
    chaser.SetBootstrapHotFinish(false);
    chaser.BootstrapPropagate();

    // store bootstrap states
    state_store chase_state; // true bootstrap state
    chase_state.t = chaser._Integrator_T;
    chase_state.x = chaser._X;
    chase_state.y = chaser._Y;
    chase_state.z = chaser._Z;
    chase_state.vx = chaser._Vx;
    chase_state.vy = chaser._Vy;
    chase_state.vz = chaser._Vz;
    chase_state.r0 = chaser._In_r0;
    chase_state.v0 = chaser._In_v0;

    state_store target_state; // target state
    target_state.t = chaser.forOrbit._Integrator_T;
    target_state.x = chaser.forOrbit._X;
    target_state.y = chaser.forOrbit._Y;
    target_state.z = chaser.forOrbit._Z;
    target_state.vx = chaser.forOrbit._Vx;
    target_state.vy = chaser.forOrbit._Vy;
    target_state.vz = chaser.forOrbit._Vz;
    target_state.r0 = chaser.forOrbit._In_r0;
    target_state.v0 = chaser.forOrbit._In_v0;


    state_store cw_state;
    vector<double> dS_ECI {chase_state.r0[0] - target_state.r0[0], chase_state.r0[1] - target_state.r0[1], chase_state.r0[2] - target_state.r0[2], chase_state.v0[0] - target_state.v0[0], chase_state.v0[1] - target_state.v0[1], chase_state.v0[2] - target_state.v0[2]};
    vector<double> dr_ECI  {dS_ECI[0], dS_ECI[1], dS_ECI[2]}; // position difference
    vector<double> dv_ECI  {dS_ECI[3], dS_ECI[4], dS_ECI[5]}; // velocity difference

    pair<vector<double>,vector<double>> drAndDv_LVLH = orb2lvlh(target_state.r0,target_state.v0, dr_ECI, dv_ECI);
    vector<double> dr_lvlh = drAndDv_LVLH.first; // position difference in LVLH frame
    vector<double> dv_lvlh = drAndDv_LVLH.second; // velocity difference in LVLH frame
    vector<double> ds_Lvlh = {dr_lvlh[0], dr_lvlh[1], dr_lvlh[2], dv_lvlh[0], dv_lvlh[1], dv_lvlh[2]}; // LVLH state difference
    double n = sqrt(C_MU_EARTH / pow(a, 3)); // mean motion
    for (int i = 0; i < target_state.t.size(); i++)
    {
        double t = target_state.t[i];
        vector<double> CW_state_lvlh = ClohessyWiltshire(t, ds_Lvlh, n);
        vector<double> dr_lvlh_t = {CW_state_lvlh[0], CW_state_lvlh[1], CW_state_lvlh[2]};
        vector<double> dv_lvlh_t = {CW_state_lvlh[3], CW_state_lvlh[4], CW_state_lvlh[5]};
        //get current target state
        vector<double> r_target = {target_state.x[i], target_state.y[i], target_state.z[i]};
        vector<double> v_target = {target_state.vx[i], target_state.vy[i], target_state.vz[i]};
        // convert back to ECI
        pair<vector<double>, vector<double>> rAndV_ECI = lvlh2orb(r_target, v_target, dr_lvlh_t, dv_lvlh_t);
        vector<double> r_ECI = rAndV_ECI.first; // ECI position
        vector<double> v_ECI = rAndV_ECI.second; // ECI velocity
        vector<double> CW_state = {r_ECI[0], r_ECI[1], r_ECI[2], v_ECI[0], v_ECI[1], v_ECI[2]};
        cw_state.t.push_back(t);
        cw_state.x.push_back(CW_state[0] + target_state.x[i]);
        cw_state.y.push_back(CW_state[1] + target_state.y[i]);
        cw_state.z.push_back(CW_state[2] + target_state.z[i]);
        cw_state.vx.push_back(CW_state[3] + target_state.vx[i]);
        cw_state.vy.push_back(CW_state[4] + target_state.vy[i]);
        cw_state.vz.push_back(CW_state[5] + target_state.vz[i]);
    }
    // calculate CW hamiltonian
    vector<double> ts = cw_state.t;
    int len = ts.size();
    double CW_H;
    double CW_H0;
    double CW_dHmax = 0;
    vector<double> CW_dH(len, 0.0);

    for (int i = 0; i < len; i++)
    {
        double s[6] = {cw_state.x[i], cw_state.y[i], cw_state.z[i], cw_state.vx[i], cw_state.vy[i], cw_state.vz[i]};

        jacobiIntegral(ts[i], s, &CW_H, 200);
        if (i == 0)
        {
            CW_H0 = CW_H;
        }
        CW_dH[i] = abs((CW_H - CW_H0) / CW_H0); // Normalized hamiltonian error
        if (CW_dH[i] > CW_dHmax)
        {
            CW_dHmax = CW_dH[i];
        }
    }

    // now rerun the bootstrap without disabling it
    Orbit orbit2("EARTH", "EARTH_IAU", "J2000");
    orbit2.SetProperties(area, reflectance, mass, drag_C);
    orbit2.SetPosition0(r0);
    orbit2.SetVelocity0(v0);
    orbit2.SetIntegrationTime(t0, tf);
    orbit2.SetComputeHamiltonian();
    orbit2.SetMaxDegree(200);
    orbit2.SetTolerance(1.1e-15);
    BootstrapOrbit bootstrap(orbit, followtime);
    bootstrap.SetBootstrapHotFinish(false);
    bootstrap.BootstrapPropagate();

    // get fevals as catergories
    double partialFevalRatio = pow(bootstrap.lowDeg, 2) / pow(bootstrap.deg, 2);
    double PrepFevals = bootstrap.Feval.Prepare[0] + bootstrap.Feval.Prepare[1] * partialFevalRatio;
    double ForPIFevals = bootstrap.forOrbit.Feval.PicardIteration[0] + bootstrap.forOrbit.Feval.PicardIteration[1] * partialFevalRatio;
    double AftPIFevals = bootstrap.aftOrbit.Feval.PicardIteration[0] + bootstrap.aftOrbit.Feval.PicardIteration[1] * partialFevalRatio;
    double BootstrapFevals = bootstrap.Feval.Bootstrap[0] + bootstrap.Feval.Bootstrap[1] * partialFevalRatio + bootstrap.Feval.PicardIteration[0] + bootstrap.Feval.PicardIteration[1] * partialFevalRatio;

    // store the bootstrap states
    state_store bootstrap_state; // true bootstrap state
    bootstrap_state.t = bootstrap._Integrator_T;
    bootstrap_state.x = bootstrap._X;
    bootstrap_state.y = bootstrap._Y;
    bootstrap_state.z = bootstrap._Z;
    bootstrap_state.vx = bootstrap._Vx;
    bootstrap_state.vy = bootstrap._Vy;
    bootstrap_state.vz = bootstrap._Vz;
    bootstrap_state.r0 = bootstrap._In_r0;
    bootstrap_state.v0 = bootstrap._In_v0;

    // run with zonal only
    Orbit orbit3("EARTH", "EARTH_IAU", "J2000");
    orbit3.SetProperties(area, reflectance, mass, drag_C);
    orbit3.SetPosition0(r0);
    orbit3.SetVelocity0(v0);
    orbit3.SetIntegrationTime(t0, tf);
    orbit3.SetComputeHamiltonian();
    orbit3.SetMaxDegree(200);
    orbit3.SetTolerance(1.1e-15);
    BootstrapOrbit zonal_only(orbit3, followtime);
    zonal_only.SetBootstrapHotFinish(false);
    zonal_only.computeZonalOnly();
    zonal_only.DisableBootstrap();
    zonal_only.BootstrapPropagate();

    // store the zonal only states
    state_store zonal_state; // true bootstrap state
    zonal_state.t = zonal_only._Integrator_T;
    zonal_state.x = zonal_only._X;
    zonal_state.y = zonal_only._Y;
    zonal_state.z = zonal_only._Z;
    zonal_state.vx = zonal_only._Vx;
    zonal_state.vy = zonal_only._Vy;
    zonal_state.vz = zonal_only._Vz;
    zonal_state.r0 = zonal_only._In_r0;
    zonal_state.v0 = zonal_only._In_v0;

    // calculate the position error between chaser orbit and the bootstrap orbit and CW orbit
    vector<double> bootstrap_error;
    vector<double> CW_error;
    vector<double> zonal_error;
    for (int i = 0; i < bootstrap_state.t.size(); i++)
    {
        bootstrap_error.push_back(sqrt(pow(chase_state.x[i] - bootstrap_state.x[i], 2) + pow(chase_state.y[i] - bootstrap_state.y[i], 2) + pow(chase_state.z[i] - bootstrap_state.z[i], 2)));
        CW_error.push_back(sqrt(pow(chase_state.x[i] - cw_state.x[i], 2) + pow(chase_state.y[i] - cw_state.y[i], 2) + pow(chase_state.z[i] - cw_state.z[i], 2)));
        zonal_error.push_back(sqrt(pow(chase_state.x[i] - zonal_state.x[i], 2) + pow(chase_state.y[i] - zonal_state.y[i], 2) + pow(chase_state.z[i] - zonal_state.z[i], 2)));
    }

    outputs output = {chaser.TotalFuncEvals, chaser.dHmax, PrepFevals, ForPIFevals, AftPIFevals, BootstrapFevals, chase_state, target_state, bootstrap_state, cw_state, zonal_state, bootstrap_error, CW_error, zonal_error, bootstrap._dH, CW_dH, zonal_only._dH};
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

    vector<double> alts{120, 1200, 12000};
    // vector<double> alts = {960};
    vector<double> spacings = {0.1, 1.0, 10.0, 100.0};
    // vector<double> spacings = {0.26};
    map<pair<double, double>, outputs> DataBootstrap;
    int counter = 0;
    printProgress(0.0);
    for (double alt : alts)
    {
        for (double d : spacings)
        {
            // cout << "Running alt: " << alt << " spacing: " << d << endl;
            outputs output = run_trade_study_bootstrap(alt, d);
            DataBootstrap.insert(make_pair(make_pair(alt, d), output));
            output = run_trade_study_bootstrap(alt, d);
            counter++;
            printProgress((double)counter / (alts.size() * spacings.size()));
        }
    }
    cout << "Writing output to CSV files..." << endl;
    // for each alt spacing pair generate a csv saving for each timestep the chase position and velocity and the three errors
    for (double alt : alts)
    {
        for (double spacing : spacings)
        {
            string filename = "Position_error_alt_" + to_string(alt) + "_spacing_" + to_string(spacing) + ".csv";
            ofstream file(filename);
            file << "Time,Chase_X,Chase_Y,Chase_Z,Target_X,Target_Y,Target_Z,Bootstrap_Error,Bootstrap_Hamiltonian_Error,CW_Error,CW_Hamiltonian_Error,Zonal_Error,Zonal_Hamiltonian_Error\n";
            for (int i = 0; i < DataBootstrap[{alt, spacing}].chase_state.t.size(); i++)
            {
                file << DataBootstrap[{alt, spacing}].chase_state.t[i] << ","
                     << DataBootstrap[{alt, spacing}].chase_state.x[i] << ","
                     << DataBootstrap[{alt, spacing}].chase_state.y[i] << ","
                     << DataBootstrap[{alt, spacing}].chase_state.z[i] << ","
                     << DataBootstrap[{alt, spacing}].target_state.x[i] << ","
                     << DataBootstrap[{alt, spacing}].target_state.y[i] << ","
                     << DataBootstrap[{alt, spacing}].target_state.z[i] << ","
                     << DataBootstrap[{alt, spacing}].bootstrap_error[i] << ","
                     << DataBootstrap[{alt, spacing}].bootstrap_Hamiltonian[i] << ","
                     << DataBootstrap[{alt, spacing}].CW_error[i] << ","
                     << DataBootstrap[{alt, spacing}].CW_Hamiltonian[i] << ","
                     << DataBootstrap[{alt, spacing}].zonal_error[i] << ","
                     << DataBootstrap[{alt, spacing}].zonal_Hamiltonian[i] << "\n";
            }
            file.close();
        }
    }
    cout << "Done!" << endl;
    return 0;
}