/*
 *  AUTHORS:          David Stanley (davidms4@illinois.edu)
 *  DATE WRITTEN:     Feb 2022
 *  LAST MODIFIED:    Feb 2022
 *  AFFILIATION:      Department of Aerospace Engineering, University of Illinois Urbana-Champaign
 *  DESCRIPTION:      Class that stores orbit solution and properties
 *  REFERENCE:        Woollands, R., and Junkins, J., "Nonlinear Differential Equation Solvers
 *                    via Adaptive Picard-Chebyshev Iteration: Applications in Astrodynamics", JGCD, 2016.
 *
 */

#include <string>
#include <vector>
#include <array>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <math.h>
#include "Orbit.h"
#include "const.h"
#include "matrix_loader.h"
#include "SpiceUsr.h"
#include "Ephemeris.hpp"
#include "APC.h"
#include "flags.h"
#include "FandG.h"
#include "polydegree_segments.h"
#include "prepare_propagator.h"
#include "interpolate.h"
#include "c_functions.h"
#include "picard_iteration.h"
#include "perturbed_gravity.h"
#include "ecef2eci.h"
#include "eci2ecef.h"
#include "lunar_perturbed_gravity.h"
#include "perturbations.h"
#include "picard_error_feedback.h"
#include "reosc_perigee.h"
#include "EGM2008.h"
#include "GRGM1200b.h"

using namespace std;
// Orbit Constructors
Orbit::Orbit()
{
    // Empty Orbit Constructor
}
Orbit::Orbit(string primary, string IOframe, string epoch)
{
    // Initialize orbit with primary body and input data frame
    if (InSet(primary, C_VALID_PRIMARIES))
    {
        _primary = primary;
        _validPrimary = true;
        SetMu(_primary);
        _Req = GetPrimaryRadius();
    }
    _IOFrame = IOframe;
    _InertFrame = GetInertFrame();
    _FixedFrame = GetFixedFrame();
    _epoch = epoch;
}
Orbit::Orbit(std::vector<std::vector<double>> Solution)
{
    SetSolution(Solution);
}
Orbit::Orbit(double area, double reflectivity, double mass, double Cd, bool compute_drag, bool compute_SRP, bool compute_third_body, bool compute_hamiltonian, int id)
{
    // Set primary body for two body targets
    _primary = "Earth"; // defaults to orbiting around Earth
    _validPrimary = true;
    // Set orbital properties
    SetProperties(area, reflectivity, mass, Cd, compute_drag, compute_SRP, compute_third_body, compute_hamiltonian, id);
    suborbital = false;
}

Orbit::Orbit(double area, double reflectivity, double mass, double Cd, bool compute_drag, bool compute_SRP, bool compute_third_body, bool compute_hamiltonian, int id, string primary = "Earth")
{
    // Set primary body for two body orbits after checking if valid
    if (InSet(primary, C_VALID_PRIMARIES))
    {
        _primary = primary;
        _validPrimary = true;
        SetMu(_primary);
    }

    SetProperties(area, reflectivity, mass, Cd, compute_drag, compute_SRP, compute_third_body, compute_hamiltonian, id);
    suborbital = false;
}

// Getter Functions
std::vector<std::vector<double>> Orbit::getPosition()
{
    std::vector<std::vector<double>>::const_iterator first = Soln.begin() + 1;
    std::vector<std::vector<double>>::const_iterator last = Soln.begin() + 4;
    std::vector<std::vector<double>> Position(first, last);

    return Position;
}
std::vector<std::vector<double>> Orbit::getVelocity()
{
    std::vector<std::vector<double>>::const_iterator first = Soln.begin() + 4;
    std::vector<std::vector<double>>::const_iterator last = Soln.begin() + 7;
    std::vector<std::vector<double>> Velocity(first, last);

    return Velocity;
}

double Orbit::GetPrimaryRadius()
{
    double Req = 0.0;
    // Earth reference radius
    if (InSet(_primary, {"Earth"}))
    {
        Req = 6378.137;
    }
    // Moon reference radius
    else if (InSet(_primary, {"Moon"}))
    {
        Req = 1.738e+03;
    }
    return Req;
}

string Orbit::GetInertFrame()
{
    string InertFrame;
    if (InSet(_IOFrame, {"MOON_PA", "J2000", "EME2000"}))
    {
        InertFrame = "J2000";
    }
    return InertFrame;
}

string Orbit::GetFixedFrame()
{
    string FixedFrame;
    if (InSet(_primary, {"MOON"}))
    {
        FixedFrame = "MOON_PA";
    }
    else if (InSet(_primary, {"EARTH"}))
    {
        // FIXME:THE IAU EARTH frame is not precies enough.
        FixedFrame = "IAU_EARTH";
    }
    return FixedFrame;
}

string Orbit::GetInertCenter()
{
    string InertCenter;
    if (InSet(_primary, {"MOON"}))
    {
        InertCenter = "MOON";
    }
    else if (InSet(_primary, {"EARTH"}))
    {
        InertCenter = "EARTH";
    }
    return InertCenter;
}

// Setter Functions
void Orbit::SetSolution(std::vector<std::vector<double>> Solution)
{
    // store the solution
    Soln = Solution;
}

void Orbit::SetTimeVec(std::vector<double> time_vec)
{
    // Store user time vector and raise flag indicating that a time vector was given
    T = time_vec;
    USER_TIME = true;
}

void Orbit::SetProperties(double area, double reflectance, double mass, double Cd, bool compute_drag, bool compute_SRP, bool compute_third_body, bool compute_hamiltonian, int id, string primary)
{
    // Set satellite properties and flags for perturbations
    satproperties._Area = area;
    satproperties._Mass = mass;
    satproperties._Reflectance = reflectance;
    satproperties._Cd = Cd;
    // Perturbation flags
    Compute_Drag = compute_drag;
    Compute_SRP = compute_SRP;
    Compute_Third_Body = compute_third_body;
    // Other flags
    Compute_Hamiltonian = compute_hamiltonian;
    // ID
    ID = id;
}

void Orbit::SetSubOrbital()
{
    suborbital = true;
}

void Orbit::SetCC(std::vector<double> A, std::vector<double> B, std::vector<double> W1, std::vector<double> W2, int N, int coeff_size, std::vector<double> seg_times, double TF, double T0, int total_segs)
{
    // Set values for chebyshev coefficients and associated variables
    CC.A = A;
    CC.B = B;
    CC.W1 = W1;
    CC.W2 = W2;
    CC.N = N;
    CC.coeff_size = coeff_size;
    CC.seg_times = seg_times;
    CC.TF = TF;
    CC.T0 = T0;
    CC.total_segs = total_segs;
}

void Orbit::SetPosition0(vector<double> r0)
{
    // set initial position and the position part of initial state
    _In_r0 = r0;
    copy(r0.begin(), r0.end(), _In_s0.begin());
}

void Orbit::SetVelocity0(vector<double> v0)
{
    // set initial velocity and the velocity part of initial state
    _In_v0 = v0;
    copy(v0.begin(), v0.end(), _In_s0.begin() + 3);
}

void Orbit::SetState0(vector<double> s0)
{
    copy(s0.begin(), s0.begin() + 3, _In_r0.begin());
    copy(s0.begin() + 3, s0.end(), _In_v0.begin());
    // set state
    _In_s0 = s0;
}

void Orbit::SetMu(string primary)
{
    // lower primary name case
    std::transform(primary.begin(), primary.end(), primary.begin(),
                   [](unsigned char c)
                   { return std::tolower(c); });
    // Assign mu based on name of primary body. Defined in const.h
    if (primary.compare("earth") == 0)
    {
        _mu = C_MU_EARTH;
    }
    else if (primary.compare("moon") == 0)
    {
        _mu = C_MU_MOON;
    }
}

bool Orbit::InSet(string item, vector<string> validset)
{
    bool inSet = false;
    // convert item to lower case for consistency
    std::transform(item.begin(), item.end(), item.begin(),
                   [](unsigned char c)
                   { return std::tolower(c); });
    // Iterate through the valid set
    for (string validItem : validset)
    {
        string lowerValid = validItem;
        // convert valid item to lowercase for consistency
        std::transform(lowerValid.begin(), lowerValid.end(), lowerValid.begin(), [](unsigned char c)
                       { return std::tolower(c); });
        // Check if item is the same as the validItems
        if (item.compare(lowerValid) == 0)
        {
            // Found in valid set, end the loop
            inSet = true;
            break;
        }
    }
    return inSet;
}

void Orbit::SetIntegrationTime(double t0, double tf)
{
    _Ephemeris_t0 = t0;
    _Ephemeris_tf = tf;
    _Integrator_t0 = 0;
    _Integrator_tf = tf - t0;
}

void Orbit::SetIntegrationTime(double tf)
{
    // If only a final time is given. Assume the start time is at 0 s in J2000
    SetIntegrationTime(0, tf);
}

void Orbit::SetIntegrationTime(string date0, string datef)
{
    // convert date strings to char array for SPICE
    const char *char_date0 = date0.c_str();
    const char *char_datef = datef.c_str();
    double t0;
    double tf;
    // convert date to seconds in J2000;
    str2et_c(char_date0, &t0);
    str2et_c(char_datef, &tf);
    // set times with the set function that takes doubles.
    SetIntegrationTime(t0, tf);
}

void Orbit::SetProperties(double area, double reflectance, double mass, double Cd)
{
    satproperties._Area = area;
    satproperties._Reflectance = reflectance;
    satproperties._Mass = mass;
    satproperties._Cd = Cd;
}

bool Orbit::HasAtmosphere()
{
    // List of bodies with atmospheric models in APC
    vector<string> atmospheric_bodies = {"Earth"};
    bool hasAtmosphere = false;

    // If the orbit primary body has an atmospher then return true
    if (InSet(_primary, atmospheric_bodies))
    {
        hasAtmosphere = true;
    }
    return hasAtmosphere;
}

void Orbit::SinglePropagate()
{
    double t0 = _Integrator_t0;
    double tf = _Integrator_tf;
    vector<double> r = _In_r0;
    vector<double> v = _In_v0;
    EphemerisManager ephem = cacheEphemeris(t0, tf + 3600);
    // Initialize orbit object
    double dt = 30;
    double len = int(ceil(tf / dt));
    PropagateOrbit(r, v, t0, tf, *this, ephem);
}

void Orbit::PrintConfig()
{
    bool Kernels_were_loaded = g_KERNELS_LOADED;
    LoadKernels();
    // Print Header
    cout << "===Orbit Information Start===" << endl;
    // Print pimary body and frame
    cout << "Primary Body: " << GetPrimaryBody() << endl;
    cout << "Fixed Frame: " << GetFixedFrame() << endl;
    // Print Integration start and end time and length
    ConstSpiceChar *format = "C";
    SpiceInt prec = 0;
    SpiceInt utclen = 21;
    SpiceChar start_date[utclen];
    SpiceChar end_date[utclen];
    et2utc_c(_Ephemeris_t0, format, prec, utclen, start_date);
    et2utc_c(_Ephemeris_tf, format, prec, utclen, end_date);
    cout << "Start Time: " << start_date << endl;
    cout << "End Time: " << end_date << endl;
    cout << "Propagation Time: " << _Ephemeris_tf - _Ephemeris_t0 << " (s)" << endl;
    cout << "===Orbit Information End===" << endl;

    // Unload kernels only if kernels were not loaded before this method was run
    if (not Kernels_were_loaded)
    {
        UnloadKernels();
    }
}

void Orbit::Interpolate()
{
    // Interpolate with time vector generated by APC or user specified time vector
    if (!USER_TIME)
    {
        _Integrator_T = DefaultIntegratorTimeVec();
    }
    else
    {
        _Integrator_T = UserIntegratorTimevec();
    }

    // access relevant variable pointers from orbit
    std::vector<double> &ALPHA = CC.A;
    std::vector<double> &BETA = CC.B;
    std::vector<double> &seg_times = segment_end_times;
    std::vector<double> &time_out = _Integrator_T;
    vector<double> &X_out = _X;
    vector<double> &Y_out = _Y;
    vector<double> &Z_out = _Z;
    vector<double> &Vx_out = _Vx;
    vector<double> &Vy_out = _Vy;
    vector<double> &Vz_out = _Vz;
    vector<vector<double>> &y = _y;
    int prev_cnt = 0;

    // User specified output times
    int len;
    len = time_out.size();
    int soln_size = time_out.size();
    std::vector<double> Soln(soln_size * 6, 0.0);
    // individual coordiante vectors
    X_out.resize(soln_size, 0.0);
    Y_out.resize(soln_size, 0.0);
    Z_out.resize(soln_size, 0.0);
    Vx_out.resize(soln_size, 0.0);
    Vy_out.resize(soln_size, 0.0);
    Vz_out.resize(soln_size, 0.0);
    y.resize(soln_size, vector<double>(6, 0.0));
    double test_time = 0.0;
    // Loop through all segments
    for (int i = 1; i <= total_segs; i++)
    {
        int sz = count_if(time_out.begin(), time_out.end(), [&](int timestep)
                          { return (timestep >= seg_times[i - 1] && timestep <= seg_times[i]); });
        // Initialization
        std::vector<double> Beta(N * 3);
        std::vector<double> Alpha((N + 1) * 3);
        std::vector<double> tt(sz);
        std::vector<double> tau(sz + 1);

        double w1, w2;
        w1 = W1[i - 1];
        w2 = W2[i - 1];

        // User desired times for a given segment
        int cnt = 0;
        // printf("seg_times %f\t%f\n",seg_times[i-1],seg_times[i]);
        for (int j = 0; j < len; j++)
        {
            if (time_out[j] == seg_times[i - 1])
            {
                // printf("t1 %f\n",time_out[j]);
                tau[cnt] = -1.0;
                cnt = cnt + 1; // Number of times steps per segment
            }
            if (time_out[j] == seg_times[i])
            {
                // printf("t2 %f\n",time_out[j]);
                tau[cnt] = 1.0;
                cnt = cnt + 1;
            }
            if (time_out[j] > seg_times[i - 1] && time_out[j] < seg_times[i])
            {
                tt[cnt] = time_out[j];
                tau[cnt] = (tt[cnt] - w1) / w2;
                cnt = cnt + 1; // Number of times steps per segment
            }
        }

        // Chebyshev Velocity & Position Matrices
        std::vector<double> Tv(cnt * N);
        std::vector<double> Tp(cnt * (N + 1));

        for (int t = 1; t <= cnt; t++)
        {
            for (int kk = 0; kk <= N - 1; kk++)
            {
                // Velocity
                Tv[ID2(t, kk + 1, cnt)] = cos(kk * acos(tau[t - 1]));
            }
            for (int kk = 0; kk <= N; kk++)
            {
                // Position
                Tp[ID2(t, kk + 1, cnt)] = cos(kk * acos(tau[t - 1]));
            }
        }

        // Velocity Coefficients for a Segment
        for (int p = 1; p <= N; p++)
        {
            Beta[ID2(p, 1, N)] = BETA[ID2(p + ((i - 1) * N), 1, coeff_size)];
            Beta[ID2(p, 2, N)] = BETA[ID2(p + ((i - 1) * N), 2, coeff_size)];
            Beta[ID2(p, 3, N)] = BETA[ID2(p + ((i - 1) * N), 3, coeff_size)];
        }
        std::vector<double> v_interp;
        v_interp = matmul(Tv, Beta, cnt, N, 3, cnt, N);

        // sanity check
        int check = ID2(cnt + prev_cnt, 6, soln_size);
        if (check >= soln_size * 6)
        {
            std::cout << "exceeding soln size\n";
        }

        // Velocity
        for (int p = 1; p <= cnt; p++)
        {
            y[p - 1 + prev_cnt][3] = v_interp[ID2(p, 1, cnt)];
            y[p - 1 + prev_cnt][4] = v_interp[ID2(p, 2, cnt)];
            y[p - 1 + prev_cnt][5] = v_interp[ID2(p, 3, cnt)];
            Vx_out[p - 1 + prev_cnt] = v_interp[ID2(p, 1, cnt)];
            Vy_out[p - 1 + prev_cnt] = v_interp[ID2(p, 2, cnt)];
            Vz_out[p - 1 + prev_cnt] = v_interp[ID2(p, 3, cnt)];
        }

        // Position Coefficients for a Segment
        for (int p = 1; p <= N + 1; p++)
        {
            Alpha[ID2(p, 1, N + 1)] = ALPHA[ID2(p + ((i - 1) * (N + 1)), 1, coeff_size)];
            Alpha[ID2(p, 2, N + 1)] = ALPHA[ID2(p + ((i - 1) * (N + 1)), 2, coeff_size)];
            Alpha[ID2(p, 3, N + 1)] = ALPHA[ID2(p + ((i - 1) * (N + 1)), 3, coeff_size)];
        }
        std::vector<double> x_interp;
        x_interp = matmul(Tp, Alpha, cnt, N + 1, 3, cnt, N + 1);
        // Position
        for (int p = 1; p <= cnt; p++)
        {
            y[p - 1 + prev_cnt][0] = x_interp[ID2(p, 1, cnt)];
            y[p - 1 + prev_cnt][1] = x_interp[ID2(p, 2, cnt)];
            y[p - 1 + prev_cnt][2] = x_interp[ID2(p, 3, cnt)];
            X_out[p - 1 + prev_cnt] = x_interp[ID2(p, 1, cnt)];
            Y_out[p - 1 + prev_cnt] = x_interp[ID2(p, 2, cnt)];
            Z_out[p - 1 + prev_cnt] = x_interp[ID2(p, 3, cnt)];
            // printf("Soln %f\t%f\t%f\n",Soln[ID2(p+prev_cnt,1,soln_size)],Soln[ID2(p+prev_cnt,2,soln_size)],Soln[ID2(p+prev_cnt,3,soln_size)]);
        }
        // printf("tfdt %f\t%f\n",tf,dt);
        prev_cnt = prev_cnt + cnt; // Counter to track position in Soln array
    }
    return;
}

vector<double> Orbit::UserIntegratorTimevec()
{
    if (USER_TIME)
    {
        // generate user time vector for integration where t0=0;
        int len = T.size();
        vector<double> time_vec(len, 0.0);
        for (int jj = 0; jj < len; jj++)
        {
            time_vec[jj] = T[jj] - _Ephemeris_t0;
        }
        return time_vec;
    }
    else
    {
        throw invalid_argument("User time vector not set");
    }
}

vector<double> Orbit::DefaultIntegratorTimeVec()
{
    // Integrator always uses t0=0.
    // Generate default time vector for interpolation
    double t0 = _Integrator_t0;
    double tf = _Integrator_tf;
    double dt = _dt;
    double len = int(ceil((tf - t0) / dt));
    vector<double> time_vec(len, 0.0);
    for (int jj = 0; jj < len; jj++)
    {
        double time = jj * dt;
        if (time > tf)
        {
            time = tf;
        }
        time_vec[jj] = time;
    }
    return time_vec;
}

vector<double> Orbit::DefaultTimeVec()
{
    // Generate default time vector for output
    double t0 = _Ephemeris_t0;
    double tf = _Ephemeris_tf;
    double dt = _dt;
    double len = int(ceil(tf / dt));
    vector<double> time_vec(len, 0.0);
    for (int jj = 0; jj < len; jj++)
    {
        double time = jj * dt;
        if (time > tf)
        {
            time = tf;
        }
        time_vec[jj] = time;
    }
    return time_vec;
}

vector<double> Orbit::GetTime()
{
    if (USER_TIME)
    {
        return T;
    }
    else
    {
        return DefaultTimeVec();
    }
}

void Orbit::HamiltonianCheck()
{
    int len = _y.size();
    double H;
    double H0;
    dHmax = 0;
    _dH.resize(len, 0.0);
    vector<double> ts = GetTime();
    for (int i = 0; i < len; i++)
    {

        jacobiIntegral_GRGM1200b(ts[i], &_y[i][0], &H, deg, *this);
        if (i == 0)
        {
            H0 = H;
        }
        _dH[i] = abs((H - H0) / H0); // Normalized hamiltonian error
        if (_dH[i] > dHmax)
        {
            dHmax = _dH[i];
        }
    }
    if (g_VERBOSE)
    {
        cout << name << ": Maximum Hamiltonian Error: " << dHmax << endl;
    }
    return;
}

// Bootstrap Orbit

BootstrapOrbit::BootstrapOrbit(const Orbit &orbit, double followTime) : forOrbit(_primary, _IOFrame, _epoch), aftOrbit(_primary, _IOFrame, _epoch), Orbit(orbit)
{
    SetName("Bootstrap Orbit");
    // Initialize orbit objects list
    orbitslist = {&forOrbit, &aftOrbit, this};
    // Initialize orbit object initial conditions
    double z0[6] = {_In_r0[0], _In_r0[1], _In_r0[2], _In_v0[0], _In_v0[1], _In_v0[2]};
    double zfor[6] = {0.0};
    double zaft[6] = {0.0};
    // Generate the state ahead and behind of the target orbit
    FandG(z0, zfor, followTime, _mu);
    FandG(z0, zaft, -followTime, _mu);
    // assign values to the forward orbit
    forOrbit.SetName("Forward Orbit");
    forOrbit.SetPosition0({zfor[0], zfor[1], zfor[2]});
    forOrbit.SetVelocity0({zfor[3], zfor[4], zfor[5]});
    forOrbit.SetIntegrationTime(_Ephemeris_t0, _Ephemeris_tf);
    forOrbit.SetProperties(satproperties._Area, satproperties._Reflectance, satproperties._Mass, satproperties._Cd);
    forOrbit.SetComputeHamiltonian(Compute_Hamiltonian);
    forOrbit.SetComputeSRP(Compute_SRP);
    forOrbit.SetComputeThirdBody(Compute_Third_Body);
    forOrbit.SetComputeDrag(Compute_Drag);
    forOrbit.deg = this->deg;
    // assign values to the aft orbit
    aftOrbit.SetName("Aft Orbit");
    aftOrbit.SetPosition0({zaft[0], zaft[1], zaft[2]});
    aftOrbit.SetVelocity0({zaft[3], zaft[4], zaft[5]});
    aftOrbit.SetIntegrationTime(_Ephemeris_t0, _Ephemeris_tf);
    aftOrbit.SetProperties(satproperties._Area, satproperties._Reflectance, satproperties._Mass, satproperties._Cd);
    aftOrbit.SetComputeHamiltonian(Compute_Hamiltonian);
    aftOrbit.SetComputeSRP(Compute_SRP);
    aftOrbit.SetComputeThirdBody(Compute_Third_Body);
    aftOrbit.SetComputeDrag(Compute_Drag);
    aftOrbit.deg = this->deg;
}

BootstrapOrbit::BootstrapOrbit()
{
}

void BootstrapOrbit::BootstrapPropagate()
{
    double t0 = _Integrator_t0;
    double tf = _Integrator_tf;
    EphemerisManager ephem = cacheEphemeris(t0, tf + 3600);
    // Initialize orbit object
    double *r0 = &_In_r0[0];
    double *v0 = &_In_v0[0];
    Bootstrap_Adaptive_Picard_Chebyshev(ephem);
    double FullEvals = 0; // Total number of full function evaluations
    double PartialEvals = 0;    // Total number of partial function evaluations
    //Sum up all function evaluations
    for(Orbit *orbitRef : orbitslist)
    {
        FullEvals += orbitRef->Feval.Prepare[0]+orbitRef->Feval.PicardIteration[0]+orbitRef->Feval.Bootstrap[0];
        PartialEvals += orbitRef->Feval.Prepare[1]+orbitRef->Feval.PicardIteration[1]+orbitRef->Feval.Bootstrap[1];
    }
    double partialRatio = pow(lowDeg, 2) / pow(deg, 2);
    TotalFuncEvals = FullEvals + PartialEvals* partialRatio;
    if (g_VERBOSE)
    {
        cout << "Total Function Evaluations: " << TotalFuncEvals << endl;
    }
    if (Compute_Hamiltonian)
    {
        for (Orbit *orbitRef : orbitslist)
        {
            orbitRef->HamiltonianCheck();
        }
    }

    return;
}

void BootstrapOrbit::Bootstrap_Adaptive_Picard_Chebyshev(class EphemerisManager ephem)
{
    // 1. DETERMINE DEGREE/SEGMENTATION SCHEME
    // Compute the polynomial degree and number of segments per orbit that will
    // result in a solution that satisfies the user specified tolerance.
    double* Feval = &this->Feval.Prepare[0];
    polydegree_segments(*this, Feval);
    // copy scheme to forward and aft orbits
    aftOrbit.seg = seg;
    aftOrbit.N = N;
    aftOrbit.tp = tp;
    forOrbit.seg = seg;
    forOrbit.N = N;
    forOrbit.tp = tp;
    // Array size for coefficients and solution
    double &tf = _Integrator_tf;
    coeff_size = int((tf / Period + 1.0) * (seg + 2.0) * (N + 1));
    aftOrbit.coeff_size = coeff_size;
    forOrbit.coeff_size = coeff_size;

    // 2. PREPARE PROPAGATOR
    // Compute and store the begin and end times for each segment (based on true
    // anomaly segmentation) and load the constant matrices corresponding to N.
    prep_HS = -1; // Hot start switch condition
    aftOrbit.prep_HS = prep_HS;
    forOrbit.prep_HS = prep_HS;
    prepare_propagator(*this);
    // FIXME: the quadrature matrices don't really need to be copied but it is easier to do this so that other code doesn't need to be modified.
    aftOrbit.T2 = T2; // [(M+1)x(N+1)]
    aftOrbit.P2 = P2; // [(N+1)xN]
    aftOrbit.T1 = T1; // [(M+1)xN]
    aftOrbit.P1 = P1; // [Nx(N-1)]
    aftOrbit.Ta = Ta; // [(M+1)x(N-1)]
    aftOrbit.A = A;   // [(N-1)x(M+1)]
    aftOrbit.t_orig = t_orig;
    aftOrbit.tvec = tvec;
    forOrbit.T2 = T2; // [(M+1)x(N+1)]
    forOrbit.P2 = P2; // [(N+1)xN]
    forOrbit.T1 = T1; // [(M+1)xN]
    forOrbit.P1 = P1; // [Nx(N-1)]
    forOrbit.Ta = Ta; // [(M+1)x(N-1)]
    forOrbit.A = A;   // [(N-1)x(M+1)]
    forOrbit.t_orig = t_orig;
    forOrbit.tvec = tvec;

    // 3. PICARD-CHEBYSHEV PROPAGATOR
    // Propagate from t0 to tf, iterating on each segment (Picard Iteration), until
    // completion. Stores solution as chebyshev nodes for each segment
    M = N; // Number of Chebyshev nodes
    CC.A.resize((coeff_size * 3), 0.0);
    CC.B.resize((coeff_size * 3), 0.0);
    CC.total_segs = 0; // running count of orbital path segments
    aftOrbit.M = M;
    aftOrbit.CC.A.resize((coeff_size * 3), 0.0);
    aftOrbit.CC.B.resize((coeff_size * 3), 0.0);
    aftOrbit.CC.total_segs = 0; // running count of orbital path segments
    forOrbit.M = M;
    forOrbit.CC.A.resize((coeff_size * 3), 0.0);
    forOrbit.CC.B.resize((coeff_size * 3), 0.0);
    forOrbit.CC.total_segs = 0; // running count of orbital path segments
    // reserve space in vectors for interpolation and solution
    int sz = int(ceil(tf / Period) * (seg + 1)); // ensure sufficient space by overestimating the number of segments
    segment_end_times.resize(sz, 0.0);
    W1.resize(sz, 0.0);
    W2.resize(sz, 0.0);
    aftOrbit.segment_end_times.resize(sz, 0.0);
    aftOrbit.W1.resize(sz, 0.0);
    aftOrbit.W2.resize(sz, 0.0);
    forOrbit.segment_end_times.resize(sz, 0.0);
    forOrbit.W1.resize(sz, 0.0);
    forOrbit.W2.resize(sz, 0.0);
    Bootstrap_Picard_Chebyshev_Propagator(ephem);

    // 4. INTERPOLATE SOLUTION
    // The Chebyshev coefficients from each of the orbit segments are used to compute
    // the solution (position & velocity) at the user specified times.
    for (Orbit *orbit : orbitslist)
    {
        // run interpolation using the output time vector for each orbit object
        orbit->segment_end_times = segment_end_times;
        orbit->total_segs = total_segs;
        orbit->Interpolate();
    }
    return;
}

void BootstrapOrbit::Bootstrap_Picard_Chebyshev_Propagator(EphemerisManager ephem)
{

    if (g_DEBUG_PICARD)
    {
        cout << "Entering Bootstrap_Picard_chebyshev_Propagator" << std::endl;
    }

    int loop = 0;              // Break loop condition
    k = 0;                     // Counter: segments per orbit
    int &seg_cnt = total_segs; // Counter: total segments
    double mu = GetPrimaryGravitationalParameter();
    double w1, w2, tf;
    vector<double> &segment_times = segment_end_times;

    // Set segment initial conditions
    for (int i = 0; i < 3; i++)
    {
        r0_seg[i] = _In_r0[i];
        v0_seg[i] = _In_v0[i];
        aftOrbit.r0_seg[i] = aftOrbit._In_r0[i];
        aftOrbit.v0_seg[i] = aftOrbit._In_v0[i];
        forOrbit.r0_seg[i] = forOrbit._In_r0[i];
        forOrbit.v0_seg[i] = forOrbit._In_v0[i];
    }

    // PROPAGATION LOOP
    while (loop == 0)
    {
        if (g_DEBUG_PICARD)
        {
            cout << "===" << std::endl;
            cout << "Starting Segment " << seg_cnt + 1 << "." << std::endl;
            double debug_altitude = 0.0;
            double debug_velocity = 0.0;
            for (int j = 0; j < 3; j++)
            {
                debug_altitude += r0_seg[j] * r0_seg[j];
                debug_velocity += v0_seg[j] * v0_seg[j];
            }
            debug_altitude = sqrt(debug_altitude);
            debug_velocity = sqrt(debug_velocity);
            cout << "debug_altitude: " << debug_altitude << std::endl;
            cout << "debug_velocity: " << debug_velocity << std::endl;
        }

        // A. CALCULATE CHEBYSHEV NODES
        //  Compute cosine time vector for the current segment
        double t0 = tvec[k];
        tf = tvec[k + 1];
        while (tf == 0.0)
        {
            k++;
            t0 = tvec[k];
            tf = tvec[k + 1];
        }
        if (tf > _Integrator_tf)
        {
            tf = _Integrator_tf;
        }
        w1 = (tf + t0) / 2.0;
        w2 = (tf - t0) / 2.0;
        W1[seg_cnt] = w1;
        W2[seg_cnt] = w2;
        aftOrbit.W1[seg_cnt] = w1;
        aftOrbit.W2[seg_cnt] = w2;
        forOrbit.W1[seg_cnt] = w1;
        forOrbit.W2[seg_cnt] = w2;

        //  Compute Chebyshev nodes for the current segment
        for (Orbit *orbitRef : orbitslist)
        {
            Orbit &orbit = *orbitRef;
            orbit.tau.resize(M + 1, 0.0);
            orbit.times_seg.resize(M + 1, 0.0);
            orbit.X_seg.resize((M + 1) * 3, 0.0);
            orbit.V_seg.resize((M + 1) * 3, 0.0);
            orbit.Beta_seg.resize(N * 3, 0.0);
            orbit.Alpha_seg.resize((N + 1) * 3, 0.0);

            // Keplerian warm start
            // Initial state
            for (int j = 0; j < 3; j++)
            {
                orbit.z0_seg[j] = orbit.r0_seg[j];
                orbit.z0_seg[j + 3] = orbit.v0_seg[j];
            }
            double z[6] = {0.0};
            for (int cnt = 0; cnt <= M; cnt++)
            {
                // compute cosine node
                orbit.tau[cnt] = -cos(C_PI * cnt / M);
                orbit.times_seg[cnt] = orbit.tau[cnt] * w2 + w1;
                // compute F and G state at node
                FandG(orbit.z0_seg, z, orbit.times_seg[cnt] - t0, mu);
                // store initial segment guess
                orbit.X_seg[ID2(cnt + 1, 1, M + 1)] = z[0];
                orbit.X_seg[ID2(cnt + 1, 2, M + 1)] = z[1];
                orbit.X_seg[ID2(cnt + 1, 3, M + 1)] = z[2];
                orbit.V_seg[ID2(cnt + 1, 1, M + 1)] = z[3];
                orbit.V_seg[ID2(cnt + 1, 2, M + 1)] = z[4];
                orbit.V_seg[ID2(cnt + 1, 3, M + 1)] = z[5];
            }
            // End orbit list iteration
        }
        // Perform picard iteration for each orbit and move on when for and aft orbits converge
        Bootstrap_Picard_Iteration(ephem);

        //  STORE TRAJECTORY COEFFICIENTS
        for (Orbit *orbitRef : orbitslist)
        {
            std::vector<double> &BETA = orbitRef->CC.B;
            std::vector<double> &ALPHA = orbitRef->CC.A;
            std::vector<double> &Beta = orbitRef->Beta_seg;
            std::vector<double> &Alpha = orbitRef->Alpha_seg;
            for (int i = 1; i <= N; i++)
            {
                BETA[ID2(i + (seg_cnt * N), 1, coeff_size)] = Beta[ID2(i, 1, N)];
                BETA[ID2(i + (seg_cnt * N), 2, coeff_size)] = Beta[ID2(i, 2, N)];
                BETA[ID2(i + (seg_cnt * N), 3, coeff_size)] = Beta[ID2(i, 3, N)];
            }
            for (int i = 1; i <= N + 1; i++)
            {
                // Store X and V points
                ALPHA[ID2(i + seg_cnt * (N + 1), 1, coeff_size)] = Alpha[ID2(i, 1, N + 1)];
                ALPHA[ID2(i + seg_cnt * (N + 1), 2, coeff_size)] = Alpha[ID2(i, 2, N + 1)];
                ALPHA[ID2(i + seg_cnt * (N + 1), 3, coeff_size)] = Alpha[ID2(i, 3, N + 1)];
            }
        }

        // ASSIGN NEXT SEGMENT ICs
        //  Assign the initial conditions for the next segment
        for (Orbit *orbitRef : orbitslist)
        {
            Orbit &orbit = *orbitRef;
            for (int j = 0; j < 3; j++)
            {
                orbit.r0_seg[j] = orbit.X_seg[ID2(M + 1, j + 1, M + 1)];
                orbit.v0_seg[j] = orbit.V_seg[ID2(M + 1, j + 1, M + 1)];
            }
        }

        // Reosc at next if current segment passes through perigee
        orb_end = 0.0;
        // FIXME: Modify forward and aft next initial conditions when the perigee reosculates. If the bootstrap orbit reosculates then the time vector is modified as is the next segments initial conditions.
        reosc_perigee(tf, *this);
        if (orb_end != 0.0)
        {
            segment_times[seg_cnt + 1] = orb_end;
        }
        else
        {
            segment_times[seg_cnt + 1] = tf;
        }
        // debug print segment start and end time
        if (g_DEBUG_PICARD)
        {
            cout << "Segment end time: " << segment_times[seg_cnt + 1] << std::endl;
        }
        // Total segments counter
        seg_cnt = seg_cnt + 1;
        k = k + 1;
        // Check if the end of the integration time has been reached
        if (fabs(tf - _Integrator_tf) / tf < 1e-12)
        {
            loop = 1;
        }
    }
}
/**
 * @brief This function performs picard iteration on a segment for all of the orbit objects in the orbit list.
 * Each iteration, all orbits are computed before moving on to the next iteration. This is done to ensure that
 * the gravity calculations from the for and aft orbit can be used with the bootstrap orbit.
 */
void BootstrapOrbit::Bootstrap_Picard_Iteration(EphemerisManager &ephem)
{
    int MaxIt = 300; // Maximum number of iterations
    // Propagator values that are consistent between each orbit object
    int N = this->N;
    int M = this->M;
    int hot = this->prep_HS;
    double tol = this->tol;
    int deg = this->deg;
    vector<double> times = this->times_seg;

    // Quadrature matrices
    vector<double> P1 = this->P1;
    vector<double> P2 = this->P2;
    vector<double> T1 = this->T1;
    vector<double> T2 = this->T2;
    vector<double> A = this->A;

    // Initialized node variables that may be reused at each step
    bool suborbital = false;
    double alt = 0.0;
    double xI[3] = {0.0};
    double vI[3] = {0.0};
    double xPrimaryFixed[3] = {0.0};
    double vPrimaryFixed[3] = {0.0};
    double aPrimaryFixed[3] = {0.0};
    double aI[3] = {0.0};
    double del_X[3] = {0.0};
    double del_aECI[3] = {0.0};
    double drag_aECEF[3] = {0.0};
    double SRP_aI[3] = {0.0};
    double third_body_aI[3] = {0.0};

    // initialized segment variables that may be reused for each iteration (do not carry necessary info between iterations)
    std::vector<double> beta(N * 3, 0.0);
    std::vector<double> alpha((N + 1) * 3, 0.0);
    std::vector<double> gamma(N * 3, 0.0);
    std::vector<double> kappa((N + 1) * 3, 0.0);
    std::vector<double> Xorig;
    std::vector<double> Vorig;
    std::vector<double> Xnew;
    std::vector<double> Vnew;
    std::vector<double> xECEFp((M + 1) * 3, 0.0);
    std::vector<double> xECIp((M + 1) * 3, 0.0);
    std::vector<double> del_a((M + 1) * 3, 0.0);
    double w2 = (times[M] - times[0]) / 2;
    vector<double> del_G_normal;
    if (g_DEBUG_BOOTSRAP)
    {
        del_G_normal.resize((Nmax + 1) * 3, 0.0);
    }
    // Initialize orbit specific variables
    for (Orbit *orbitRef : orbitslist)
    {
        Orbit &orbit = *orbitRef;
        orbit.G.resize((M + 1) * 3, 0.0);
        orbit.del_G.resize((Nmax + 1) * 3);
        fill(orbit.del_G.begin(), orbit.del_G.end(), 0.0);
        orbit.itr = 0;
        orbit.ITRs = IterCounters(); // reset all counters to 0
        if (hot == 1)
        {
            orbit.err = 1e-2;
        }
        else
        {
            orbit.err = 10;
        }
        orbit.converged = false;
    }
    this -> Exit_Bootstrap = false;
    bool notConverged = true;
    // Start the iteration loop
    while (notConverged)
    {
        // Each iteration compute for each orbit object
        for (Orbit *orbitRef : orbitslist)
        {
            Orbit &currOrbit = *orbitRef;
            if (currOrbit.converged)
            {
                continue;
            }
            vector<double> &X = currOrbit.X_seg;
            vector<double> &V = currOrbit.V_seg;
            vector<double> &G = currOrbit.G;
            double *del_G = &currOrbit.del_G[0];
            double *Xint = &currOrbit.r0_seg[0];
            double *Vint = &currOrbit.v0_seg[0];
            vector<double> &Beta = currOrbit.Beta_seg;
            vector<double> &Alpha = currOrbit.Alpha_seg;
            IterCounters &ITRs = currOrbit.ITRs;
            int &itr = currOrbit.itr;
            double &err = currOrbit.err;
            double* Feval;
            if(&currOrbit == this && Bootstrap_On && !this->Exit_Bootstrap)
            {
                //point to the counter for the bootstrap orbit
                Feval = &currOrbit.Feval.Bootstrap[0];
            }
            else
            {
                // point to the counter for the picard iteration counter
                Feval = &currOrbit.Feval.PicardIteration[0];
            }
            // cout << currOrbit.name << " Iteration  " << itr << " error:" << err << endl;
            for (int i = 1; i <= M + 1; i++) // Get forces at each node on the segment in inertial frame
            {

                for (int j = 1; j <= 3; j++)
                {
                    xI[j - 1] = X[ID2(i, j, M + 1)];
                    vI[j - 1] = V[ID2(i, j, M + 1)];
                }

                // Exit loop early if xI or vI is NaN
                if (isnan(xI[0]) || isnan(xI[1]) || isnan(xI[2]) || isnan(vI[0]) || isnan(vI[1]) || isnan(vI[2]))
                {
                    // print error message
                    cout << "Error: NaN in Picard Iteration" << endl;
                    // end picard iteration
                    return;
                }
                // Convert from ECI to ECEF
                // InertialToBodyFixed(xI,vI,xPrimaryFixed,vPrimaryFixed,times[i-1],orbit);
                eci2ecef(et(times[i - 1]), xI, vI, xPrimaryFixed, vPrimaryFixed);
                // Compute Variable Fidelity Gravity for the for and aft satellites only

                if (orbitRef == this && Bootstrap_On && !this->Exit_Bootstrap)
                {
                    // Compute low fidelity
                    lunar_Grav_Approx_Function(times[i - 1], xPrimaryFixed, aPrimaryFixed, Feval, lowDeg);
                    // DEBUG check what the normal calculation would have been

                    // Compute approximation for the bootstrap orbit
                    for (int j = 0; j <= 2; j++)
                    {
                        del_G[ID2(i, j + 1, Nmax + 1)] = 0.5 * aftOrbit.del_G[ID2(i, j + 1, Nmax + 1)] + 0.5 * forOrbit.del_G[ID2(i, j + 1, Nmax + 1)];
                        aPrimaryFixed[j] += del_G[ID2(i, j + 1, Nmax + 1)];
                    }
                    // if (g_DEBUG_BOOTSRAP)
                    // {
                    //     double aPrimaryFixed_normal[3] = {0.0};
                    //     lunar_perturbed_gravity(times[i - 1], xPrimaryFixed, err, i, M, deg, hot, aPrimaryFixed_normal, tol, &itr, Feval, ITRs, &del_G_normal[0], lowDeg);

                    //     double delg_G_mag = 0.0;
                    //     double del_G_err_norm[3] = {0.0};
                    //     for (int j = 0; j < 3; j++)
                    //     {
                    //         delg_G_mag += del_G[j] * del_G[j];
                    //     }
                    //     delg_G_mag = sqrt(delg_G_mag);
                    //     for (int j = 0; j < 3; j++)
                    //     {
                    //         double del_G_err = (del_G[j] - del_G_normal[j]);
                    //         if (delg_G_mag == 0.0 && del_G_err != 0.0)
                    //         {
                    //             del_G_err_norm[j] = INFINITY;
                    //         }
                    //         else if (delg_G_mag == 0.0 && del_G_err == 0.0)
                    //         {
                    //             del_G_err_norm[j] = 0.0;
                    //         }

                    //         else
                    //         {
                    //             del_G_err_norm[j] = del_G_err / delg_G_mag;
                    //         }
                    //     }
                    //     cout << "del_G_err: " << del_G_err_norm[0] << "\t" << del_G_err_norm[1] << "\t" << del_G_err_norm[2] << endl;
                    // }
                }
                else
                {
                    // run normal APC procedures
                    lunar_perturbed_gravity_error(times[i - 1], xPrimaryFixed, err, i, M, deg, hot, aPrimaryFixed, tol, &itr, Feval, ITRs, del_G, lowDeg);
                }

                // Calculate acceleration from drag
                Perturbed_Drag(xPrimaryFixed, vPrimaryFixed, currOrbit, drag_aECEF);

                // sum pertubed gravity and drag accelerations
                for (int k = 0; k < 3; k++)
                {
                    aPrimaryFixed[k] = aPrimaryFixed[k] + drag_aECEF[k];
                }

                // Convert acceleration vector from ECEF to ECI
                // BodyFixedAccelerationToInertial(aPrimaryFixed,aI,times[i-1],orbit);
                ecef2eci(et(times[i - 1]), aPrimaryFixed, aI);
                // calculate SRP and Third Body
                Perturbed_SRP(times[i - 1], xI, currOrbit, ephem, SRP_aI);
                Perturbed_three_body_moon(times[i - 1], xI, currOrbit, ephem, third_body_aI);
                // Add perturbations to acceleration
                for (int k = 0; k < 3; k++)
                {
                    aI[k] = aI[k] + SRP_aI[k] + third_body_aI[k];
                }
                for (int j = 1; j <= 3; j++)
                {
                    G[ID2(i, j, M + 1)] = aI[j - 1];
                    xECIp[ID2(i, j, M + 1)] = xI[j - 1];
                }
            } // Forces found for each node

            // Perform quadrature for velocity then position in inertial frame
            // Velocity
            std::vector<double> tmp1;
            std::vector<double> tmp2;
            tmp1 = matmul(A, G, N - 1, M + 1, 3, N - 1, M + 1);
            tmp2 = matmul(P1, tmp1, N, N - 1, 3, N, N - 1);
            for (int i = 1; i <= N; i++)
            {
                for (int j = 1; j <= 3; j++)
                {
                    beta[ID2(i, j, N)] = w2 * tmp2[ID2(i, j, N)];
                    if (i == 1)
                    {
                        beta[ID2(i, j, N)] = beta[ID2(i, j, N)] + Vint[j - 1];
                    }
                }
            }
            Vorig = matmul(T1, beta, M + 1, N, 3, M + 1, N);

            // Position
            std::vector<double> tmp3;
            tmp3 = matmul(P2, beta, N + 1, N, 3, N + 1, N);
            for (int i = 1; i <= N + 1; i++)
            {
                for (int j = 1; j <= 3; j++)
                {
                    alpha[ID2(i, j, N + 1)] = w2 * tmp3[ID2(i, j, N + 1)];
                    if (i == 1)
                    {
                        alpha[ID2(i, j, N + 1)] = alpha[ID2(i, j, N + 1)] + Xint[j - 1];
                    }
                }
            }
            Xorig = matmul(T2, alpha, M + 1, N + 1, 3, M + 1, N + 1);

            for (int i = 1; i <= M + 1; i++)
            {
                for (int j = 1; j <= 3; j++)
                {
                    xI[j - 1] = Xorig[ID2(i, j, M + 1)];
                    vI[j - 1] = Vorig[ID2(i, j, M + 1)];
                }
                // Linear Error Correction Position
                for (int j = 1; j <= 3; j++)
                {
                    del_X[j - 1] = xI[j - 1] - xECIp[ID2(i, j, M + 1)];
                }
                // Convert from inertial to body fixed frame
                // InertialToBodyFixed(xI,vI,xPrimaryFixed,vPrimaryFixed,times[i-1],orbit);
                eci2ecef(et(times[i - 1]), xI, vI, xPrimaryFixed, vPrimaryFixed);
                // Linear Error Correction Acceleration
                double del_aECEF[3] = {0.0};
                picard_error_feedback_GRGM1200b(xPrimaryFixed, del_X, del_aECEF);
                // Convert from ECEF to ECI
                ecef2eci(et(times[i - 1]), del_aECEF, del_aECI);
                // BodyFixedAccelerationToInertial(del_aECEF,del_aECI,times[i-1],orbit);

                for (int j = 1; j <= 3; j++)
                {
                    del_a[ID2(i, j, M + 1)] = del_aECI[j - 1];
                }
            }

            // Linear Error Correction Velocity Coefficients
            std::vector<double> tmp4;
            std::vector<double> tmp5;
            tmp4 = matmul(A, del_a, N - 1, M + 1, 3, N - 1, M + 1);
            tmp5 = matmul(P1, tmp4, N, N - 1, 3, N, N - 1);
            for (int i = 1; i <= N; i++)
            {
                for (int j = 1; j <= 3; j++)
                {
                    gamma[ID2(i, j, N)] = w2 * tmp5[ID2(i, j, N)];
                }
            }

            // Corrected Velocity
            for (int i = 1; i <= N; i++)
            {
                for (int j = 1; j <= 3; j++)
                {
                    Beta[ID2(i, j, N)] = beta[ID2(i, j, N)];
                    if (err < 1e-13)
                    {
                        Beta[ID2(i, j, N)] = beta[ID2(i, j, N)] + gamma[ID2(i, j, N)];
                    }
                }
            }
            Vnew = matmul(T1, Beta, M + 1, N, 3, M + 1, N);

            // Corrected Position
            std::vector<double> tmp6;
            tmp6 = matmul(P2, gamma, N + 1, N, 3, N + 1, N);
            for (int i = 1; i <= N + 1; i++)
            {
                for (int j = 1; j <= 3; j++)
                {
                    kappa[ID2(i, j, N + 1)] = w2 * tmp6[ID2(i, j, N + 1)];
                    Alpha[ID2(i, j, N + 1)] = alpha[ID2(i, j, N + 1)];
                    if (err < 1e-13)
                    {
                        Alpha[ID2(i, j, N + 1)] = alpha[ID2(i, j, N + 1)] + kappa[ID2(i, j, N + 1)];
                    }
                }
            }
            Xnew = matmul(T2, Alpha, M + 1, N + 1, 3, M + 1, N + 1);

            // Non-dimensional Error
            double tmp = 0.0;
            double curr_err = 0.0;
            for (int i = 1; i <= M + 1; i++)
            {
                for (int j = 1; j <= 6; j++)
                {
                    if (j <= 3)
                    {
                        tmp = fabs(Xnew[ID2(i, j, M + 1)] - X[ID2(i, j, M + 1)]) / DU;
                    }
                    if (j > 3)
                    {
                        tmp = fabs(Vnew[ID2(i, j - 3, M + 1)] - V[ID2(i, j - 3, M + 1)]) / DU * TU;
                    }
                    if (tmp > curr_err)
                    {
                        curr_err = tmp;
                    }
                }
            }
            err = curr_err;

            // check for convergence
            if (err < tol)
            {
                if (ITRs.DID_FULL)
                // if(true)
                {
                    currOrbit.converged = true;
                    if (g_DEBUG_PICARD)
                    {
                        cout << currOrbit.name << " converged in " << itr << " iterations." << endl;
                    }
                }
            }
            // Update
            X = Xnew;
            V = Vnew;
            itr++;
            // Iteration Counter


            // End of orbit object for loop
        }
        if (aftOrbit.converged && forOrbit.converged && (this->converged || !Bootstrap_To_Convergence))
            {
                notConverged = false;
            }
        if(aftOrbit.converged && forOrbit.converged)
        {
            this->Exit_Bootstrap = true;
        }
        // End of iteration while loop
    }
    if (g_DEBUG_PICARD)
    {
        for (Orbit *orbitRef : orbitslist)
        {
            Orbit &orbit = *orbitRef;
            segment seg;
            seg.itrs = orbit.itr;
            seg.err = orbit.err;
            seg.converged = orbit.converged;
            orbit.DebugData.segments.push_back(seg);
        }
    }
    return;
}