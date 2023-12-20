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

#include "SpiceUsr.h"
#include "Ephemeris.hpp"
#include "APC.h"
#include "flags.h"
#include "FandG.h"

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

// Bootstrap Orbit

BootstrapOrbit::BootstrapOrbit(const Orbit &orbit, double followTime) : forOrbit(_primary,_IOFrame,_epoch), aftOrbit(_primary, _IOFrame, _epoch), Orbit(orbit)
{
    double z0[6] = {_In_r0[0], _In_r0[1], _In_r0[2], _In_v0[0], _In_v0[1], _In_v0[2]};
    double zfor[6] = {0.0};
    double zaft[6] = {0.0};
    // Generate the state ahead and behind of the target orbit
    FandG(z0, zfor, followTime, _mu);
    FandG(z0, zaft, -followTime, _mu);
    // assign values to the forward orbit
    forOrbit.SetPosition0({zfor[0], zfor[1], zfor[2]});
    forOrbit.SetVelocity0({zfor[3], zfor[4], zfor[5]});
    forOrbit.SetIntegrationTime(_Ephemeris_t0, _Ephemeris_tf);
    forOrbit.SetProperties(satproperties._Area, satproperties._Reflectance, satproperties._Mass, satproperties._Cd);
    forOrbit.SetComputeHamiltonian(Compute_Hamiltonian);
    forOrbit.SetComputeSRP(Compute_SRP);
    forOrbit.SetComputeThirdBody(Compute_Third_Body);
    forOrbit.SetComputeDrag(Compute_Drag);
    //assign values to the aft orbit
    aftOrbit.SetPosition0({zaft[0], zaft[1], zaft[2]});
    aftOrbit.SetVelocity0({zaft[3], zaft[4], zaft[5]});
    aftOrbit.SetIntegrationTime(_Ephemeris_t0, _Ephemeris_tf);
    aftOrbit.SetProperties(satproperties._Area, satproperties._Reflectance, satproperties._Mass, satproperties._Cd);
    aftOrbit.SetComputeHamiltonian(Compute_Hamiltonian);
    aftOrbit.SetComputeSRP(Compute_SRP);
    aftOrbit.SetComputeThirdBody(Compute_Third_Body);
    aftOrbit.SetComputeDrag(Compute_Drag);
}

BootstrapOrbit::BootstrapOrbit()
{
}

BootstrapOrbit::BootstrapPropagate()
{
}
