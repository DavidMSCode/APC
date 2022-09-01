/*
*  AUTHORS:          David Stanley (davidms4@illinois.edu)
*  DATE WRITTEN:     Feb 2022
*  LAST MODIFIED:    Feb 2022
*  AFFILIATION:      Department of Aerospace Engineering, University of Illinois Urbana-Champaign
*  DESCRIPTION:      Methods that are acessible from python and binding code
*  REFERENCE:        Woollands, R., and Junkins, J., "Nonlinear Differential Equation Solvers
*                    via Adaptive Picard-Chebyshev Iteration: Applications in Astrodynamics", JGCD, 2016.
*/

#ifndef __APCMAIN__
#define __APCMAIN__

#include <vector>

#include "Orbit.h"
#include "Ephemeris.hpp"

/**
 * @brief 
 * 
 * @param t0 
 * @param tf 
 * @return EmphemerisManager 
 */
EphemerisManager cacheEphemeris(double t0, double tf);

/**
 * @brief Returns propagated orbit solution given initial conditions.
 * 
 * The inputs of this function are given in an Earth centered inertial frame and returned in an Earth
 * centered inertial frame. The output vectors are stored in a single vector such that each index corresponds
 * to a single output vector.
 * t=output[0],
 * position (x,y,z) = output[1:3],
 * velocity (u,v,w) = output[4:6],
 * Hamiltonian = output[7].
 * 
 * @param r Initial position of satellite (km)
 * @param v Initial velocity of satellite (km/s)
 * @param t0 Epoch start time of satellite (s)
 * @param tf Epoch end time of satellite (s)
 * @param area "Cannonball" area of satellite (m^2)
 * @param reflectance Coefficient of relectance of satellite
 * @param mass Mass of satellite (kg)
 * @param drag_C Drag coefficient of satellite
 * @param compute_drag Boolean to toggle drag computation on or off
 * @param compute_SRP Boolean to toggle solar radiation pressure computation on or off
 * @param compute_third_body Boolean to toggle third body perturbation computation on or off
 * @return std::vector<std::vector<double> > The solutions stored in a vector of vectors {t,x,y,z,u,v,w,H}
 */
std::vector<std::vector<double> > PropagateICs(std::vector<double> r, std::vector<double> v, double t0, double tf, Orbit orbit, EphemerisManager ephem);

/**
 * @brief Returns propagated orbit solution as an Orbit object given initial conditions.
 * 
 * The inputs of this function are given in an Earth centered inertial frame and returned in an Earth
 * centered inertial frame. The solution is stored in Orbit class object.
 * @see Orbit
 * 
 * @param r Initial position of satellite (km)
 * @param v Initial velocity of satellite (km/s)
 * @param t0 Epoch start time of satellite (s)
 * @param tf Epoch end time of satellite (s)
 * @param area "Cannonball" area of satellite (m^2)
 * @param reflectance Coefficient of relectance of satellite
 * @param mass Mass of satellite (kg)
 * @param drag_C Drag coefficient of satellite
 * @param compute_drag Boolean to toggle drag computation on or off
 * @param compute_SRP Boolean to toggle solar radiation pressure computation on or off
 * @param compute_third_body Boolean to toggle third body perturbation computation on or off
 * @return Orbit Orbit object that stores the solution and other relevant parameters
 */
class Orbit PropagateOrbit(std::vector<double> r, std::vector<double> v, double t0, double tf, Orbit orbit, EphemerisManager ephem);

/**
 * @brief Structure containing the position and velocity of some orbit state.
 * 
 * @param r Satellite position (km)
 * @param v Satellite velocity (km/s)
 */
struct SatState
{
    std::vector<double> r;
    std::vector<double> v;
};

/**
 * @brief Takes a satellite initial condition and returns 13 derived satellite initial conditions that were perturbed in velocty or position
 * 
 * @param r Initial satellite position (km)
 * @param v Initial satellite velocity (km/s)
 * @param pos_error Position perturbation (km)
 * @param vel_error Velocity perturbation (km/s)
 * @return std::vector<SatState> The output vector of 13 satellite initial conditions
 */
std::vector<SatState> GenSigma13(std::vector<double> r, std::vector<double> v, double pos_error, double vel_error);

/**
 * @brief Function writes the positions and velocities and id of a list of satellite states to stdout
 * 
 * @param statelist A vector contatining satellite states
 */
void printStateList(std::vector<SatState> statelist);

/**
 * @brief Writes a satellite state to stdout
 * 
 * @param state A struct of a satellites position and velocity
 * @param id An integer id 
 */
void printState(SatState state, int id);


/**
 * @brief Propagates multiple satellite orbits in parallel
 * 
 * @param StateList A vector of satellite initial conditions in J2000
 * @param t0 Epoch start time of satellite (s)
 * @param tf Epoch end time of satellite (s)
 * @param area "Cannonball" area of satellite (m^2)
 * @param reflectance Coefficient of relectance of satellite
 * @param mass Mass of satellite (kg)
 * @param drag_C Drag coefficient of satellite
 * @param compute_drag Boolean to toggle drag computation on or off
 * @param compute_SRP Boolean to toggle solar radiation pressure computation on or off
 * @param compute_third_body Boolean to toggle third body perturbation computation on or off
 * @return std::vector<class Orbit> Returns a vector of orbit objects that each describe the satellite properties and the propagation solution
 */
std::vector<Orbit> ParallelPropagate(std::vector<SatState> StateList, double t0, double tf, double area, double reflectance, double mass, double drag_C, bool compute_drag, bool compute_SRP, bool compute_third_body);

/**
 * @brief Propagates a single satellite orbit
 * 
* @param r Initial position of satellite (km)
 * @param v Initial velocity of satellite (km/s)
 * @param t0 Epoch start time of satellite (s)
 * @param tf Epoch end time of satellite (s)
 * @param area "Cannonball" area of satellite (m^2)
 * @param reflectance Coefficient of relectance of satellite
 * @param mass Mass of satellite (kg)
 * @param drag_C Drag coefficient of satellite
 * @param compute_drag Boolean to toggle drag computation on or off
 * @param compute_SRP Boolean to toggle solar radiation pressure computation on or off
 * @param compute_third_body Boolean to toggle third body perturbation computation on or off
 * @return std::vector<class Orbit> Returns an orbit object that describes the satellite properties and the propagation solution
 */
class Orbit SinglePropagate(std::vector<double> r, std::vector<double> v, double t0, double tf, double area, double reflectance, double mass, double drag_C, bool compute_drag, bool compute_SRP, bool compute_third_body);

/**
 * @brief Propagates a single satellite orbit using a user defined time vector
 * 
* @param r Initial position of satellite (km)
 * @param v Initial velocity of satellite (km/s)
 * @param t0 Epoch start time of satellite (s)
 * @param tf Epoch end time of satellite (s)
 * @param area "Cannonball" area of satellite (m^2)
 * @param reflectance Coefficient of relectance of satellite
 * @param mass Mass of satellite (kg)
 * @param drag_C Drag coefficient of satellite
 * @param compute_drag Boolean to toggle drag computation on or off
 * @param compute_SRP Boolean to toggle solar radiation pressure computation on or off
 * @param compute_third_body Boolean to toggle third body perturbation computation on or off
 * @return std::vector<class Orbit> Returns an orbit object that describes the satellite properties and the propagation solution
 */
class Orbit SinglePropagate(std::vector<double> r, std::vector<double> v, std::vector<double> time_vec, double area, double reflectance, double mass, double drag_C, bool compute_drag, bool compute_SRP, bool compute_third_body);

/**
 * @brief Tests the thread safety of the ephemeris manager by reading data from the EphemerisManager and writing to stdout in parallel
 * 
 * @param ephem EphemerisManager object containing ephermeris data for solar system objects over a given interval
 * @param t0  Beginning of EphemerisManager time interval (s)
 * @param tf End of EphemerisManager time interval (s)
 */
void MPGetTest(EphemerisManager ephem, double t0, double tf);

/**
 * @brief Calculates 1000 orbits to benchmark parallel processing performance
 * 
 * @return std::pair<int,double> Returns a pair of values indicating the number of threads used and the time to finish the benchmark respectively.
 */
std::pair<int,double>  Benchmark1000(int max_threads);

/**
 * @brief Queries MATRICES_LAODED flag to find out if EGM2008 matrices are in memory
 * 
 */
void printMatrixState();

/**
 * @brief Returnes true if EGM2008 matrices have been loaded.
 * 
 */
bool MatricesLoaded();
#endif

