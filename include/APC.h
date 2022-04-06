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
#include <Orbit.h>
#include <Ephemeris.hpp>

/**
 * @brief 
 * 
 * @param t0 
 * @param tf 
 * @return EmphemerisManager 
 */
class EphemerisManager cacheEphemeris(double t0, double tf);

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
 */
struct SatState
{
    std::vector<double> r;
    std::vector<double> v;
};

/**
 * @brief 
 * 
 * @param r 
 * @param v 
 * @param pos_error 
 * @param vel_error 
 * @return std::vector<SatState> 
 */
std::vector<SatState> GenSigma13(std::vector<double> r, std::vector<double> v, double pos_error, double vel_error);

/**
 * @brief 
 * 
 * @param r 
 * @param v 
 * @param pos_error 
 * @param vel_error 
 * @return std::vector<SatState> 
 */
std::vector<SatState> GenSigma3(std::vector<double> r, std::vector<double> v, double pos_error, double vel_error);

/**
 * @brief 
 * 
 * @param sigma13 
 */
void printStateList(std::vector<SatState> sigma13);

/**
 * @brief 
 * 
 * @param state 
 * @param i 
 */
void printState(SatState state, int i);


/**
 * @brief 
 * 
 * @param r 
 * @param v 
 * @param t0 
 * @param tf 
 * @param area 
 * @param reflectance 
 * @param mass 
 * @param drag_C 
 * @param compute_drag 
 * @param compute_SRP 
 * @param compute_third_body 
 * @return std::vector<class Orbit> 
 */
std::vector<Orbit> ParallelPropagate(std::vector<SatState> StateList, double t0, double tf, double area, double reflectance, double mass, double drag_C, bool compute_drag, bool compute_SRP, bool compute_third_body);

/**
 * @brief 
 * 
 * @param r 
 * @param v 
 * @param t0 
 * @param tf 
 * @param area 
 * @param reflectance 
 * @param mass 
 * @param drag_C 
 * @param compute_drag 
 * @param compute_SRP 
 * @param compute_third_body 
 * @return std::vector<class Orbit> 
 */
class Orbit SinglePropagate(std::vector<double> r, std::vector<double> v, double t0, double tf, double area, double reflectance, double mass, double drag_C, bool compute_drag, bool compute_SRP, bool compute_third_body);

/**
 * @brief 
 * 
 * @param ephem 
 * @param t0 
 * @param tf 
 */
void MPGetTest(EphemerisManager ephem, double t0, double tf);
#endif