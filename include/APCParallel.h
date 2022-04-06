/**
 * @file APCParallel.h
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2022-03-23
 * 
 * @copyright Copyright (c) 2022
 * 
 */


#ifndef __PARALLEL__
#define __PARALLEL__

#include <vector>
#include <Orbit.h>

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
void ParallelPropagateTest(std::vector<SatState> StateList, double t0, double tf, double area, double reflectance, double mass, double drag_C, bool compute_drag, bool compute_SRP, bool compute_third_body);


#endif