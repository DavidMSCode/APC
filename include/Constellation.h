/**
 * @file Constellation.h
 * @author David Stanley (davidms4@illinois.edu)
 * @brief 
 * @version 0.1
 * @date 2023-04-04
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef __CONSTELLATION__
#define __CONSTELLATION__
#include <vector>
#include "Orbit.h"
using namespace std;

/**
 * Takes a list of orbit class objects and organizes them for simultaneous integration using APC. Constellation is a child class of the orbit class and adopts all its methods.
*/
class Constellation: public Orbit
{
    private:
        vector<Orbit> _Orbit_List;  //List of all orbit objects in constellation
                
    public:
        Constellation();
        Constellation(vector<Orbit> Orbit_List);
};

#endif