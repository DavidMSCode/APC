/**
 * @file Constellation.h
 * @author David Stanley (davidms4@illinois.edu)
 * @ Modified by: Your name
 * @ Modified time: 2023-07-06 12:49:44
 * @date 2023-03-10
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

                
    public:

        vector<Orbit> _Orbit_List;  //List of all orbit objects in constellation

        //Constructors
        Constellation();
        Constellation(vector<Orbit> Orbit_List);
        //Methods
        void SharedPropagate();
        
};

#endif