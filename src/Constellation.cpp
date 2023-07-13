/**
 * @file Constellation.cpp
 * @author David Stanley (davidms4@illinois.edu)
 * @brief 
 * @version 0.1
 * @date 2023-04-04
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include <vector>
#include "Constellation.h"
#include "Orbit.h"

using namespace std;


Constellation::Constellation()
{
    
}

Constellation::Constellation(vector<Orbit> Orbit_list)
{
    _Orbit_List = Orbit_list;
    //TODO: Iterate through orbit list and ensure each orbit is using the same values as the overall constellation
}

void Constellation::SharedPropagate(){

}
