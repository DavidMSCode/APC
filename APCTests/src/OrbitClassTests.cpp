/**
 * @file OrbitClassTests.cpp
 * @author David Stanley (davidms4@illinois.edu)
 * @brief Runs tests on various properties and methods of the Orbit class. 
 * @version 0.1
 * @date 2023-03-10
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include <iostream>
#include "Orbit.h"
#include "OrbitClassTests.h"
using namespace std;

int main() {
    std::cout << "Orbit Class Tests started!" << std::endl;
    int primaryTestFailures = testValidPrimaries();           //Primary body name tests
    return primaryTestFailures;
};

int testValidPrimaries(){
    int failures = 0;
    //Valid Inputs for primary body
    Orbit orbit(10,.8,100,2.2,false,false,false,false,1,"Earth");
    if(!orbit.ValidPrimary()){
        failures++;
        cout << "\"Earth\" should be a valid primary body name, but Orbit::ValidPrimary returned False." << endl;
    }
    Orbit orbit2(10,.8,100,2.2,false,false,false,false,1,"Moon");
    if(!orbit2.ValidPrimary()){
        failures++;
        cout << "\"Moon\" should be a valid primary body name, but Orbit::ValidPrimary returned False." << endl;
    }
    Orbit orbit3(10,.8,100,2.2,false,false,false,false,1,"earth");
    if(!orbit3.ValidPrimary()){
        failures++;
        cout << "\"earth\" should be a valid primary body name, but Orbit::ValidPrimary returned False." << endl;
    }
    Orbit orbit4(10,.8,100,2.2,false,false,false,false,1,"moon");
    if(!orbit4.ValidPrimary()){
        failures++;
        cout << "\"moon\" should be a valid primary body name, but Orbit::ValidPrimary returned False." << endl;
    }
    Orbit orbit5(10,.8,100,2.2,false,false,false,false,1);
    if(!orbit5.ValidPrimary()){
        failures++;
        cout << "Primary body name should default to \"Earth\" and be a valid name, but Orbit::ValidPrimary returned False." << endl;
    }
    return failures;
};