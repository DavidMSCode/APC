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
    int initTestFailures = testInitOrbit();
    int initStateTestFailures = testSetState();
    return primaryTestFailures+initTestFailures+initStateTestFailures;
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
}
int testInitOrbit()
{
    int failures;
    //Initialize orbit around the moon, in moon primary axis frame with J2000 epoch.
    string primary = "MOON";
    string IOFrame = "MOON_PA";
    string epoch = "J2000";
    Orbit orbit(primary,IOFrame,epoch);
    //check if primary was set properly
    if (orbit.InSet(orbit.GetPrimaryBody(),{primary})+
        orbit.InSet(orbit.GetIOFrame(),{IOFrame})+
        orbit.InSet(orbit.GetEpoch(),{epoch})!=3){
        failures+=1;
        cout << "Orbit initialization with primary axis frame centered on the moon did not initialize properly"<<endl;
    }
    return failures;
}

int testSetState()
{
    int failures=0;
    vector<double> r0 = {1838.,         0.,         0.};    // Initial Position (km)
    vector<double> v0 = {0.   , 0.74844211, 1.49688755};
    vector<double> s0;
    s0.reserve( r0.size() + v0.size() ); // preallocate memory
    s0.insert( s0.end(), r0.begin(), r0.end() );
    s0.insert( s0.end(), v0.begin(), v0.end() );
    string primary = "MOON";
    string IOFrame = "MOON_PA";
    string epoch = "J2000";
    Orbit orbit(primary,IOFrame,epoch);
    orbit.SetState0(s0);
    if(orbit.getPosition0()!=r0){
        failures+=1;
    }
    if(orbit.getVelocity0()!=v0){
        failures+=1;
    }
    if(orbit.getState0()!=s0){
        failures+=1;
    }

    if (failures>0){
        cout<<"Orbit state initialization test failed."<<endl;
    }
    return failures;
}
