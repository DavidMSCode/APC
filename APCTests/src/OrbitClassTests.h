

#ifndef _OrbitTests_
#define _OrbitTests_

/**
 * @brief Checks whether or not the Orbit class method PrimaryValid() returns proper truth values for different scenarios.
 * 
 * @return int, returns the number of failed checks
 */
int OrbitClassTests(int argc, char** argv);
int testValidPrimaries();

int testInitOrbit();

int testSetState();

int testSetProperties();

int testSetIntegrationTime();

#endif