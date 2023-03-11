

#ifndef __EPHEMTESTS__
#define __EPHEMTESTS__
#include "Ephemeris.hpp"

/**
 * @brief Tests the thread safety of the ephemeris manager by reading data from the EphemerisManager and writing to stdout in parallel
 * 
 * @param ephem EphemerisManager object containing ephermeris data for solar system objects over a given interval
 * @param t0  Beginning of EphemerisManager time interval (s)
 * @param tf End of EphemerisManager time interval (s)
 */
int MPGetTest(EphemerisManager &ephem, double t0, double tf);

#endif