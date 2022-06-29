/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com), Alex Pascarella (alexp3@illinois.edu)
*  DATE WRITTEN:     Dec 2020
*  LAST MODIFIED:    Feb 2021
*  AFFILIATION:      University of Illinois at Urbana-Champaign (UIUC)
*  DESCRIPTION:      Ephemeris computation using Chebyshev polynomials
*/

#ifndef _EPHEM_H_
#define _EPHEM_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <cassert>
#include <list>
#include <map>

#include "const.hpp"
#include "c_functions2.hpp"
#include "c_functions.h"
#include "Ephemeris.hpp"
#include "SpiceUsr.h"

using namespace std;

// Class for storing the ephemeris of a single celestial body using Chebyshev polynomials
class ChebyshevEphemeris
{
    private:
        int N;
        double  initialEpoch, finalEpoch, segmentLength;
        vector< double > w1, w2;
        vector< vector < double > > coefficients;
        string body, center, frame;

    public:
        // Default constructor
        ChebyshevEphemeris();

        // Overloaded constructor
        // INPUTS:
        // spkFile
        // lskFile
        // initialEpoch: initial ephemeris epoch, in seconds after J2000
        // finalEpoch: final ephemeris epoch, in seconds after J2000
        // body: celestial body for which the ephemeris is computed
        // center: origin of the frame used for the ephemeris computation
        // frame: frame used for the ephemeris computation
        ChebyshevEphemeris(string spkFile,
                           string lskFile,
                           double initialEpoch,
                           double finalEpoch,
                           string body,
                           string center,
                           string frame);

        // Function to compute the ephemeris of the body at a given epoch
        void getState(double* state, double epoch);

        // Getters
        string getBody();
        string getCenter();
        string getFrame();
        double getInitialEpoch();
        double getFinalEpoch();

        // Static methods
        static void fitChebyshevPolynomials(double s, int N, int M, std::vector<double> &T, std::vector<double> &A);
        static void getChebyshevPolynomials(double s, int N, int M, int arg, std::vector<double> &T);
};

// Class for managing the ephemeris of multiple celestial bodies
class EphemerisManager
{
    private:
        map<string, ChebyshevEphemeris> ephemeris;

    public:
        // Overloaded constructor
        // INPUTS:
        // spkFile
        // lskFile
        // initialEpoch: initial ephemeris epoch, in seconds after J2000
        // finalEpoch: final ephemeris epoch, in seconds after J2000
        // bodies: list of celestial body for which the ephemerides are computed
        // center: origin of the frame used for the ephemeris computation
        // frame: frame used for the ephemeris computation
        EphemerisManager();
        EphemerisManager(string spkFile,
                         string lskFile,
                         double initialEpoch,
                         double finalEpoch,
                         list<string> bodies,
                         string center,
                         string frame);

        // Function to retrieve the ephemeris of a body at a given epoch
        std::vector<double> getState(string body, double epoch);
};

#endif