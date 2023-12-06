
#include "const.h"
#include "APC.h"

#include <fstream>
#include <omp.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <map>
#include <random>
#include <chrono>
#include <sstream>

// Function to calculate E using Newton's method
double calculateE(double M, double e, double tol = 1e-16)
{
    double E = M;
    double g = 1;
    int iters = 0;
    int max_iter = 100000;

    while (std::abs(g) > tol && iters < max_iter)
    {
        g = E - e * std::sin(E) - M;
        double dgde = 1 - e * std::cos(E);
        E = E - g / dgde;
        iters++;
    }

    return E;
}

// Function to convert orbital elements to position and velocity vectors
std::vector<std::vector<double>> elms2rv(double a, double e, double inc, double Om, double w, double M, double mu)
{
    double tol = 1e-10;

    // Calculate E
    double E = calculateE(M, e, tol);

    // Calculate true anomaly
    double f = 2 * std::atan(std::sqrt((1 + e) / (1 - e)) * std::tan(E / 2));

    // Calculate radius and speed
    double h = std::sqrt(mu * a * (1 - e * e));
    double r = a * (1 - e * std::cos(E));
    double T = w + f;

    // Calculate position vector in ECI coordinates
    std::vector<double> r_vec = {
        r * (std::cos(Om) * std::cos(T) - std::sin(Om) * std::sin(T) * std::cos(inc)),
        r * (std::sin(Om) * std::cos(T) + std::cos(Om) * std::sin(T) * std::cos(inc)),
        r * (std::sin(T) * std::sin(inc))};

    // Calculate velocity vector
    std::vector<double> v_vec = {
        mu / h * (-(std::cos(Om) * (std::sin(T) + e * std::sin(w))) - std::sin(Om) * (std::cos(T) + e * std::cos(w)) * std::cos(inc)),
        mu / h * (-(std::sin(Om) * (std::sin(T) + e * std::sin(w))) + std::cos(Om) * (std::cos(T) + e * std::cos(w)) * std::cos(inc)),
        mu / h * ((std::cos(T) + e * std::cos(w)) * std::sin(inc))};

    return {r_vec, v_vec};
}

// Function to run APC propagator
void runAPCPropagator(double t0, double tf, std::vector<double> r0, std::vector<double> v0, int deg = 4)
{
    double mass = 212;        // sat mass (kg)
    double area = 10;         // sat wetted area (m^2)
    double reflectance = 1.5; // sat refelction absorption ratio
    double drag_C = 2.0;      // sat coefficient of drag

    // Perturbation calc flags
    bool compute_drag = false;       // atmostpheric drag toggle
    bool compute_SRP = false;        // Solar radiation pressure toggle
    bool compute_third_body = false; // Third body gravity toggle
    bool compute_hamiltonian = true; // whether or not the hamiltonian should be compuited for the output
    // Ephemeris
    string spk = "de440.bsp";
    string lsk = "naif0012.tls";
    list<string> bodies = {"SUN", "EARTH"};
    string center = "MOON";
    string frame = "J2000";

    Orbit orbit("MOON", "MOON_PA", "J2000");
    orbit.SetProperties(area, reflectance, mass, drag_C);
    orbit.SetPosition0(r0);
    orbit.SetVelocity0(v0);
    orbit.SetIntegrationTime(t0, tf);
    orbit.SetGravityApproximationDegree(deg);
    orbit.SinglePropagate();
}

int main()
{
    int numLongs = 8;
    int numIncs = 7;
    int numIterations = 100; // Define the number of iterations
    double numOrbits = 2;    // Define the number of orbits to propagate per iteration
    int maxDegree = 40;      // Define the maximum degree of the gravity field
    std::stringstream stream;
    if (std::fmod(numOrbits, 1) == 0)
    {
        stream << std::fixed << std::setprecision(2) << int(numOrbits);
    }
    else
    {
        stream << std::fixed << std::setprecision(2) << numOrbits;
    }
    std::string orbStr = stream.str();
    stream.str("");
    stream << std::fixed << std::setprecision(2) << numIterations;
    string iterStr = stream.str();
    stream.str("");
    stream << std::fixed << std::setprecision(2) << maxDegree;
    string maxDegreeStr = stream.str();
    // Define the filename
    string filename = "degree_trade_study_" + orbStr + "orbits_" + iterStr + "iterations_" + maxDegreeStr + "maxDeg.txt";
    // generate a list of longitude in radians from 0 to 2pi
    std::vector<double> longitudes;
    for (int i = 0; i < numLongs; i++)
    {
        longitudes.push_back(i * 2 * C_PI / numLongs);
    }
    // generate a list of inclinations in radians from 0 to 2pi (we want retro orbits too)
    std::vector<double> inclinations;
    for (int i = 0; i < numIncs; i++)
    {
        inclinations.push_back(i * C_PI / (numIncs - 1));
    }

    // Define orbital elements
    double e = 0.0;        // Eccentricity
    double w = 0.0;        // Argument of periapsis (radians)
    double M = 0.0;        // Mean anomaly (radians)
    double mu = C_MU_MOON; // Gravitational parameter (km^3/s^2)
    // Define the range of altitudes
    std::vector<double> altitudes = {30.0, 120.0};
    //log spaced degrees from 2 to maxDegree
    std::vector<double> gravitationalDegrees = {2,3,4,5,6,7,8,9,11,13,15,18,21,24,29,34,39};
    // Create a list of altitude-degree pairs
    std::vector<std::pair<double, int>> altitudeDegreePairs;
    for (double altitude : altitudes)
    {
        for (int gravitationalDegree : gravitationalDegrees)
        {
            altitudeDegreePairs.push_back(std::make_pair(altitude, gravitationalDegree));
        }
    }

    // Print the randomized order of altitude-degree pairs
    std::cout << "The trade study will be performed in this order: " << std::endl;
    for (const auto &pair : altitudeDegreePairs)
    {
        std::cout << "Altitude: " << pair.first << ", Gravitational Degree: " << pair.second << std::endl;
    }

    // Start the results file
    //  write results to file in table format with a header describing the orbital elements
    std::ofstream myfile;
    myfile.open(filename);
    // write full precision
    myfile.precision(16);
    // list orbital elements from their variable (except semi major axis) in header
    myfile << "min Semi-major Axis (km) = " << altitudes[0] + C_Rmoon << "\t";
    myfile << "max Semi-major Axis (km) = " << altitudes[altitudes.size() - 1] + C_Rmoon << "\t";
    myfile << "Argument of Periapsis (rad) = " << w << "\t";
    myfile << "Mean Anomaly (rad) = " << M << "\t";
    myfile << "Gravitational Parameter (km^3/s^2) = " << C_MU_MOON << "\t";
    myfile << "Number of Iterations = " << numIterations << "\t";
    myfile << "Number of Orbits = " << numOrbits << std::endl;
    // write column header for results table (alt, deg, longitude, inclination, avg time, std dev)
    myfile << "Altitude (km)\tGravitational Degree\tStarting Longitude (deg)\tInclination(deg)\tAverage Time (s)\tStandard Deviation (s)" << std::endl;
    myfile.close();
    int count = 0;

    auto progress_eta_start_time = std::chrono::high_resolution_clock::now();
    for (double longitude : longitudes)
    {
        // perform once per inclination
        for (double inclination : inclinations)
        {
            // last two orbital elements
            double Om = longitude;
            double inc = inclination;
            // Calculate the average computation time for each altitude and gravitational degree
            std::map<std::pair<double, int>, double> averageTimes;
            std::map<std::pair<double, int>, double> standardDeviations;
            // Shuffle the list of altitude-degree pairs to reduce bias
            std::random_device rd;
            std::mt19937 gen(rd());
            std::shuffle(altitudeDegreePairs.begin(), altitudeDegreePairs.end(), gen);
            for (const auto &pair : altitudeDegreePairs)
            {

                // Calculate the initial state
                // calculate the initial position and velocity
                double a = pair.first + C_Rmoon; // Semi-major axis (km)

                // Calculate position and velocity using elms2rv function
                std::vector<std::vector<double>> result = elms2rv(a, e, inc, Om, w, M, mu);
                std::vector<double> r_vec = result[0];
                std::vector<double> v_vec = result[1];

                // Calculate orbital period
                double T = 2 * C_PI * std::sqrt(std::pow(a, 3) / mu);

                double t0 = 0.0;           // Initial time (s)
                double tf = numOrbits * T; // Final time (s)
                // Perform all iterations for one specific altitude-degree pair
                double totalTime = 0.0;
                double times[numIterations];
                for (int i = 0; i < numIterations; i++)
                {
                    // Run one iteration of APC propagator
                    //  Measure computation time
                    auto start = std::chrono::high_resolution_clock::now();

                    // Run APC propagator
                    runAPCPropagator(t0, tf, r_vec, v_vec, pair.second);

                    // Measure computation time
                    auto end = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> duration = end - start;

                    // Accumulate total time in milliseconds
                    totalTime += duration.count();
                    times[i] = duration.count();
                }

                // Calculate average time in milliseconds
                double averageTime = totalTime / numIterations;
                // Calculate standard deviation
                double standardDeviation = 0.0;
                for (int i = 0; i < numIterations; i++)
                {
                    standardDeviation += std::pow((times[i] - averageTime), 2);
                }
                standardDeviation = std::sqrt(standardDeviation / (numIterations-1));

                // Store average time with altitude-degree pair as key
                averageTimes[pair] = averageTime;
                standardDeviations[pair] = standardDeviation;
                // iterate counter
                count++;
                // display progress
                // calculate average time to run one set
                auto progress_eta_cur_time = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> progress_eta = progress_eta_cur_time - progress_eta_start_time;
                double time_per_set = progress_eta.count() / count;
                double time_remaining = time_per_set * (altitudeDegreePairs.size() * numLongs * numIncs - count);
                int eta_hours = int(time_remaining / 3600);
                int eta_minutes = int((time_remaining - eta_hours * 3600) / 60);
                int eta_seconds = int(time_remaining - eta_hours * 3600 - eta_minutes * 60);
                std::cout << "Progress: " << count << "/" << altitudeDegreePairs.size() * numLongs * numIncs << " sets completed. ETA: " << std::setfill('0') << std::setw(2) << eta_hours << "h:" << eta_minutes << "m:" << eta_seconds << "s \t";
                // Output average time and std for pair
                std::cout << "Average time for altitude " << pair.first << " and degree " << pair.second << " is " << averageTime << " seconds with standard deviation " << standardDeviation << std::endl;
            }

            // Sort the list of altitude-degree pairs in altitude then degree order
            std::sort(altitudeDegreePairs.begin(), altitudeDegreePairs.end(), [](const auto &a, const auto &b)
                      {
        if (a.first != b.first)
            return a.first < b.first;
        return a.second < b.second; });

            // Append to results file
            myfile.open(filename, std::ios_base::app);
            for (const auto &pair : altitudeDegreePairs)
            {
                myfile << pair.first << "\t" << pair.second << "\t" << longitude * 180 / C_PI << "\t" << inclination * 180 / C_PI << "\t" << averageTimes[pair] << "\t" << standardDeviations[pair] << std::endl;
            }
            myfile.close();
        }
    }
    return 0;
}