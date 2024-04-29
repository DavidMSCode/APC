/*
 *  AUTHORS:          David Stanley (DavidMS4@Illinois.edu)
 *  DATE WRITTEN:     Nov 2023
 * @ Modified by: Your name
 * @ Modified time: 2023-12-24 17:44:44
 *  DESCRIPTION:      Set up an Adaptive-Picard-Chebyshev integration test case
 *  REFERENCE:        Woollands, R., and Junkins, J., "Nonlinear Differential Equation Solvers
 *                    via Adaptive Picard-Chebyshev Iteration: Applications in Astrodynamics", JGCD, 2016.
 */
#include <string>
#include <iostream>

#include <vector>
#include <iostream>
#include <utility>
#include <unistd.h>
#include <fstream>
#include <EGM2008.h>
#include <GRGM1200b.h>
#include <omp.h>
#include <time.h>
#include <math.h>
#include <algorithm>
using namespace std;

vector<double> randomP(double r_min, double r_max)
{
    // Generate a random v3 position around a sphere
    // set random seed
    double r = r_min + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX / (r_max - r_min)));
    double theta = static_cast<double>(rand()) / (static_cast<double>(RAND_MAX / (2 * M_PI)));
    double phi = static_cast<double>(rand()) / (static_cast<double>(RAND_MAX / (2 * M_PI)));
    vector<double> pos = {r * cos(theta) * sin(phi), r * sin(theta) * sin(phi), r * cos(phi)};
    return pos;
}

vector<vector<double>> randomP(int n, double r_min, double r_max)
{
    srand(time(NULL));
    // Generate n random positions around a sphere
    vector<vector<double>> pos;
    for (int i = 0; i < n; i++)
    {
        pos.push_back(randomP(r_min, r_max));
    }
    return pos;
}

vector<bool> compareGrav(vector<vector<double>> grav1, vector<vector<double>> grav2)
{
    // Compare the gravity vectors
    vector<bool> results;
    for (int i = 0; i < grav1.size(); i++)
    {
        results.push_back(grav1[i] == grav2[i]);
    }
    return results;
}


int main(int argc, char **argv)
{
    int deg = 200;
    // generate random test points around the Earth
    auto Earth_test_points = randomP(40, 6372.0, 6371.0 + 1000.0);
    // empty grav vector of size Earth_test_points,3
    vector<vector<double>> grav;
    grav.resize(Earth_test_points.size());
    for (int i = 0; i < Earth_test_points.size(); i++)
    {
        grav[i].resize(3);
    }
    // start timer
    double start = omp_get_wtime();
    // for each test point calculate the gravity
    for (int i = 0; i < Earth_test_points.size(); i++)
    {
        double Grav[3];
        // calculate the gravity at the test point
        GRGM1200b(&Earth_test_points[i][0], Grav, deg);
        // store grav in a vector
        grav[i][0] = Grav[0];
        grav[i][1] = Grav[1];
        grav[i][2] = Grav[2];
    }
    // end timer
    double end = omp_get_wtime();
    double elapsed_secs = end - start;
    cout << "Time taken for serial run 1: " << elapsed_secs << "s" << endl;

    // Run again and compare the results
    vector<vector<double>> grav2;
    grav2.resize(Earth_test_points.size());
    for (int i = 0; i < Earth_test_points.size(); i++)
    {
        grav2[i].resize(3);
    }
    //start timer
    double start2 = omp_get_wtime();
    for (int i = 0; i < Earth_test_points.size(); i++)
    {
        double Grav[3];
        // calculate the gravity at the test point
        GRGM1200b(&Earth_test_points[i][0], Grav,deg);
        //store grav in a vector
        grav2[i] = {Grav[0], Grav[1], Grav[2]};
    }
    //end timer
    double end2 = omp_get_wtime();
    double elapsed_secs2 = double(end2 - start2);
    cout << "Time taken for serial run 2: " << elapsed_secs2 << "s" << endl;
    // compare the results
    vector<bool> results = compareGrav(grav, grav2);
    //print pass if all results are true
    if(std::all_of(results.begin(), results.end(), [](bool v) { return v; }))
    {
        cout << "All single tests passed" << endl;
    }
    else
    {
        cout << "Some single tests failed" << endl;
    }


    //Run the loop again but in parallel with openmp

    vector<vector<double>> grav3;
    //start timer
    double start3 = omp_get_wtime();
    #pragma omp parallel
    {
        std::vector<std::vector<double>> grav3_private;
        #pragma omp for nowait schedule(static)
        for (int i = 0; i < Earth_test_points.size(); i++)
        {

            double Grav[3];
            // calculate the gravity at the test point
            GRGM1200b(&Earth_test_points[i][0], Grav,deg);
            grav3_private.push_back({Grav[0], Grav[1], Grav[2]});
        }
        #pragma omp for schedule(static) ordered
        for (int i = 0; i < omp_get_num_threads(); i++)
        {
            #pragma omp ordered
            grav3.insert(grav3.end(), grav3_private.begin(), grav3_private.end());
        }
    }
    //end timer
    double end3 = omp_get_wtime();
    double elapsed_secs3 = double(end3 - start3);
    cout << "Time taken for parallel run: " << elapsed_secs3 << "s" << endl;

    //repeat parallel but a simple pre allocated vector
    vector<vector<double>> grav4(Earth_test_points.size(),vector<double>(3));
    //start timer
    double start4 = omp_get_wtime();
    #pragma omp parallel for
    for (int i = 0; i < Earth_test_points.size(); i++)
    {
        double Grav[3];
        // calculate the gravity at the test point
        GRGM1200b(&Earth_test_points[i][0], Grav,deg);
        grav4[i] = {Grav[0], Grav[1], Grav[2]};
    }
    //end timer
    double end4 = omp_get_wtime();
    double elapsed_secs4 = double(end4 - start4);
    cout << "Time taken for parallel run with preallocated vector: " << elapsed_secs4 << "s" << endl;


    //repear parallel but with a preallocated array
    double grav5[Earth_test_points.size()][3];
    //start timer
    double start5 = omp_get_wtime();
    #pragma omp parallel for
    for (int i = 0; i < Earth_test_points.size(); i++)
    {
        double Grav[3];
        // calculate the gravity at the test point
        GRGM1200b(&Earth_test_points[i][0], Grav,deg);
        grav5[i][0] = Grav[0];
        grav5[i][1] = Grav[1];
        grav5[i][2] = Grav[2];
    }
    //end timer
    double end5 = omp_get_wtime();
    double elapsed_secs5 = double(end5 - start5);
    cout << "Time taken for parallel run with preallocated array: " << elapsed_secs5 << "s" << endl;


    // compare the results
    vector<bool> results3 = compareGrav(grav, grav4);

        // compare the results
    vector<bool> results2 = compareGrav(grav, grav3);

    // compare the results
    vector<vector<double>> grav5_vec(Earth_test_points.size(),vector<double>(3));
    for (int i = 0; i < Earth_test_points.size(); i++)
    {
        grav5_vec[i][0] = grav5[i][0];
        grav5_vec[i][1] = grav5[i][1];
        grav5_vec[i][2] = grav5[i][2];
    }
    vector<bool> results4 = compareGrav(grav, grav5_vec);

    //print pass if all results are true
    if (std::all_of(results2.begin(), results2.end(), [](bool v) { return v; }))
    {
        cout << "All parallel tests passed" << endl;
    }
    else
    {
        cout << "Some parallel tests failed" << endl;
    }

    //print pass if all results are true
    if (std::all_of(results3.begin(), results3.end(), [](bool v) { return v; }))
    {
        cout << "All parallel tests with preallocated vector passed" << endl;
    }
    else
    {
        cout << "Some parallel tests with preallocated vector failed" << endl;
    }

    //print pass if all results are true
    if (std::all_of(results4.begin(), results4.end(), [](bool v) { return v; }))
    {
        cout << "All parallel tests with preallocated array passed" << endl;
    }
    else
    {
        cout << "Some parallel tests with preallocated array failed" << endl;
    }

    //print percent speedup range
    double speedup1 = elapsed_secs / elapsed_secs3;
    double speedup2 = elapsed_secs2 / elapsed_secs3;
    double speedup3 = elapsed_secs / elapsed_secs4;
    double speedup4 = elapsed_secs2 / elapsed_secs4;
    double speedup5 = elapsed_secs / elapsed_secs5;
    double speedup6 = elapsed_secs2 / elapsed_secs5;

    cout << "Speedup range: " << min({speedup1,speedup2,speedup3,speedup4, speedup5, speedup6}) << " to " << max({speedup1,speedup2,speedup3,speedup4,speedup5,speedup6}) << " times faster." << endl;


    return 0;
}
