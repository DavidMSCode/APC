#include <omp.h>
#include <stdlib.h>
#include <SpiceUsr.h>
#include <vector>
#include <fstream>
#include <utility>
#include <string>
#include <tuple>
#include <iostream>
#include <algorithm>
#include "GravTests.h"
#include "EGM2008.h"
#include "GRGM1200b.h"
#include "lunar_perturbed_gravity.h"
#include "APC.h"
#include "perturbed_gravity.h"

using namespace std;

//Function prin

// Get the command line arguments for given flag
char *getCmdOption(char **begin, char **end, const string &option)
{
    char **itr = find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}
// Check if command line argument exists
bool cmdOptionExists(char **begin, char **end, const string &option)
{
    return find(begin, end, option) != end;
}

int GravTests(int argc, char *argv[])
{
    GravMesh(argc, argv);
    LunarGravSample(argc, argv);
    return 0;
}


int LunarGravSample(int argc, char *argv[]){
    //This function samples points all over the moon at a specific altitude and outputs the gravity at those points
    //The points are sampled on the entire lunar sphere

    //Initialize constants
    double alt;
    double N;
    double gravDeg;
    string filename;
    //Default values
    double defaultalt = 100;                                // default alt coordinate
    int defaultgravDeg = 250;                               // default degree of gravity model
    int defaultN = 10;                                     // default number of sample points in one direction
    string defaultfilename = "/Users/davidstanley/Documents/Github/APC/GRGM1200b_250x250_100km.txt"; // default name of file to write error mesh to
    // check if any command line args were passed
    if (argc >1)
    {
        // check if alt coordinate provided, if so convert string to double
        alt = (cmdOptionExists(argv, argv + argc, "-alt") ? atof(getCmdOption(argv, argv + argc, "-alt")) : defaultalt);
        // check if degree of gravity model provided, if so convert string to int
        gravDeg = (cmdOptionExists(argv, argv + argc, "-deg") ? atoi(getCmdOption(argv, argv + argc, "-deg")) : defaultgravDeg);
        // check if number of sample points provided, if so convert string to int
        N = (cmdOptionExists(argv, argv + argc, "-N") ? atoi(getCmdOption(argv, argv + argc, "-N")) : defaultN);
        // check if name of file to write sample is provided, if so convert string to char*
        string filename = (cmdOptionExists(argv, argv + argc, "-o") ? getCmdOption(argv, argv + argc, "-o") : defaultfilename);
    }
    else
    {
        // if no command line args provided, use all the default values
        alt = defaultalt;
        gravDeg = defaultgravDeg;
        N = defaultN;
        filename = defaultfilename;
    }

    //generate latitude vector
    vector<double> lats(N, 0.0);
    for (int i = 0; i < N; i++)
    {
        lats[i] = (i - N / 2) * 180.0 / N;
    }
    //generate longitude vector
    vector<double> lons(N, 0.0);
    for (int i = 0; i < N; i++)
    {
        lons[i] = (i - N / 2) * 360.0 / N;
    }
    // Create mesh
    vector<vector<double>> latmesh(N, vector<double>(N, 0.0));
    vector<vector<double>> lonmesh(N, vector<double>(N, 0.0));
    tie(latmesh, lonmesh) = meshgrid(lats, lons); // std::tie() is used to unpack the tuple returned by meshgrid

    //Calculate just the full gravity vector at each point and store in a mesh
    vector<vector<vector<double>>> FullGravMesh(N, vector<vector<double>>(N, vector<double>(3, 0.0)));
    //also store the xyz coordinates of each point after converting from lunar PA coordinates to J2000 coordinates centered on the moon
    vector<vector<vector<double>>> xyzMesh(N, vector<vector<double>>(N, vector<double>(3, 0.0)));
    // Calculate the Full (up to degree GravDeg) and approximate acceleration surface
    string lskFile = "naif0012.tls";
    furnsh_c(lskFile.c_str());
    furnsh_c("moon_pa_de440_200625.bpc");
    furnsh_c("moon_de440_220930.tf");

    //convert coordiantes to J2000 using SPICE
    double et;
    SpiceDouble rot[3][3];
    //convert to ephemeris time
    utc2et_c("2003 Dec 19 16:48:00", &et);
    //convert to J2000
    pxform_c("MOON_PA", "J2000", et, rot);

    for (int ii = 0; ii < N; ii++)
    {
        for (int jj = 0; jj < N; jj++)
        {
            vector<double> xyz = latlonalt2xyz(latmesh[ii][jj], lonmesh[ii][jj], alt);
            // convert to array
            double arrxyz[3] = {xyz[0], xyz[1], xyz[2]};
            double G[3];
            GRGM1200b(arrxyz, G, gravDeg);
            FullGravMesh[ii][jj][0] = G[0];
            FullGravMesh[ii][jj][1] = G[1];
            FullGravMesh[ii][jj][2] = G[2];
            //rotate the coordinates to J2000
            double xyzJ2000[3];
            mxv_c(rot, arrxyz, xyzJ2000);
            xyzMesh[ii][jj][0] = xyzJ2000[0];
            xyzMesh[ii][jj][1] = xyzJ2000[1];
            xyzMesh[ii][jj][2] = xyzJ2000[2];
        }
    }
    kclear_c();

    //Output data to a txt file with one line per point with tab delimited values in equal columns knowing that each value will be 16 digits long and with a set width of 30 for each value
    ofstream myfile;
    myfile.open(filename);
    //print max double digits of precision
    myfile.precision(16);
    int width = 22;
    //Write a header line tab delim
    myfile << "All coordinates are in the J2000 frame centered on the moon with an altitude of approximately 100 km\n";
    myfile << "The kernels used for the Lunar orientation data are naif0012.tls, moon_pa_de440_200625.bpc, and moon_de440_220930.tf\n";
    myfile << "The gravity model used is GRGM1200b to degree and order 250\n";
    myfile << "UTC time for the sampling is 2003 Dec 19 16:48:00\n";
    myfile << std::setw(width)<<left<<"t (s J2000)" << "  ";
    myfile << std::setw(width)<<left<<"x (km)" << "  ";
    myfile << std::setw(width)<<left<<"y (km)" << "  ";
    myfile << std::setw(width)<<left<<"z (km)" << "  ";
    myfile << std::setw(width)<<left<<"Gx (km/s^2)" << "  ";
    myfile << std::setw(width)<<left<<"Gy (km/s^2)" << "  ";
    myfile << std::setw(width)<<left<<"Gz (km/s^2)" << "\n";

    for (int ii = 0; ii < N; ii++)
    {
        for (int jj = 0; jj < N; jj++)
        {
            myfile << std::setw(width)<<left<<et << "  ";
            myfile << std::setw(width)<<left<<xyzMesh[ii][jj][0] << "  ";
            myfile << std::setw(width)<<left<<xyzMesh[ii][jj][1] << "  ";
            myfile << std::setw(width)<<left<<xyzMesh[ii][jj][2] << "  ";
            myfile << std::setw(width)<<left<<FullGravMesh[ii][jj][0] << "  ";
            myfile << std::setw(width)<<left<<FullGravMesh[ii][jj][1] << "  ";
            myfile << std::setw(width)<<left<<FullGravMesh[ii][jj][2] << "\n";

        }
    }

    myfile.close();

    return 0;
}

int GravMesh(int argc, char *argv[])
{
    double lat;      // latitude coordinate
    double lon;      // longitude coordinate
    double alt;      // altitude coordinate
    int gravDeg;     // degree of gravity model
    int N;           // number of sample points in one direction
    double range;    // range of lat and lon coordinates for meshgrid
    string filename; // name of file to write error mesh to
    string body;     // body to use for gravity model

    double defaultlat = 0;                          // default lat coordinate (greenwich, UK)
    double defaultlon = 0;                           // default lon coordinate (greenwich, UK)
    double defaultalt = 100;                         // default alt coordinate
    int defaultgravDeg = 100;                        // default degree of gravity model
    int defaultN = 100;                              // default number of sample points in one direction
    double defaultrange = 1.0;                       // default range of lat and lon coordinates for meshgrid
    string defaultfilename = "EGM2008ErrorMesh.txt"; // default name of file to write error mesh to
    string defaultBody = "Earth";                    // default body to use for gravity model

    // check if any command line args were passed
    if (argc > 1)
    {
        // check if a lat coordinate provided, if so convert string to double
        lat = (cmdOptionExists(argv, argv + argc, "-lat") ? atof(getCmdOption(argv, argv + argc, "-lat")) : defaultlat);
        // check if a long coordinate provided, if so convert string to double
        lon = (cmdOptionExists(argv, argv + argc, "-lon") ? atof(getCmdOption(argv, argv + argc, "-lon")) : defaultlon);
        // check if alt coordinate provided, if so convert string to double
        alt = (cmdOptionExists(argv, argv + argc, "-alt") ? atof(getCmdOption(argv, argv + argc, "-alt")) : defaultalt);
        // check if degree of gravity model provided, if so convert string to int
        gravDeg = (cmdOptionExists(argv, argv + argc, "-deg") ? atoi(getCmdOption(argv, argv + argc, "-deg")) : defaultgravDeg);
        // check if number of sample points provided, if so convert string to int
        N = (cmdOptionExists(argv, argv + argc, "-N") ? atoi(getCmdOption(argv, argv + argc, "-N")) : defaultN);
        // check if range of lat and lon coordinates for meshgrid provided, if so convert string to double
        range = (cmdOptionExists(argv, argv + argc, "-range") ? atof(getCmdOption(argv, argv + argc, "-range")) : defaultrange);
        // check if name of file to write error mesh to provided, if so convert string to char*
        filename = (cmdOptionExists(argv, argv + argc, "-o") ? getCmdOption(argv, argv + argc, "-o") : defaultfilename);
        // check if body to use for gravity model provided, if so store as string
        body = (cmdOptionExists(argv, argv + argc, "-body") ? getCmdOption(argv, argv + argc, "-body") : defaultBody);
    }
    else
    {
        // if no command line args provided, use all the default values
        lat = defaultlat;
        lon = defaultlon;
        alt = defaultalt;
        gravDeg = defaultgravDeg;
        N = defaultN;
        range = defaultrange;
        filename = defaultfilename;
        body = defaultBody;
    }
    // Assign full gravity function to EGM2008 if Earth or GRGM1200b if Moon
    // Assign approx gravity function to hardcoded J2-J6 if Earth or GRGM1200b 4x4 if Moon
    auto FullGrav = [gravDeg, body](double *p, double *Gxyz)
        {
            if (body == "Earth"){
            EGM2008(p, Gxyz, gravDeg);
            }
            else if (body == "Moon"){
            GRGM1200b(p, Gxyz, gravDeg);
            }
        };
    auto ApproxGrav = [body](double *p, double *Gxyz)
        {
            if (body == "Earth"){
                Grav_EarthJ2J6(p, Gxyz);
            }
            else if (body == "Moon"){
                GRGM1200b(p, Gxyz, 4);
            }
        };


    // remove any file extension from filename and store in type string
    string type = filename.substr(filename.find_last_of("."));
    filename = filename.substr(0, filename.find_last_of("."));

    // append filename with "lats" to indicate that it contains the latitudes used to generate the error mesh
    string latsfilename = filename;
    latsfilename.append("lats");
    // append filename with "lons" to indicate that it contains the longitudes used to generate the error mesh
    string lonsfilename = filename;
    lonsfilename.append("lons");
    // append all filenames with filetype
    filename.append(type);
    latsfilename.append(type);
    lonsfilename.append(type);

    // This function takes a latitude and longitude and a degree range and returns the differential gravity surface in km/s^2

    // Calculate reference acceleration at the center point
    vector<double> xyz = latlonalt2xyz(lat, lon, alt);
    // convert to array
    double arrxyz[3] = {xyz[0], xyz[1], xyz[2]};
    double refG[3];
    FullGrav(arrxyz, refG); // The "full" gravity at the central point
    double refGmag = sqrt(refG[0] * refG[0] + refG[1] * refG[1] + refG[2] * refG[2]);
    // Print the reference position xyz and acceleration
    cout << "refX: " << arrxyz[0] << endl;
    cout << "refY: " << arrxyz[1] << endl;
    cout << "refZ: " << arrxyz[2] << endl;
    cout << "refGmag: " << refGmag << endl;
    // Calculate low fidelity gravity with J2-J6 terms
    double refDelG[3];
    double refLowG[3];
    ApproxGrav(arrxyz, refLowG); // The low order gravity at the central point
    for (int i = 0; i < 3; i++)
    {
        refDelG[i] = refG[i] - refLowG[i]; // The high order gravity at the central point (no low order terms)
    }

    // lats and lons
    vector<double> lats(N, 0.0);
    vector<double> lons(N, 0.0);
    for (int i = 0; i < N; i++)
    {
        lats[i] = lat + (i - N / 2) * range / N;
        lons[i] = lon + (i - N / 2) * range / N;
    }
    // Create mesh
    vector<vector<double>> latmesh(N, vector<double>(N, 0.0));
    vector<vector<double>> lonmesh(N, vector<double>(N, 0.0));
    tie(latmesh, lonmesh) = meshgrid(lats, lons); // std::tie() is used to unpack the tuple returned by meshgrid

    // Calculate the Full (up to degree GravDeg) and approximate acceleration surface
    vector<vector<double>> FullGravMesh(N, vector<double>(N, 0.0));
    // And the approximate gravity mesh using the low order terms from the central point
    vector<vector<double>> ApproxGravMesh(N, vector<double>(N, 0.0));
    // and the difference between the two as an error mesh
    vector<vector<double>> ErrorGravMesh(N, vector<double>(N, 0.0));
    // Calculate the error mesh multithreaded
    omp_set_num_threads(omp_get_max_threads());
    #pragma omp parallel for collapse(2)
    for (int ii = 0; ii < N; ii++)
    {
        for (int jj = 0; jj < N; jj++)
        {
            //output the thread #,ii and jj values to stdout
            printf("Thread %d: ii = %d, jj = %d\n", omp_get_thread_num(), ii, jj);
            vector<double> xyz = latlonalt2xyz(latmesh[ii][jj], lonmesh[ii][jj], alt);
            // convert to array
            double arrxyz[3] = {xyz[0], xyz[1], xyz[2]};
            double G[3];
            FullGrav(arrxyz, G);
            FullGravMesh[ii][jj] = sqrt(G[0] * G[0] + G[1] * G[1] + G[2] * G[2]);
            // Get low order terms at current coordinate
            double LowG[3];
            ApproxGrav(arrxyz, LowG);
            // Add low order terms to high order terms from central point
            double approxG[3];
            for (int i = 0; i < 3; i++)
            {
                approxG[i] = refDelG[i] + LowG[i]; // approximate acceleration vector
            }
            // store approximate magnitude
            ApproxGravMesh[ii][jj] = sqrt(approxG[0] * approxG[0] + approxG[1] * approxG[1] + approxG[2] * approxG[2]);
            // store error (length of error vector)
            ErrorGravMesh[ii][jj] = sqrt((G[0] - approxG[0]) * (G[0] - approxG[0]) + (G[1] - approxG[1]) * (G[1] - approxG[1]) + (G[2] - approxG[2]) * (G[2] - approxG[2]));
        }
    }

    // Write error mesh to filename using double precision floating point format with 16 digits of precision using cout and ifstream
    ofstream myfile;
    myfile.open(filename);
    myfile.precision(16);
    for (int ii = 0; ii < N; ii++)
    {
        for (int jj = 0; jj < N; jj++)
        {
            myfile << ErrorGravMesh[ii][jj] << " ";
        }
        myfile << "\n";
    }
    myfile.close();
    // write latmesh to latsfilename using double precision floating point format with 16 digits of precision using cout and ifstream
    myfile.open(latsfilename);
    myfile.precision(16);
    for (int ii = 0; ii < N; ii++)
    {
        for (int jj = 0; jj < N; jj++)
        {
            myfile << latmesh[ii][jj] << " ";
        }
        myfile << "\n";
    }
    myfile.close();
    // write lonmesh to lonsfilename using double precision floating point format with 16 digits of precision using cout and ifstream
    myfile.open(lonsfilename);
    myfile.precision(16);
    for (int ii = 0; ii < N; ii++)
    {
        for (int jj = 0; jj < N; jj++)
        {
            myfile << lonmesh[ii][jj] << " ";
        }
        myfile << "\n";
    }
    myfile.close();

    return 0;
}

// vector<double> latlonalt2xyz(double lat, double lon, double alt)
// {
//     // This function takes latitude longitude and altitude over Earth and returns the x y z coordinates in km

//     // Initialize constants
//     double a = 6378.1370;                                   // kilometers
//     double b = 6356.7523142;                                // kilometers
//     double e = sqrt(1.0 - (b * b) / (a * a));               // eccentricity
//     double N = a / sqrt(1.0 - e * e * sin(lat) * sin(lat)); // radius of curvature
//     double x = (N + alt) * cos(lat) * cos(lon);             // x coordinate
//     double y = (N + alt) * cos(lat) * sin(lon);             // y coordinate
//     double z = (N * (1.0 - e * e) + alt) * sin(lat);        // z coordinate

//     vector<double> xyz(3, 0.0); // initialize vector to store x y z
//     xyz[0] = x;
//     xyz[1] = y;
//     xyz[2] = z;

//     return xyz;
// }
vector<double> latlonalt2xyz(double lat, double lon, double alt)
{
    double R = 6378.1370+alt; // radius of position
    lat = lat * M_PI / 180.0; // convert to radians
    lon = lon * M_PI / 180.0; // convert to radians

    //This function takes lat lona nd altitude around the earth and returns the x y z coordinates in km assuming a spherical earth

    double x = R * cos(lat) * cos(lon);
    double y = R * cos(lat) * sin(lon);
    double z = R * sin(lat);

    vector<double> xyz(3, 0.0); // initialize vector to store x y z
    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;

    return xyz;
}

tuple<vector<vector<double>>, vector<vector<double>>> meshgrid(vector<double> x, vector<double> y)
{
    // This function takes two vectors and returns two matrices that are the result of a meshgrid operation
    // The first matrix is the result of the meshgrid operation on the first vector
    // The second matrix is the result of the meshgrid operation on the second vector

    int m = x.size();
    int n = y.size();

    vector<vector<double>> X(m, vector<double>(n, 0.0));
    vector<vector<double>> Y(m, vector<double>(n, 0.0));

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            X[i][j] = x[i];
            Y[i][j] = y[j];
        }
    }

    return make_pair(X, Y);
}