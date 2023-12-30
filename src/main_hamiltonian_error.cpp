#include "GRGM1200b.h"
#include "const.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
using namespace std;
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
vector<double> geomSpace(double start, double end, int num)
{
    vector<double> result(num);
    double ratio = pow(end / start, 1.0 / (num - 1));
    for (int i = 0; i < num; i++)
    {
        result[i] = start * pow(ratio, i);
    }
    return result;
}
int main(int argc, char **argv)
{
    // Print help message with default values
    if (cmdOptionExists(argv, argv + argc, "-help") || cmdOptionExists(argv, argv + argc, "-h"))
    {
        cout << "Usage: ./main_hamiltonian_error [-alt_min <altitude in km>] [-alt_max <altitude in km>] [-steps <number of steps>] [-alts <alt1,alt2,...,altN>] [-body <earth or moon>] [-filename <filename>]" << endl;
        cout << "Default values:" << endl;
        cout << "alt_min = 200 km" << endl;
        cout << "alt_max = alt_min" << endl;
        cout << "steps = 10 if alt_max is different from alt_min, else 1" << endl;
        cout << "alts = [alt_min, alt_min + (alt_max - alt_min) / (steps - 1), ..., alt_max]" << endl;
        cout << "body = earth" << endl;
        cout << "filename = H_error_<body>_<alt_min>km_to_<alt_max>km_<steps>.csv" << endl;
        return 0;
    }
    // set default values
    double alt_min_default = 200;
    double alt_min;
    double alt_max;
    int steps_default = 10;
    int steps;
    double MU_default = C_MU_EARTH;
    double MU;
    double R_default = C_Req;
    double R;
    string filename;
    string body_default = "Earth";
    string body;
    // check for altitude argument
    if (cmdOptionExists(argv, argv + argc, "-alt_min"))
    {
        char *alt_arg = getCmdOption(argv, argv + argc, "-alt_min");
        alt_min = atof(alt_arg);
    }
    else
    {
        alt_min = alt_min_default;
    }
    if (cmdOptionExists(argv, argv + argc, "-alt_max"))
    {
        char *alt_arg = getCmdOption(argv, argv + argc, "-alt_max");
        alt_max = atof(alt_arg);
    }
    else
    {
        alt_max = alt_min;
    }
    // number of steps
    if (alt_max == alt_min)
    {
        steps = 1;
    }
    else if (cmdOptionExists(argv, argv + argc, "-steps"))
    {
        char *steps_arg = getCmdOption(argv, argv + argc, "-steps");
        steps = atoi(steps_arg);
    }
    else
    {
        steps = steps_default;
    }
    //check for alts argument
    vector<double> alts;
    if (cmdOptionExists(argv, argv + argc, "-alts"))
    {
        char *alts_arg = getCmdOption(argv, argv + argc, "-alts");
        char *p = strtok(alts_arg, ",");
        while (p != 0)
        {
            alts.push_back(atof(p));
            p = strtok(NULL, ",");
        }
        alt_min = alts[0];
        alt_max = alts[alts.size() - 1];
        steps = alts.size();
    }

    // check for body argument
    if (cmdOptionExists(argv, argv + argc, "-body"))
    {
        char *body_arg = getCmdOption(argv, argv + argc, "-body");
        if (strcmp(body_arg, "earth") == 0)
        {
            body = "Earth";
            MU = C_MU_EARTH;
            R = C_Req;
        }
        else if (strcmp(body_arg, "moon") == 0)
        {
            body = "Moon";
            MU = C_MU_MOON;
            R = C_Rmoon;
        }
        else
        {
            cout << "Invalid body argument. Defaulting to Earth." << endl;
            body = "Earth";
            MU = C_MU_EARTH;
            R = C_Req;
        }
    }
    else
    {
        body = body_default;
        MU = MU_default;
        R = R_default;
    }
    // check for filename argument
    if (cmdOptionExists(argv, argv + argc, "-filename"))
    {
        char *filename_arg = getCmdOption(argv, argv + argc, "-filename");
        filename = filename_arg;
    }
    else
    {
        // generate filename from arguments
        filename = "H_error_" + body + "_" + to_string(int(alt_min)) + "km_to_" + to_string(int(alt_max)) + "km_" + to_string(steps) + ".csv";
    }

    // Calculate altitudes vector
    vector<double> altitudes;
    if(cmdOptionExists(argv, argv + argc, "-alts")){
    altitudes=alts;
    }
    else if (steps == 1)
    {
        altitudes.push_back(alt_min);
    }
    else
    {
        for (int j = 0; j < steps; j++)
        {
            double alt = alt_min + (alt_max - alt_min) * j / (steps - 1);
            altitudes.push_back(alt);
        }
    }
    // Run for each altitude
    //  generate position offset vector from 1e-15 to 1 km
    vector<double> r_offset = geomSpace(1e-15, 1, 100);
    vector<vector<double>> H(steps, vector<double>(r_offset.size()));
    vector<vector<double>> dH(steps, vector<double>(r_offset.size()));
    for (int ii = 0; ii < altitudes.size(); ii++)
    {
        double alt = altitudes[ii];
        // Calulcate initial position and velocity
        double r0[3] = {alt + R, 0, 0};
        double v0 = sqrt(MU / r0[0]);
        double PE;
        GRGM1200bPot(r0, &PE, 70);
        double KE = 0.5 * v0 * v0;
        double H0 = KE + PE;

        for (int jj = 0; jj < r_offset.size(); jj++)
        {
            double dr = r_offset[jj];
            double r[3] = {r0[0] + dr, r0[1], r0[2]};
            GRGM1200bPot(r, &PE, 70);
            double KE = 0.5 * v0 * v0;
            H[ii][jj] = KE + PE;
            dH[ii][jj] = abs(H[ii][jj] - H0) / H0;
        }
    }
    // write to csv where the columns are the altitudes and the rows are the position offsets using full precision

    ofstream myfile;
    myfile.open(filename);
    // set precision
    myfile.precision(16);
    myfile << "Altitude Offset(km),";
    for (double alt : altitudes)
    {
        myfile << alt;
        if (alt < alt_max)
        {
            myfile << ",";
        }
    }
    myfile << endl;
    for (int i = 0; i < r_offset.size(); i++)
    {
        myfile << r_offset[i] << ",";

        for (int j = 0; j < altitudes.size(); j++)
        {
            myfile << dH[j][i];
            if (j < altitudes.size() - 1)
            {
                myfile << ",";
            }
        }
        myfile << endl;
    }
    myfile.close();

    return 0;
}