


#ifndef CLOHESSYWILTSHIRE_H
#define CLOHESSYWILTSHIRE_H

#include <vector>

using namespace std;

vector<double> ClohessyWiltshire(double t, vector<double> ds0, double n);

pair<vector<double>, vector<double>> orb2lvlh(vector<double> r, vector<double> v, vector<double> dr, vector<double> dv);

pair<vector<double>, vector<double>> lvlh2orb(vector<double> r, vector<double> v, vector<double> dr, vector<double> dv);

#endif // CLOHESSYWILTSHIRE_H
