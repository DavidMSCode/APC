#ifndef F521A02F_A4CB_4CB3_BD71_E894A1D1B26A
#define F521A02F_A4CB_4CB3_BD71_E894A1D1B26A



using namespace std;

int GravTests(int argc, char *argv[]);
int LunarGravSample(int argc, char *argv[]);

int GravMesh(int argc, char *argv[]);
vector<double> latlonalt2xyz(double lat, double lon, double alt);
tuple<vector<vector<double>>, vector<vector<double>>> meshgrid(vector<double> x, vector<double> y);
#endif /* F521A02F_A4CB_4CB3_BD71_E894A1D1B26A */
