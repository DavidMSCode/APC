/**
 * @file CCQ_tests.cpp
 * @brief Simples tests for the CCQ function outputs
 * 
 * This file contains tests for the CCQ function outputs
 */

#include "chebyshev.h"
#include "const.h"
#include "c_functions.h"
#include "clenshaw_curtis_ivpII.h"
#include "lsq_chebyshev_fit.h"

#include <iostream>
using namespace std;



int main(int argc, char **argv){
cout<<"Running CCQ_tests"<<endl;
// Test Chebyshev
int N = 5;
int M = 5;
int arg = 2;
// Compute Clenshaw-Curtis Quadrature Constant Matrices
std::vector<double> T2((M+1)*(N+1),0.0);
//memset( T2, 0.0, ((M+1)*(N+1)*sizeof(double)));
std::vector<double> P2((N+1)*N,0.0);
//memset( P2, 0.0, ((N+1)*N*sizeof(double)));
std::vector<double> T1((M+1)*N,0.0);
//memset( T1, 0.0, ((M+1)*N*sizeof(double)));
std::vector<double> P1(N*(N-1),0.0);
//memset( P1, 0.0, (N*(N-1)*sizeof(double)));
std::vector<double> Ta((M+1)*(N-1),0.0);
//memset( Ta, 0.0, ((M+1)*(N-1)*sizeof(double)));
std::vector<double> A((N-1)*(M+1),0.0);
//memset( A, 0.0, ((N-1)*(M+1)*sizeof(double)));
clenshaw_curtis_ivpII(N,M,T2,P2,T1,P1,Ta,A);

return 0;
}