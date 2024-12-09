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
std::vector<double> Ta((M+1)*(N+1),0.0);
std::vector<double> A((N+1)*(M+1),0.0);
lsq_chebyshev_fit(-1.0, N, M, Ta, A);

// print A
cout<<"A: "<<endl;
pretty_print_matrix(A,N+1,16);
cout<<"Ta: "<<endl;
pretty_print_matrix(Ta,M+1,16);
return 0;
}