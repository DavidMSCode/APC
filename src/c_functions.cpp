/*
 *  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
 *  DATE WRITTEN:     May 2017
 *  LAST MODIFIED:    May 2017
 *  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
 *  DESCRIPTION:      Performs some simple vector-matrix operations
 */

#include <math.h>

#include <vector>
#include <iostream>

#include "c_functions.h"
// Cross Product in 3D
// INPUT:  vector a & vector b
// OUTPUT: vector c, (c = axb)
void cross_product_3D(double *a, double *b, double *c)
{
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}
// Dot Product in 3D
// INPUT:  vector a & vector b
// OUTPUT: dot product
double Cdot(double *a, double *b)
{
  double product = 0.0;
  for (int i = 0; i < 3; i++)
  {
    product += a[i] * b[i];
  }
  return product;
}
/**
 * @brief returns the inner dot product of a and b
 *
 * @param a a 3 vector of doubles
 * @param b a 3 vector of doubles
 * @return double the resulting dot product
 *
 */
double Cdot(std::vector<double> a, std::vector<double> b)
{
  double product = 0.0;
  for (int i = 0; i < 3; i++)
  {
    product += a[i] * b[i];
  }
  return product;
}

// Length of 3d vector
// INPUT:  vector a
// OUTPUT: n
void Cnorm(double *a, double &n)
{
  n = sqrt(pow(a[0], 2) + pow(a[1], 2) + pow(a[2], 2));
}

// Maximum of a 2D array
// INPUT: vector a, size of vector a
// OUTPUT: max value in vector a
void Cmax(double *a, int size, double *max)
{
  *max = a[0];
  for (int i = 0; i <= size; i++)
  {
    if (a[i] > *max)
    {
      *max = a[i];
    }
  }
}

// Minimum of a 2D array
// INPUT: vector a, size of vector a
// OUTPUT: min value in vector a
void Cmin(double *a, int size, double *min)
{
  *min = a[0];
  for (int i = 0; i <= size; i++)
  {
    if (a[i] < *min)
    {
      *min = a[i];
    }
  }
}

/*!
 * \brief Matrix Multiplication (Brent Macomber)
 * This is a simple matrix multiplication function
 *
 * \param[in] A Vector representation of matrix A (size m x n)
 * \param[in] B Vector representation of matrix B (size n x q)
 * \param[in] m Column dimension of A
 * \param[in] n Shared dimension of A and B
 * \param[in] q Row dimension of B
 * \param[out] C Matrix Output (size m x q)
 */
std::vector<double> matmul(std::vector<double> A, std::vector<double> B,
                           const int m, const int n, const int q,
                           const int ldA, const int ldB)
{
  int Ai;
  int Bi;
  double Av;
  double Bv;
  std::vector<double> C(m * q, 0.0);
  const int ldC = ldA;
  // Stand Alone Method
  double sum;
  int ii, jj, kk;
  for (ii = 0; ii < m; ii++)
  {
    for (jj = 0; jj < q; jj++)
    {
      sum = 0.0;
      for (kk = 0; kk < n; kk++)
      {
        Ai = ID2(ii + 1, kk + 1, ldA);
        Bi = ID2(kk + 1, jj + 1, ldB);
        Av = A[Ai];
        Bv = B[Bi];
        sum += A[ID2(ii + 1, kk + 1, ldA)] * B[ID2(kk + 1, jj + 1, ldB)];
      }
      C[ID2(ii + 1, jj + 1, ldC)] = sum;
    } // end for jj
  } // end for ii
  return C;
}
