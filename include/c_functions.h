/*
 *  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
 *  DATE WRITTEN:     Feb 2017
 *  LAST MODIFIED:    Feb 2017
 *  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
 *  DESCRIPTION:      Header file
 */

#ifndef _C_FUNC_H_
#define _C_FUNC_H_
#include <vector>
#include <iostream>
// Macro for looking up array indices in a 2-d array, using "1"-based indexing (Brent Macomber).
// INPUT: i = desired row number, j = desired col number, ld = number of rows
#define ID2(i, j, ld) ((((j) - 1) * (ld)) + ((i) - 1))

// Cross Product in 3D
// INPUT:  vector a & vector b
// OUTPUT: vector c, (c = axb)

/**
 * @brief Takes the 3d cross product of a and b and returns it as c
 *
 * @param a a 3 array of doubles
 * @param b a 3 array of doubles
 * @param c output: a 3 array of doubles
 */
void cross_product_3D(double *a, double *b, double *c);

/**
 * @brief Takes a 3d vector and returns its L2 norm
 *
 * @param a a 3 array of doubles
 * @param n the length of a
 */
void Cnorm(double *a, double &n);

/**
 * @brief returns the inner dot product of a and b
 *
 * @param a a 3 array of doubles
 * @param b a 3 array of doubles
 * @return double the resulting dot product
 */
double Cdot(double *a, double *b);
/**
 * @brief returns the inner dot product of a and b
 * 
 * @param a a 3 vector of doubles
 * @param b a 3 vector of doubles
 * @return double the resulting dot product
 */
double Cdot(std::vector<double> a, std::vector<double> b);

void find_less(double *a, double l, int len, double *ind);

void find_less_equal(double *a, double *le, int len, double *ind);

void find_more(double *a, double m, int len, double *ind);

void find_more_equal(double *a, double me, int len, double *ind);

void find_equal(double *a, double e, int len, double *res);

void Cmax(double *a, int size, double *max);

void Cmin(double *a, int size, double *min);

void Cabs(double *a, int size, double *b);

void matadd(double *a, double *b, int size, double *c);

void matsub(double *a, double *b, int size, double *c);

/*!
 * \brief Matrix Multiplication
 * This is a simple matrix multiplication function
 *
 * \param[in] A Vector representation of matrix A (size m x n)
 * \param[in] B Vector representation of matrix B (size n x q)
 * \param[in] m Column dimension of A
 * \param[in] n Shared dimension of A and B
 * \param[in] q Row dimension of B
 * \param[out] B Matrix Output (size m x q)
 */

std::vector<double> matmul(std::vector<double> A, std::vector<double> B,
                           const int m, const int n, const int q,
                           const int ldA, const int ldB);

void pretty_print_matrix(std::vector<double> A, int lr, int prec=6);

#endif
