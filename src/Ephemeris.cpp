/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com), Alex Pascarella (alexp3@illinois.edu)
*  DATE WRITTEN:     Dec 2020
*  LAST MODIFIED:    Feb 2021
*  AFFILIATION:      University of Illinois at Urbana-Champaign (UIUC)
*  DESCRIPTION:      Ephemeris computation using Chebyshev polynomials
*/


#include <stdexcept>
#include <vector>
#include <string>
#include <map>
#include <list>

#include "SpiceUsr.h"

#include "Ephemeris.hpp"
#include "c_functions.h"
#include "c_functions2.hpp"
#include "const.h"
using namespace std;

ChebyshevEphemeris::ChebyshevEphemeris(){}
ChebyshevEphemeris::ChebyshevEphemeris(string spkFile,
                                       string lskFile,
                                       double initialEpoch,
                                       double finalEpoch,
                                       string body,
                                       string center,
                                       string frame)
{
    this->initialEpoch = initialEpoch;
    this->finalEpoch = finalEpoch;
    this->body = body;
    this->center = center;
    this->frame = frame;


    if ( (finalEpoch - initialEpoch) < (86400. * 300.) )
    {
        this->segmentLength = finalEpoch - initialEpoch;
    }

    //TODO: allow user to select the segment length
    else
    {
        this->segmentLength = 86400. * 300.;
    }

    // Number of ephemeris segments
    int segs = (int) ceil((finalEpoch - initialEpoch)/segmentLength);

    // Number of coefficients
    N = (int) (segmentLength/30/86400) * 40;
    if (N == 0){ N=40; }

    // Ensure that only one thread at a time accesses Spice
    #pragma omp critical
    {
        double t0 = initialEpoch;
        double t1;
        vector<double> coeff;

        furnsh_c (spkFile.c_str());
        furnsh_c (lskFile.c_str());

        for (int idx=0; idx<segs; idx++)
        {
            t1 = t0 + segmentLength;
            w1.push_back( (t1 + t0)/2.0 );
            w2.push_back( (t1 - t0)/2.0 );
            coeff.resize((N+1)*6);
            fill(coeff.begin(), coeff.end(), 0.0);

            // Chebyshev Least Squares Matrices
            double s = -1.0;

            // Least Squares Operator (A)
            std::vector<double> Ta((N+1)*(N+1),0.0);
            // memset( Ta, 0.0, ((N+1)*(N+1)*sizeof(double)));
            std::vector<double> A((N+1)*(N+1),0.0);
            // memset( A, 0.0, ((N+1)*(N+1)*sizeof(double)));

            fitChebyshevPolynomials(s, N, N, Ta, A);

            std::vector<double> tau(N+1,0.0);
            // memset( tau, 0.0, ((N+1)*sizeof(double)));
            std::vector<double> time((N+1),0.0);
            // memset( time, 0.0, ((N+1)*sizeof(double)));

            for (int i=0; i<=N; i++)
            {
                tau[i]   = -cos(i*C_PI/N);
                time[i] = tau[i] * w2[idx] + w1[idx];
            }

            // Compute States
            double state[6] = {0.0};
            double owlt = 0.0;
            std::vector<double> X(6*(N+1),0.0);
            // memset( X, 0.0, (6*(N+1)*sizeof(double)));

            for (int i=0; i<=N; i++)
            {
                spkezr_c(body.c_str(), time[i], frame.c_str(), "LT+S", center.c_str(), state, &owlt);
                for (int j=0; j<=5; j++)
                {
                    X[ID2(i+1,j+1,N+1)] = state[j];
                }
            }

            // Compute Chebyshev Coefficients
            coeff = matmul(A, X, N+1, N+1, 6, N+1, N+1);
            this->coefficients.push_back(coeff);
            t0 = t1;
        }

        unload_c(spkFile.c_str());
        unload_c(lskFile.c_str());
    }
}

// Getters
string ChebyshevEphemeris::getBody(){ return body; }
string ChebyshevEphemeris::getCenter(){ return center; }
string ChebyshevEphemeris::getFrame(){ return frame; }
double ChebyshevEphemeris::getInitialEpoch(){ return initialEpoch; }
double ChebyshevEphemeris::getFinalEpoch(){ return finalEpoch; }


void ChebyshevEphemeris::getState(double* state, double epoch)
{
    // Check that the query epoch is within the allowed bounds
    if ( (epoch < initialEpoch) || (epoch > finalEpoch) )
    {
        string error_msg = "Epoch " + std::to_string(epoch) + " is not within the allowed epoch range (" +
                            std::to_string(initialEpoch) + ", " + to_string(finalEpoch) + ")";
        throw std::invalid_argument( error_msg );
    }

    // Compute segment index
    int idx = (int) floor((epoch - initialEpoch)/segmentLength);

    // Compute tau
    double tau = (epoch - w1[idx])/w2[idx];

    // Interpolate Chebyshev Polynomials
    std::vector<double> T((N+1),0.0);
    // memset( T, 0.0, ( (N+1)*sizeof(double) ) );

    for (int i=0; i<=N; i++)
    {
        T[ID2(1,i+1,1)] = cos(i*acos(tau));
    }

    // Compute state
    state = &matmul(T, coefficients[idx], 1, N+1, 6, 1, N+1)[0];
}

void ChebyshevEphemeris::fitChebyshevPolynomials(double s, int N, int M, std::vector<double> &T, std::vector<double> &A)
{
    // Generate Chebyshev Polyniomials
    getChebyshevPolynomials(s, N, M, 2, T);

    // Weight Matrix
    std::vector<double> W((M+1)*(M+1),0.0);
    // memset( W, 0.0, ( (M+1)*(M+1)*sizeof(double) ) );
    for (int i=1; i<=M+1; i++)
    {
        for (int j=1; j<=M+1; j++)
        {
            if (i == j)
            {
                W[ID2(i,j,M+1)] = 1.0;
            }
        }
    }

    W[ID2(1,1,M+1)] = 0.5;
    W[ID2(M+1,M+1,M+1)] = 0.5;

    // V Matrix
    std::vector<double> V((N+1)*(N+1),0.0);
    // memset( V, 0.0, ((N+1)*(N+1)*sizeof(double)));
    for (int i=1; i<=N+1; i++)
    {
        for (int j=1; j<=N+1; j++)
        {
            if (i == j)
            {
                V[ID2(i,j,N+1)] = 2.0/M;
            }
        }
    }

    if (M == N)
    {
        V[ID2(1,1,N+1)] = 1.0/M;
        V[ID2(N+1,N+1,N+1)] = 1.0/M;
    }

    if (M > N)
    {
        V[ID2(1,1,N+1)] = 1.0/M;
    }

    // T Transpose
    std::vector<double> TT((N+1)*(M+1),0.0);
    // memset( TT, 0.0, ((N+1)*(M+1)*sizeof(double)));
    for (int i=1; i<=N+1; i++)
    {
        for (int j=1; j<=M+1; j++)
        {
            TT[ID2(i,j,N+1)] = T[ID2(j,i,M+1)];
        }
    }

    // Least Squares Operator
    std::vector<double> TTW((M+1)*(M+1),0.0);
    // memset( TTW, 0.0, ((M+1)*(M+1)*sizeof(double)));
    
    TTW = matmul(TT,W,N+1,M+1,M+1,N+1,M+1);
    A = matmul(V,TTW,N+1,N+1,M+1,N+1,N+1);
}

void ChebyshevEphemeris::getChebyshevPolynomials(double s, int N, int M, int arg, std::vector<double> &T)
{

  // Cosine Sample Points
  std::vector<double> tau((M+1),0.0);
//   memset( tau, 0.0, ((M+1)*sizeof(double)));
  for (int i=0; i<=M; i++){
    tau[i] = s*cos(i*C_PI/M);
  }

  if (arg == 1){
    // Chebyshev Polynomials (Recursive Formulation)
    for (int j=1; j<=N+1; j++){
      for (int i=1; i<=M+1; i++){
        if (j == 1){
          T[ID2(i,j,M+1)] = 1.0;
        }
        if (j == 2){
          T[ID2(i,j,M+1)] = tau[i-1];
        }
        if (j > 2){
          T[ID2(i,j,M+1)] = 2*tau[i-1]*T[ID2(i,j-1,M+1)] - T[ID2(i,j-2,M+1)];
        }
      }
    }
  }

  if (arg == 2){
    // Chebyshev Polyniomials (Trigonometric Formulation)
    for (int i=0; i<=M; i++){
      for (int j=0; j<=N; j++){
        T[ID2(i+1,j+1,M+1)] = cos(j*acos(tau[i]));
      }
    }
  }

}

EphemerisManager::EphemerisManager(){}
EphemerisManager::EphemerisManager(string spkFile,
                                   string lskFile,
                                   double initialEpoch,
                                   double finalEpoch,
                                   list<string> bodies,
                                   string center,
                                   string frame)
{
    for (list<string>::iterator it = bodies.begin(); it != bodies.end(); it++)
    {
        ephemeris[*it] = ChebyshevEphemeris(spkFile, lskFile, initialEpoch, finalEpoch, *it, center, frame);
    }
}

std::vector<double> EphemerisManager::getState(string body, double epoch)
{
    double state[6] = {0.0};
    ephemeris[body].getState(state, epoch);
    std::vector<double> result = {state[0], state[1], state[2], state[3], state[4], state[5]};
    return result;
}