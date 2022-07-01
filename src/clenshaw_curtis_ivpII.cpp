/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     May 2017
*  LAST MODIFIED:    May 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Generates constant matrices for second order Clensahe & Curtis Quadrature
*
* INPUT:
*    N   -- Chebyshev polynomial order
*    M   -- Number of sample points
*
* OUTPUTS:
*    T2  -- Chebyshev Matrix [(M+1)x(N+1)]
*    P2  -- Picard Iteration Operator [(N+1)xN]
*    T1  -- Chebyshev Matrix [(M+1)xN]
*    P1  -- Picard Integration Operator [Nx(N-1)]
*    Ta  -- Chebyshev Matrix [(M+1)x(N-1)]
*    A   -- Least Squares Operator [(N-1)x(M+1)]
*/



#include <vector>

#include "lsq_chebyshev_fit.h"
#include "c_functions.h"
#include "chebyshev.h"
#include "const.h"


void clenshaw_curtis_ivpII( int N, int M, std::vector<double> &T2, std::vector<double> &P2, std::vector<double> &T1, std::vector<double> &P1, std::vector<double> &Ta, std::vector<double> &A ){

  // Least Squares Operator (A) [(N-1)x(M+1)]
  lsq_chebyshev_fit(-1.0,N-2,M,Ta,A);

  // Compute "Velocity" Constants of Integration (i.e. evaluated T at tau = -1).
  std::vector<double> Lv(N*N,0.0);
  //memset( Lv, 0.0, ((N*N)*sizeof(double)));
  for (int k=0; k<=N-1; k++){
    Lv[ID2(1,k+1,N)] = cos(k*acos(-1));   // Const of Integration (CoI)
  }

  // Compute "Position" Constants of Integration (i.e. evaluated T at tau = -1).
  std::vector<double> Lp((N+1)*(N+1),0.0);
  //memset( Lp, 0.0, ((N+1)*(N+1)*sizeof(double)));
  for (int k=0; k<=N; k++){
    Lp[ID2(1,k+1,N+1)] = cos(k*acos(-1)); // Const of Integration (CoI)
  }

  // S Matrix for Velocity (Sv) [Nx(N-1)]
  std::vector<double> temp1(N*N,0.0);
  //memset( temp1, 0.0, ((N*N)*sizeof(double)));
  temp1[ID2(1,1,N)] = 1.0;
  for (int i=1; i<N; i++){
    for (int j=1; j<N; j++){
      if (i == j){
        temp1[ID2(i+1,j+1,N)] = 1.0/(2.0*i);
      }
    }
  }

  std::vector<double> temp2((N-1)*(N-1),0.0);
  //memset( temp2, 0.0, ((N-1)*(N-1)*sizeof(double)));
  for (int i=1; i<=N-1; i++){
    for (int j=1; j<=N-1; j++){
      if (i == j){
        temp2[ID2(i,j,N-1)] = 1.0;
      }
    }
  }

  std::vector<double> temp3((N-1)*(N-1),0.0);
  //memset( temp3, 0.0, ((N-1)*(N-1)*sizeof(double)));
  for (int i=1; i<=N-1; i++){
    for (int j=3; j<=N-1; j++){
      if (i == j-2){
        temp3[ID2(i,j,N-1)] = -1.0;
      }
    }
  }

  std::vector<double> temp4(N*(N-1),0.0);
  //memset( temp4, 0.0, (N*(N-1)*sizeof(double)));
  for (int i=2; i<=N; i++){
    for (int j=1; j<=N-1; j++){
      temp4[ID2(i,j,N)] = temp2[ID2(i-1,j,N-1)] + temp3[ID2(i-1,j,N-1)];
    }
  }

  std::vector<double> Sv;
  //memset( Sv, 0.0, (N*(N-1)*sizeof(double)));
  Sv = matmul(temp1,temp4,N,N,N-1,N,N);
  Sv[ID2(1,1,N)] = 0.25;
  Sv[ID2(2,1,N)] = 1.0;

  // Picard Integration Operator (Acceleration to Velocity)
  std::vector<double> temp5(N*N,0.0);
  //memset( temp5, 0.0, (N*N*sizeof(double)));
  for (int i=1; i<=N; i++){
    for (int j=1; j<=N; j++){
      temp5[ID2(i,j,N)] = -Lv[ID2(i,j,N)];
      if (i == j){
        temp5[ID2(i,j,N)] = temp5[ID2(i,j,N)] + 1.0;
      }
    }
  }

  P1 = matmul(temp5,Sv,N,N,N-1,N,N);  // [Nx(N-1)]

  std::vector<double> temp6((N+1)*(N+1),0.0);
  //memset( temp6, 0.0, ((N+1)*(N+1)*sizeof(double)));
  temp6[ID2(1,1,N+1)] = 1.0;
  for (int i=1; i<N+1; i++){
    for (int j=1; j<N+1; j++){
      if (i == j){
        temp6[ID2(i+1,j+1,N+1)] = 1.0/(2.0*i);
      }
    }
  }

  std::vector<double> temp7((N)*(N),0.0);
  //memset( temp7, 0.0, ((N)*(N)*sizeof(double)));
  for (int i=1; i<=N; i++){
    for (int j=1; j<=N; j++){
      if (i == j){
        temp7[ID2(i,j,N)] = 1.0;
      }
    }
  }

  std::vector<double> temp8((N)*(N),0.0);
  //memset( temp8, 0.0, ((N)*(N)*sizeof(double)));
  for (int i=1; i<=N; i++){
    for (int j=3; j<=N; j++){
      if (i == j-2){
        temp8[ID2(i,j,N)] = -1.0;
      }
    }
  }

  std::vector<double> temp9((N+1)*N,0.0);
  //memset( temp9, 0.0, ((N+1)*N*sizeof(double)));
  for (int i=2; i<=N+1; i++){
    for (int j=1; j<=N; j++){
      temp9[ID2(i,j,N+1)] = temp7[ID2(i-1,j,N)] + temp8[ID2(i-1,j,N)];
    }
  }

  // S Matrix for Position (Sp) [(N+1)xN]
  std::vector<double> Sp((N+1)*N,0.0);
  //memset( Sp, 0.0, ((N+1)*N*sizeof(double)));
  Sp = matmul(temp6,temp9,N+1,N+1,N,N+1,N+1);
  Sp[ID2(1,1,N+1)] = 0.25;
  Sp[ID2(2,1,N+1)] = 1.0;

  // Picard Integration Operator (Velocity to Position)
  std::vector<double> temp10((N+1)*(N+1),0.0);
  //memset( temp10, 0.0, ((N+1)*(N+1)*sizeof(double)));
  for (int i=1; i<=N+1; i++){
    for (int j=1; j<=N+1; j++){
      temp10[ID2(i,j,N+1)] = -Lp[ID2(i,j,N+1)];
      if (i == j){
        temp10[ID2(i,j,N+1)] = temp10[ID2(i,j,N+1)] + 1.0;
      }
    }
  }

  P2 = matmul(temp10,Sp,N+1,N+1,N,N+1,N+1);  // [(N+1)xN]

  // Chebyshev Matrix (interpolate "velocity" coefficients)
  chebyshev(-1.0,N-1,M,2,T1);    // arg4 = 2 -> Trig Cheby Poly

  // // Chebyshev Matrix (interpolate "position" coefficients)
  chebyshev(-1.0,N,M,2,T2);    // arg4 = 2 -> Trig Cheby Poly

}
