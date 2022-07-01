/*
*  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
*  DATE WRITTEN:     May 2017
*  LAST MODIFIED:    May 2017
*  AFFILIATION:      Department of Aerospace Engineering, Texas A&M University, College Station, TX
*  DESCRIPTION:      Interpolate solution onto user specified times
*
* INPUT:
*    ALPHA      -- Position Coefficients
*    BETA       -- Velocity Coefficients
*    soln_size  -- Length of solution
*    coeff_size -- Length of Coefficients
*    N          -- Degree of polynomial
*    seg_times  -- t0 and tf for each segment
*    t0         -- Initial time (s)
*    tf         -- Final time (s)
*    dt         -- User desired output time interval
*    total_segs -- Total number of segments
*
* OUTPUTS:
*    Soln       -- Output Solution (Position (km) & Velocity (km/s))
*
* COMMENTS:
*
*/
#include <vector>
#include <iostream>

#include "interpolate.h"
#include "c_functions.h"


std::vector<double> interpolate(std::vector<double>  ALPHA, std::vector<double> BETA, int soln_size, int coeff_size, int N, std::vector<double> seg_times,
  std::vector<double> W1, std::vector<double> W2, double t0, double tf, double dt, int total_segs){
  
  int prev_cnt = 0;
  std::vector<double> Soln(soln_size*6,0.0);
  // User specified output times
  int len;
  len = int(ceil(tf/dt));
  std::vector<double> time_out(len+1,0.0);
  //memset( time_out, 0.0, (len*sizeof(double)));
  time_out[0] = t0;
  for (int ii=1; ii<=len; ii++){
    time_out[ii] = time_out[ii-1] + dt;
  }
  double test_time = 0.0;
  // Loop through all segments
  for (int i=1; i<=total_segs; i++){
    int sz = int(ceil((seg_times[i]-seg_times[i-1])/dt));
    // Initialization
    std::vector<double> Beta(N*3);
    std::vector<double> Alpha((N+1)*3);
    std::vector<double> tt(sz);
    std::vector<double> tau(sz);

    double w1, w2;
    w1 = W1[i-1];
    w2 = W2[i-1];

    // User desired times for a given segment
    int cnt = 0;
    // printf("seg_times %f\t%f\n",seg_times[i-1],seg_times[i]);
    for (int j=0; j<=len; j++){
      if (time_out[j] == seg_times[i-1]){
        //printf("t1 %f\n",time_out[j]);
        tau[cnt] = -1.0;
        cnt      = cnt + 1;   // Number of times steps per segment
      }
      if (time_out[j] == seg_times[i]){
        // printf("t2 %f\n",time_out[j]);
        tau[cnt] = 1.0;
        cnt = cnt + 1;
      }
      if (time_out[j] > seg_times[i-1] && time_out[j] < seg_times[i]){
        tt[cnt]  = time_out[j];
        tau[cnt] = (tt[cnt] - w1)/w2;
        cnt      = cnt + 1;   // Number of times steps per segment
      }
    }
  
    // Chebyshev Velocity & Position Matrices
    std::vector<double> Tv(cnt*N);
    std::vector<double> Tp(cnt*(N+1));

    for (int t=1; t<=cnt; t++){
      for (int kk=0; kk<=N-1; kk++){
        // Velocity
        Tv[ID2(t,kk+1,cnt)] = cos(kk*acos(tau[t-1]));
      }
      for (int kk=0; kk<=N; kk++){
        // Position
        Tp[ID2(t,kk+1,cnt)] = cos(kk*acos(tau[t-1]));
      }
    }

    // Velocity Coefficients for a Segment
    for (int p=1; p<=N; p++){
      Beta[ID2(p,1,N)] = BETA[ID2(p+((i-1)*N),1,coeff_size)];
      Beta[ID2(p,2,N)] = BETA[ID2(p+((i-1)*N),2,coeff_size)];
      Beta[ID2(p,3,N)] = BETA[ID2(p+((i-1)*N),3,coeff_size)];
    }
    std::vector<double> v_interp;
    v_interp = matmul(Tv,Beta,cnt,N,3,cnt,N);

    //sanity check
    int check = ID2(cnt+prev_cnt,6,soln_size);
    if (check >=soln_size*6){
        std::cout << "exceeding soln size\n";
    }

    // Velocity
    for (int p=1; p<=cnt; p++){
      Soln[ID2(p+prev_cnt,4,soln_size)] = v_interp[ID2(p,1,cnt)];
      Soln[ID2(p+prev_cnt,5,soln_size)] = v_interp[ID2(p,2,cnt)];
      Soln[ID2(p+prev_cnt,6,soln_size)] = v_interp[ID2(p,3,cnt)];
    }

    // Position Coefficients for a Segment
    for (int p=1; p<=N+1; p++){
      Alpha[ID2(p,1,N+1)] = ALPHA[ID2(p+((i-1)*(N+1)),1,coeff_size)];
      Alpha[ID2(p,2,N+1)] = ALPHA[ID2(p+((i-1)*(N+1)),2,coeff_size)];
      Alpha[ID2(p,3,N+1)] = ALPHA[ID2(p+((i-1)*(N+1)),3,coeff_size)];
    }
    std::vector<double> x_interp;
    x_interp = matmul(Tp,Alpha,cnt,N+1,3,cnt,N+1);
    // Position
    for (int p=1; p<=cnt; p++){
      Soln[ID2(p+prev_cnt,1,soln_size)] = x_interp[ID2(p,1,cnt)];
      Soln[ID2(p+prev_cnt,2,soln_size)] = x_interp[ID2(p,2,cnt)];
      Soln[ID2(p+prev_cnt,3,soln_size)] = x_interp[ID2(p,3,cnt)];
      // printf("Soln %f\t%f\t%f\n",Soln[ID2(p+prev_cnt,1,soln_size)],Soln[ID2(p+prev_cnt,2,soln_size)],Soln[ID2(p+prev_cnt,3,soln_size)]);
    }
    // printf("tfdt %f\t%f\n",tf,dt);
    prev_cnt = prev_cnt + cnt;  // Counter to track position in Soln array
  }return Soln;
}
