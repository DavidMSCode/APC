//This code uses Runge Kutta 45 method to integrate the orbit of a satellite around the moon using a high fidelity lunar gravity model.
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <queue>
#include "lunar_perturbed_gravity.h"
#include "TwoBody.h"
#include "const.h"
#include "GRGM1200b.h"
#include "EphemerisRotation.h"
#include "APC.h"

using namespace std;

vector<double> mult(double a, vector<double> b)
{
    // This function multiplies a scalar with a vector and returns the product as a vector.
    vector<double> c;
    for (int i = 0; i < b.size(); i++)
    {
        c.push_back(a * b[i]);
    }
    return c;
}

vector<double> add(vector<double> a, vector<double> b)
{
    // This function adds two vectors and returns the sum as a vector.
    vector<double> c;
    for (int i = 0; i < a.size(); i++)
    {
        c.push_back(a[i] + b[i]);
    }
    return c;
}

vector<double> polyinterp(deque<vector<double>> yq, deque<double> tq, double t)
{
    // This func uses Neville's algorithm to interpolate the state vector at time t given the state vector and time dequeues of the same length
    int n = yq.size();
    double Q[n][n];
    vector<double> y = yq[0];
    // iterate through each state variable
    for (int ll = 0; ll < yq[0].size(); ll++)
    {
        // zero the Q array
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                Q[i][j] = 0;
            }
        }
        // set the first column of Q to the llth coordinate values
        for (int i = 0; i < n; i++)
        {
            Q[i][0] = yq[i][ll];
        }

        // Neville's Method
        for (int i = 1; i < n; i++)
        {
            for (int j = 1; j <= i; j++)
            {
                Q[i][j] = ((t - tq[i - j]) * Q[i][j - 1] - (t - tq[i]) * Q[i - 1][j - 1]) / (tq[i] - tq[i - j]);
            }
        }
        // Store the last element of the last row of Q in the state vector
        y[ll] = Q[n - 1][n - 1];
    }
    // return solution
    return y;
}

class integrator
{
public:
    pair<vector<vector<double>>, vector<double>> integrate(double t0, double tf, vector<double> y0, double dt = 30);                                                                                 // The RK45 integrator
    pair<vector<vector<double>>, vector<double>> integrate(double t0, double tf, vector<double> y0, function<vector<double>(double t, vector<double>)> EOMFunc, double dt = 30, double epsi = 1e-9); // The RK45 integrator with arbitrary EOM function
    virtual vector<double> f(double t, vector<double> y) = 0;                                                                                                                                        // ODE equations of motion to be replaced in child class
    double TE(vector<double> k0, vector<double> k1, vector<double> k2, vector<double> k3, vector<double> k4, vector<double> k5);                                                                     // The estimated truncation error function
    // RK45 coefficients
    const double A[6] = {0, 1. / 4, 3. / 8, 12. / 13, 1, 1. / 2};                                                                                                                                                                                          // The time steps
    const double B[6][5] = {{0, 0, 0, 0, 0}, {1. / 4, 0, 0, 0, 0}, {3. / 32, 9. / 32, 0, 0, 0}, {1932. / 2197, -7200. / 2197, 7296. / 2197, 0, 0}, {439. / 216, -8, 3680. / 513, -845. / 4104, 0}, {-8. / 27, 2, -3544. / 2565, 1859. / 4104, -11. / 40}}; // The Butcher tableau for alpha2 = 3/8
    const double CH[6] = {16. / 135, 0, 6656. / 12825, 28561. / 56430, -9. / 50, 2. / 55};                                                                                                                                                                 // The weights for the 5th order solution
    const double CT[6] = {-1. / 360, 0, 128. / 4275, 2197. / 75240, -1. / 50, -2. / 55};                                                                                                                                                                   // The weights for the truncation error
};
double integrator::TE(vector<double> k0, vector<double> k1, vector<double> k2, vector<double> k3, vector<double> k4, vector<double> k5)
{
    // This function takes the 6 k vectors and returns the truncation error.
    vector<double> kt = add(mult(CT[0], k0), add(mult(CT[1], k1), add(mult(CT[2], k2), add(mult(CT[3], k3), add(mult(CT[4], k4), mult(CT[5], k5))))));
    double TE = 0;
    for (int i = 0; i < 6; i++)
    {
        if (i < 3)
        {
            TE += pow(kt[i], 2) / DU;
        }
        else
        {
            TE += pow(kt[i], 2) / VU;
        }
    }
    TE = sqrt(TE);
    return TE;
}

pair<vector<vector<double>>, vector<double>> integrator::integrate(double t0, double tf, vector<double> y0, double dt)
{
    double epsi = 1e-15;
    // This function uses the runge kutta 45 method to find the final state vector given the initial time, final time and initial state vector. And stores the solution every dt seconds.
    pair<vector<vector<double>>, vector<double>> sol = integrate(
        t0, tf, y0, [this](double t, vector<double> y)
        {
            return f(t, y); // The class built in equations of motion
        },
        dt, epsi);
    return sol; // Return the solution
};

pair<vector<vector<double>>, vector<double>> integrator::integrate(double t0, double tf, vector<double> yinit, function<vector<double>(double, vector<double>)> f, double dt, double epsi)
{
    // This function uses the runge kutta 45 method to find the final state vector given the initial time, final time and initial state vector. And stores the solution every dt seconds.
    vector<double> ts = {t0};
    vector<vector<double>> ys = {yinit};
    double next_dt = dt;                   // The next time at which the solution is to be stored
    double h = 0.1;                        // The initial time step
    double t = t0;                         // The current time
    double te;                             // The truncation error
    vector<double> y = yinit;              // The current state vector
    vector<double> k0, k1, k2, k3, k4, k5; // The RK45 estimates
    vector<double> y0, y1, y2, y3, y4, y5; // The intermediate state vectors
    int N = 6;                             // number of sample points for output vector

    // This dequeue stores the last 4 state vectors
    deque<vector<double>> yq;
    // This dequeue stores the last 4 time vectors
    deque<double> tq;
    for (int i = 0; i < N; i++)
    {
        yq.push_back(yinit);
        tq.push_back(t0);
    }

    // Counts the number of timesteps past the desired solution time
    int steps_past = 0;
    bool stored_final = false;
    while (not(t > tf && stored_final))
    {
        y0 = y;
        k0 = mult(h, f(t + A[0] * h, y0));
        y1 = add(y, mult(B[1][0], k0));
        k1 = mult(h, f(t + A[1] * h, y1));
        y2 = add(y, add(mult(B[2][0], k0), mult(B[2][1], k1)));
        k2 = mult(h, f(t + A[2] * h, y2));
        y3 = add(y, add(mult(B[3][0], k0), add(mult(B[3][1], k1), mult(B[3][2], k2))));
        k3 = mult(h, f(t + A[3] * h, y3));
        y4 = add(y, add(mult(B[4][0], k0), add(mult(B[4][1], k1), add(mult(B[4][2], k2), mult(B[4][3], k3)))));
        k4 = mult(h, f(t + A[4] * h, y4));
        y5 = add(y, add(mult(B[5][0], k0), add(mult(B[5][1], k1), add(mult(B[5][2], k2), add(mult(B[5][3], k3), mult(B[5][4], k4))))));
        k5 = mult(h, f(t + A[5] * h, y5));
        vector<double> yh = add(y, add(mult(CH[0], k0), add(mult(CH[1], k1), add(mult(CH[2], k2), add(mult(CH[3], k3), add(mult(CH[4], k4), mult(CH[5], k5))))))); // The 5th order solution
        te = TE(k0, k1, k2, k3, k4, k5);                                                                                                                           // The truncation error
        if (te > epsi)
        {
            h = 0.9 * h * pow(epsi / te, 1. / 5); // The new time step
        }
        else
        {
            t = t + h; // The new time
            y = yh;    // The new state vector
            // add new values to dequeus and remove the first elements
            yq.push_back(y);
            yq.pop_front();
            tq.push_back(t);
            tq.pop_front();

            if (t > next_dt)
            {
                // Interpolate solution if N/2 steps past desired solution time ()
                steps_past += 1; // increment the number of steps past the desired solution time
                if (steps_past == N / 2)
                {
                    // interpolate the output vector
                    vector<double> y_interp = polyinterp(yq, tq, next_dt);
                    steps_past = 0;         // reset steps past desired solution counter
                    ts.push_back(next_dt);  // Store the time
                    ys.push_back(y_interp); // Store the interpolated state vector

                    if (next_dt == tf)
                    {
                        stored_final = true; // check if the final solution has been stored
                    }
                    else
                    {
                        next_dt += dt; // Update the next time at which the solution is to be stored
                        if (next_dt > tf)
                            next_dt = tf; // make sure the last solution is stored at tf
                    }
                }
            }
            h = 0.9 * h * pow(epsi / te, 1. / 5); // The new time step
        }
    };
    pair<vector<vector<double>>, vector<double>> sol(ys, ts); // The solution
    return sol;                                               // Return the solution
};

// Earth Orbit is a child class of integrator class.
class EarthOrbit : public integrator
{
public:
    vector<double> f(double t, vector<double> y);
};

vector<double> EarthOrbit::f(double t, vector<double> y)
{
    // This function takes the 6 dimension state y (km,km,km,km/s,km/s,km/s) and returns the change in position and velocity due to Earth's gravity in the inertial frame.
    double mu = 398600.4418;
    double r = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
    vector<double> f = {y[3], y[4], y[5], -mu * y[0] / (r * r * r), -mu * y[1] / (r * r * r), -mu * y[2] / (r * r * r)};
    return f;
}

// J2EarthOrbit is a child class of integrator class.
class J2EarthOrbit : public integrator
{
public:
    vector<double> f(double t, vector<double> y);
};
vector<double> J2EarthOrbit::f(double t, vector<double> y)
{
    // This function takes the 6 dimension state y (km,km,km,km/s,km/s,km/s) and returns the change in position and velocity due to Earth's gravity with J2 perturbation in the inertial frame.
    double mu = 398600.4418;
    double r = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
    double J2 = 0.0010826269;
    double Re = 6378.137;
    vector<double> f = {y[3], y[4], y[5], -mu * y[0] / (r * r * r) + 1.5 * mu * J2 * Re * Re / r / r / r / r * (y[0] / r * (5 * y[2] * y[2] / r / r - 1)), -mu * y[1] / (r * r * r) + 1.5 * mu * J2 * Re * Re / r / r / r / r * (y[1] / r * (5 * y[2] * y[2] / r / r - 1)), -mu * y[2] / (r * r * r) + 1.5 * mu * J2 * Re * Re / r / r / r / r * (y[2] / r * (5 * y[2] * y[2] / r / r - 3))};
    return f;
}

class HOLunarOrbit : public integrator
{
public:
    double Feval[2]={0.0};
    vector<double> f(double t, vector<double> y);
};
double gravdegree=100;
vector<double> HOLunarOrbit::f(double t, vector<double> y)
{
    // double Feval[2] = {0.0};
    double acc[3] = {0.0};
    double tol = 1e-15;
    double deg = gravdegree;
    vector<double> XI = {y[0], y[1], y[2]};
    // Get acceleration in moon PA frame
    // vector<double> FState = InertialToMoonPA(y, t);
    vector<double> FState = y;
    double *XF = &FState[0];
    lunar_Grav_Full(t, XF, acc, tol, deg, Feval);
    // convert to inertial frame
    // vector<double> aI = MoonPAAccelerationToInertial(acc, t);
    std::vector<double> aI = {acc[0], acc[1], acc[2]};
    return {y[3], y[4], y[5], aI[0], aI[1], aI[2]};
}
void jacobiIntegral_GRGM1200b(double t, vector<double> SolN, double *H, int Deg)
{
    double xI[3] = {0.0};
    double vI[3] = {0.0};
    double vF[3] = {0.0};

    xI[0] = SolN[0];
    xI[1] = SolN[1];
    xI[2] = SolN[2];
    vI[0] = SolN[3];
    vI[1] = SolN[4];
    vI[2] = SolN[5];
    // Convert from inertial frame to body fixed frame
    // vector<double> stateF = InertialToMoonPA(SolN, t);
    vector<double> stateF = SolN;
    double *xF = &stateF[0];
    double KE, PE;

    KE = 0.5 * (vI[0] * vI[0] + vI[1] * vI[1] + vI[2] * vI[2]); // Kinetic energy in inertial frame
    GRGM1200bPot(xF, &PE, Deg);                                 // Potential energy from body fixed position
    PE = -PE;
    *H = PE + KE; // Hamiltonian

    // printf("KE: %e\tPE: %e\tRT: %e\tSum: %e\n ",KE,PE,RotTerm,*H);
    // getchar();
}

vector<double> calcHamiltonians(pair<vector<vector<double>>, vector<double>> sol)
{
    vector<vector<double>> ys = sol.first;
    vector<double> ts = sol.second;
    // preallocate Hs
    vector<double> Hs;
    Hs.reserve(ts.size());
    double H0 = 0.0;
    // iterate through each state vector
    for (int i = 0; i < ys.size(); i++)
    {
        double t = ts[i];
        vector<double> SolN = ys[i];
        // calculate the hamiltonian
        double H = 0.0;
        jacobiIntegral_GRGM1200b(t, SolN, &H, 100);
        if (i == 0)
        {
            H0 = H;
        }
        // store the hamiltonian
        Hs.push_back(fabs((H - H0) / H0));
    }
    return Hs;
}

int main()
{
    double alt = 100; // km
    double a = C_Rmoon + alt;
    double e = 0.0;
    double i = 0.0;
    double raan = 0.0;
    double aop = 0.0;
    double ta = 0.0;
    vector<vector<double>> states = elms2rv(a, e, i, raan, aop, ta, C_MU_MOON);
    // integrates the orbit of a LEO satellite with a 45 deg inclination for 1 day.
    vector<double> r0 = states[0];
    vector<double> v0 = states[1];
    vector<double> y0 = {r0[0], r0[1], r0[2], v0[0], v0[1], v0[2]};
    double t0 = 0;
    double T = 2*C_PI*sqrt(pow(a,3)/C_MU_MOON);                             //Orbital period (s)
    double tf = t0+T;
    // Integrate the orbit without any perturbations
    LoadKernels();
    HOLunarOrbit HO;
    pair<vector<vector<double>>, vector<double>> sol = HO.integrate(t0, tf, y0);
    vector<double> Hs = calcHamiltonians(sol);
    cout<<"Hmax:  "<<*max_element(Hs.begin(), Hs.end())<<endl;
    cout<<"Total Full Gravity Evaluations: "<<HO.Feval[0]<<endl;
    UnloadKernels();
    return 0;
}