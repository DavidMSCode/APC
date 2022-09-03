/*
*  AUTHORS:          David Stanley (davidms4@illinois.edu)
*  DATE WRITTEN:     Feb 2022
*  LAST MODIFIED:    Feb 2022
*  AFFILIATION:      Department of Aerospace Engineering, University of Illinois Urbana-Champaign
*  DESCRIPTION:      Orbit class for storing orbit properties and solution
*  REFERENCE:        Woollands, R., and Junkins, J., "Nonlinear Differential Equation Solvers
*                    via Adaptive Picard-Chebyshev Iteration: Applications in Astrodynamics", JGCD, 2016.
*/


#ifndef ORBIT_H
#define ORBIT_H

#include <vector>
struct ChebyshevCoefficients{
            std::vector<double> A;          //Position coefficients
            std::vector<double> B;          //Velocity coefficients
            std::vector<double> W1;         //Timescale factors 1 and 2
            std::vector<double> W2; 
            int N;                          //Polydegree
            int coeff_size;                 //length of coefficients
            std::vector<double> seg_times;  //t0 and tf for each segment
            double TF;                      //latest calculated time
            double T0;                      //earliest calcualted time
            int total_segs;                 //Total number of segs
        };

class Orbit
{
    private:
        std::vector<std::vector<double> > Soln;
        struct SatProperties{
            double Area;
            double Reflectance;
            double Mass;
            double Cd;
        };

    public:
        bool Compute_Drag;
        bool Compute_SRP;
        bool Compute_Third_Body;
        bool suborbital;
        int ID;
        std::vector<double> T;              //user defined time vector
        bool USER_TIME = false;             //flag if user time is used
        //Constructors
        Orbit();
        Orbit(std::vector<std::vector<double > > Solution);
        Orbit(double area, double reflectivity, double mass, double Cd, bool compute_drag, bool compute_SRP, bool compute_third_body, int id);
        //Setters
        void SetSolution(std::vector<std::vector<double > > Solution);
        void SetProperties(double area, double reflectance, double mass, double Cd, bool compute_drag, bool compute_SRP, bool compute_third_body, int id);
        void SetSubOrbital(); 
        void SetTimeVec(std::vector<double> time_vec);
        void SetCC(std::vector<double> A, std::vector<double> B, std::vector<double> W1, std::vector<double> W2, int N, int coeff_size, std::vector<double> seg_times, double TF, double T0,int total_segs);
        //Getters
        std::vector<double> getTimes(){return Soln[0];};
        std::vector<double> getPositionX(){return Soln[1];};
        std::vector<double> getPositionY(){return Soln[2];};
        std::vector<double> getPositionZ(){return Soln[3];};
        std::vector<double> getVelocityX(){return Soln[4];};
        std::vector<double> getVelocityY(){return Soln[5];};
        std::vector<double> getVelocityZ(){return Soln[6];};
        std::vector<double> getHamiltonian(){return Soln[7];};
        std::vector<std::vector<double> > getPosition();        //rest of method in cpp
        std::vector<std::vector<double> > getVelocity();        //rest of method in cpp
        std::vector<std::vector<double> > getSolution(){return Soln;};
        double GetMass(){return satproperties.Mass;};
        double GetArea(){return satproperties.Area;};
        double GetDragCoefficient(){return satproperties.Cd;};
        double GetReflectance(){return satproperties.Reflectance;};
        //Internal properties struct declaration
        struct SatProperties satproperties;
        struct ChebyshevCoefficients CC;
        
};



#endif
