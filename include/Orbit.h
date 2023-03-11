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

#include <string>
#include <vector>

using namespace std;

struct ChebyshevCoefficients{
            std::vector<double> A;          //Position coefficients
            std::vector<double> B;          //Velocity coefficients
            std::vector<double> W1;         //Timescale factors 1 and 2
            std::vector<double> W2; 
            int N;                          //Polydegree
            int coeff_size;                 //length of coefficients
            std::vector<double> seg_times;  //t0 and tf for each segment
            double TF;                      //latest calculated time
            double T0;                      //earliest calculated time
            int total_segs;                 //Total number of segs
        };

class Orbit
{
    private:
        std::vector<std::vector<double> > Soln;
        struct SatProperties{
            double _Area;               
            double _Reflectance;
            double _Mass;
            double _Cd;
        };
            double _mu;                                 //Gravitational parameter of primary body
            string _primary;                            //Primary body
            string _frame;                              //frame for input
            string _epoch = "J2000";                    //epoch for time and date, defaults to J2000
            bool _validOrbit = false;                   //bool to indicate whether the orbit object is valid for APC integration
            bool _validPrimary = false;
            bool _validPosition = false;
            bool _validVelocity = false;
            bool _validTimes = false;
            bool _validDrag = false;
            bool _validSRP = false;
    public:
        bool Compute_Drag;
        bool Compute_SRP;
        bool Compute_Third_Body;
        bool Compute_Hamiltonian;
        bool suborbital;
        int ID;

        std::vector<double> T;              //user defined time vector
        bool USER_TIME = false;             //flag if user time is used
        //Constructors
        Orbit();
        Orbit(string primary, string frame);
        Orbit(std::vector<std::vector<double>> Solution);
        Orbit(double area, double reflectivity, double mass, double Cd, bool compute_drag, bool compute_SRP, bool compute_third_body, bool compute_hamiltonian, int id);
        Orbit(double area, double reflectivity, double mass, double Cd, bool compute_drag, bool compute_SRP, bool compute_third_body, bool compute_hamiltonian, int id, string primary);
        //Setters
        void SetSolution(std::vector<std::vector<double > > Solution);
        void SetProperties(double area, double reflectance, double mass, double Cd, bool compute_drag, bool compute_SRP, bool compute_third_body, bool compute_hamiltonian, int id, string primary = "Earth");
        void SetSubOrbital(); 
        void SetTimeVec(std::vector<double> time_vec);
        void SetCC(std::vector<double> A, std::vector<double> B, std::vector<double> W1, std::vector<double> W2, int N, int coeff_size, std::vector<double> seg_times, double TF, double T0,int total_segs);
        /**
         * @brief Sets gravitational parameter for two body gravity based on the primary body name
         * 
         * @param primary The primary body for integration
         */
        void SetMu(string primary);
        bool InSet(string item, vector<string> validset);
        // Getters
        std::vector<double> getTimes(){return Soln[0];};
        std::vector<double> getPositionX(){return Soln[1];};
        std::vector<double> getPositionY(){return Soln[2];};
        std::vector<double> getPositionZ(){return Soln[3];};
        std::vector<double> getVelocityX(){return Soln[4];};
        std::vector<double> getVelocityY(){return Soln[5];};
        std::vector<double> getVelocityZ(){return Soln[6];};
        std::vector<double> getHamiltonian(){return Soln[7];};
        /**
         * @brief Get the Position of the satellite at the interpolated timesteps
         * 
         * @return std::vector<std::vector<double> > 
         */
        std::vector<std::vector<double> > getPosition();        //rest of method in cpp
        /**
         * @brief Get the Velocities of the satellite at the interpolated timesteps
         * 
         * @return std::vector<std::vector<double> > 
         */
        std::vector<std::vector<double> > getVelocity();        //rest of method in cpp
        /**
         * @brief Get the full interpolated Solution of the propagation
         * 
         * @return std::vector<std::vector<double> > 
         */
        std::vector<std::vector<double> > getSolution(){return Soln;};
        /**
         * @brief Get the Mass of the satellite
         * 
         * @return double 
         */
        double GetMass(){return satproperties._Mass;};
        /**
         * @brief Get the Area of the satellite
         * 
         * @return double 
         */
        double GetArea(){return satproperties._Area;};
        /**
         * @brief Get the Drag Coefficient of the satellite
         * 
         * @return double 
         */
        double GetDragCoefficient(){return satproperties._Cd;};
        /**
         * @brief Get the Reflectance of the satellite
         * 
         * @return double 
         */
        double GetReflectance(){return satproperties._Reflectance;};
        /**
         * @brief Get the Primary Body Gravitational Parameter 
         * 
         * @return double 
         */
        double GetPrimaryGravitationalParameter(){return _mu;};
        /**
         * @brief Get the Primary Body object
         * 
         * @return string 
         */
        string GetPrimaryBody(){return _primary;};
        bool ValidPrimary(){return _validPrimary;};
        //Internal properties struct declaration
        /**
         * @brief Struct containing Mass, Drag Coeficicent, Drag Area and Relflectance of the satellite
         * 
         */
        struct SatProperties satproperties;
        /**
         * @brief Struct storing the chebyshev coefficients for the solution to the propagation
         * 
         */
        struct ChebyshevCoefficients CC;
        
};



#endif
