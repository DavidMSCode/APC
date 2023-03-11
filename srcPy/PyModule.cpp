/*
*  AUTHORS:          David Stanley (davidms4@illinois.edu)
*  DATE WRITTEN:     Feb 2022
*  LAST MODIFIED:    Feb 2022
*  AFFILIATION:      Department of Aerospace Engineering, University of Illinois Urbana-Champaign
*  DESCRIPTION:      Pybind11 code for making methods that are acessible from python.
*  REFERENCE:        Woollands, R., and Junkins, J., "Nonlinear Differential Equation Solvers
*                    via Adaptive Picard-Chebyshev Iteration: Applications in Astrodynamics", JGCD, 2016.
*/

#include <pybind11/pybind11.h>      //Pybind11 headers must be first include
#include <pybind11/stl.h>
#include <Orbit.h>
#include <linktest.h>
#include <APC.h>
#include "matrix_loader.h"


using namespace pybind11::literals;
namespace py = pybind11;


PYBIND11_MODULE(APC, m) {
      m.doc() = "Test plugin for adaptive picard chebychev integrator";
      using namespace py::literals;
      m.def("PropagateICs", &PropagateICs, "takes satellite state around Earth and returns the orbital state vectors for a given time interval");
      m.def("PropagateOrbit", &PropagateOrbit, "Takes initial conditions and returns an orbit class object");
      m.def("Linktest",&Linktest,"testing linking against cspice");
      m.def("GenSigma13",&GenSigma13,"Generate 13 perturbed variations of an orbit");
      m.def("ParallelPropagate",&ParallelPropagate,"Propagates multiple orbits in parallel");
      m.def("SinglePropagate",static_cast<class Orbit (*)(std::vector<double>, std::vector<double>, double, double, double, double, double, double, bool, bool, bool, bool )>(&SinglePropagate),"Propagates a single orbit","r0"_a,"v0"_a,"t0"_a,"tf"_a,"area"_a,"reflectance"_a,"mass"_a,"Cd"_a,"compute_drag"_a=false,"compute_srp"_a=false,"compute_third_body"_a=false, "compute_hamiltonian"_a=false);
      m.def("SinglePropagate",static_cast<class Orbit (*)(std::vector<double>, std::vector<double>, std::vector<double>, double, double, double, double, bool, bool, bool, bool) >(&SinglePropagate),"Propagates a single orbit with user specified time vector","r0"_a,"v0"_a,"time_vec"_a,"area"_a,"reflectance"_a,"mass"_a,"Cd"_a,"compute_drag"_a=false,"compute_srp"_a=false,"compute_third_body"_a=false, "compute_hamiltonian"_a=false);
      m.def("Benchmark1000",&Benchmark1000,"Returns the number of threads used and the time (s) to complete 1000 orbit propagations");
      py::class_<Orbit>(m, "Orbit")
            .def("getTimes", &Orbit::getTimes)
            .def("getPositionX", &Orbit::getPositionX)
            .def("getPositionY", &Orbit::getPositionY)
            .def("getPositionZ", &Orbit::getPositionZ)
            .def("getVelocityX", &Orbit::getVelocityX)
            .def("getVelocityY", &Orbit::getVelocityY)
            .def("getVelocityZ", &Orbit::getVelocityZ)
            .def("getPosition", &Orbit::getPosition)
            .def("getVelocity", &Orbit::getVelocity)
            .def("getHamiltonian", &Orbit::getHamiltonian)
            .def_readwrite("CC",&Orbit::CC)
            .def_readwrite("T",&Orbit::T);
      py::class_<ChebyshevCoefficients>(m, "ChebyshevCoefficients")
            .def(py::init<>())
            .def_readwrite("A",&ChebyshevCoefficients::A)
            .def_readwrite("B",&ChebyshevCoefficients::B)
            .def_readwrite("W1",&ChebyshevCoefficients::W1)
            .def_readwrite("W2",&ChebyshevCoefficients::W2)
            .def_readwrite("N",&ChebyshevCoefficients::N)
            .def_readwrite("coeff_size",&ChebyshevCoefficients::coeff_size)
            .def_readwrite("seg_times",&ChebyshevCoefficients::seg_times)
            .def_readwrite("TF",&ChebyshevCoefficients::TF)
            .def_readwrite("T0",&ChebyshevCoefficients::T0)
            .def_readwrite("total_segs",&ChebyshevCoefficients::total_segs);
      py::class_<SatState>(m, "SatState")
            .def(py::init<>())
            .def_readwrite("r", &SatState::r)
            .def_readwrite("v", &SatState::v);
      //Functions for debugging    
      m.def("matrix_loader",&matrix_loader,"Loads Picard iteration matrices into memory.");
      m.def("MatricesLoaded",&MatricesLoaded,"Returns true if Picard iteration matrices have been loaded.");
}
