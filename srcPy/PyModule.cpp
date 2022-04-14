/*
*  AUTHORS:          David Stanley (davidms4@illinois.edu)
*  DATE WRITTEN:     Feb 2022
*  LAST MODIFIED:    Feb 2022
*  AFFILIATION:      Department of Aerospace Engineering, University of Illinois Urbana-Champaign
*  DESCRIPTION:      Pybind11 code for making methods that are acessible from python.
*  REFERENCE:        Woollands, R., and Junkins, J., "Nonlinear Differential Equation Solvers
*                    via Adaptive Picard-Chebyshev Iteration: Applications in Astrodynamics", JGCD, 2016.
*/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <Orbit.h>
#include <linktest.h>
#include <APC.h>

namespace py = pybind11;

PYBIND11_MODULE(APC, m) {
  m.doc() = "Test plugin for adaptive picard chebychev integrator";
  using namespace py::literals;
  m.def("PropagateICs", &PropagateICs, "takes satellite state around Earth and returns the orbital state vectors for a given time interval");
  m.def("PropagateOrbit", &PropagateOrbit, "Takes initial conditions and returns an orbit class object");
  m.def("Linktest",&Linktest,"testing linking against cspice");
  m.def("GenSigma13",&GenSigma13,"Generate 13 perturbed variations of an orbit");
  m.def("ParallelPropagate",&ParallelPropagate,"Propagates multiple orbits in parallel");
  m.def("SinglePropagate",&SinglePropagate,"Propagates a single orbit");
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
        .def("getHamiltonian", &Orbit::getHamiltonian);
  py::class_<SatState>(m, "SatState")
        .def_readwrite("r", &SatState::r)
        .def_readwrite("v", &SatState::v);
}
