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
  py::class_<Orbit>(m, "Orbit")
        .def("getPositionX", &Orbit::getPositionX)
        .def("getPositionY", &Orbit::getPositionY)
        .def("getPositionZ", &Orbit::getPositionZ)
        .def("getVelocityX", &Orbit::getVelocityX)
        .def("getVelocityY", &Orbit::getVelocityY)
        .def("getVelocityZ", &Orbit::getVelocityZ);
}
