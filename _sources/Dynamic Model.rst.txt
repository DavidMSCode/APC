Dynamic Model
=============
APC's dynamic model includes a variable fidelity gravity model. The gravity perturbations are calculated using the EGM2008 spherical harmonic gravity model up to a maximum degree of 100. This model is used to demonstrate the propagator's efficiency with high order gravity models. The model also includes simple cannonball drag model for LEO orbits as well as a solar radiation pressure perturbation and solilunar third body gravity perturbations for MEO & GEO orbits. These models are not necessary for APC to function and may be easily replaced by substituting their functions with alternatives or editing the picard_iteration.cpp file. Though to make full use of APC the dynamic model should be capable of calculating an approximation of the gravitational model or dynamic model during iterations where the full calculation is unecessary (see perturbed_gravity.cpp).

Perturbed Gravity Functions
===========================

.. doxygenfile:: perturbed_gravity.h
   :project: APC

Other perturbations
===================

.. doxygenfile:: perturbations.h
   :project: APC
