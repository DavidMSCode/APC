include_directories(${APC_SOURCE_DIR}/include ${APC_SOURCE_DIR}/extern/cspice/include)
SET(CMAKE_OSX_ARCHITECTURES x86_64)
add_library(APC 
FandG.cpp
EGM2008.cc
APC.cpp
EGM2008.cpp
Orbit.cpp
adaptive_picard_chebyshev.cpp
c_functions.cpp
chebyshev.cpp
clenshaw_curtis_ivpII.cpp
ecef2eci.cpp
eci2ecef.cpp
interpolate.cpp
linktest.cpp
lsq_chebyshev_fit.cpp
matrix_loader.cpp
perigee_approx.cpp
perturbations.cpp
perturbed_gravity.cpp
picard_chebyshev_propagator.cpp
picard_error_feedback.cpp
picard_iteration.cpp
polydegree_segments.cpp
prepare_propagator.cpp
radial_gravity.cpp
reosc_perigee.cpp
rv2elm.cpp)

target_link_libraries(APC ${APC_SOURCE_DIR}/extern/cspice/lib/cspice.a)