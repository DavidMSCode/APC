include_directories(${APC_SOURCE_DIR}/include ${APC_SOURCE_DIR}/extern/cspice/include ${APC_SOURCE_DIR}/extern/pugixml) #add header file locations for APC, CSPICE, and pugixml
find_package(OpenMP REQUIRED)

# Source files must be listed here to be compiled into the APC library
add_library(APC 
FandG.cpp
EGM2008.cc
APC.cpp
EGM2008.cpp
Orbit.cpp
adaptive_picard_chebyshev.cpp
c_functions.cpp
chebyshev.cpp
Ephemeris.cpp
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
preferences.cpp
picard_chebyshev_propagator.cpp
picard_error_feedback.cpp
picard_iteration.cpp
polydegree_segments.cpp
prepare_propagator.cpp
radial_gravity.cpp
reosc_perigee.cpp
rv2elm.cpp
../extern/pugixml/pugixml.cpp)

if (WIN32)
    target_link_libraries(APC ${APC_SOURCE_DIR}/extern/cspice/lib/cspice.lib) #link cspice library
else ()
    target_link_libraries(APC ${APC_SOURCE_DIR}/extern/cspice/lib/cspice.a) #link cspice library
endif()
target_link_libraries(APC OpenMP::OpenMP_CXX)   #link openmp library