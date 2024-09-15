include_directories(${APC_SOURCE_DIR}/include ${APC_SOURCE_DIR}/extern/cspice/include) #add header file locations for APC and CSPICE
find_package(OpenMP REQUIRED)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
 

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
Constellation.cpp
Ephemeris.cpp
EphemerisRotation.cpp
clenshaw_curtis_ivpII.cpp
ecef2eci.cpp
eci2ecef.cpp
flags.cpp
Gravity.cpp
GRGM1200b.cc
GRGM1200b.cpp
interpolate.cpp
lsq_chebyshev_fit.cpp
lunar_perturbed_gravity.cpp
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
rv2elm.cpp
TwoBody.cpp) 

add_library(Grav
EGM2008.cc
EGM2008.cpp
GRGM1200b.cpp
GRGM1200b.cc
radial_gravity.cpp)
# use c++14
set_property(TARGET Grav PROPERTY CXX_STANDARD 14)


if (WIN32)
    target_link_libraries(APC ${APC_SOURCE_DIR}/extern/cspice/lib/cspice.lib) #link cspice library
else ()
    target_link_libraries(APC ${APC_SOURCE_DIR}/extern/cspice/lib/cspice.a) #link cspice library
endif()
target_link_libraries(APC OpenMP::OpenMP_CXX)   #link openmp library
# target_link_libraries (APC Eigen3::Eigen)

add_executable(main_bootstrap main_bootstrap.cpp)
target_link_libraries(main_bootstrap APC)

add_executable(main_degree_trade_study main_degree_trade_study.cpp TwoBody.cpp)
target_link_libraries(main_degree_trade_study APC)
set_property(TARGET main_degree_trade_study PROPERTY CXX_STANDARD 14)
target_link_libraries(main_degree_trade_study OpenMP::OpenMP_CXX)   #link openmp library

add_executable(main_hamiltonian_error main_hamiltonian_error.cpp TwoBody.cpp)
target_link_libraries(main_hamiltonian_error APC)
set_property(TARGET main_hamiltonian_error PROPERTY CXX_STANDARD 14)

add_executable(RK45 RK45.cpp)
target_link_libraries(RK45 APC)
set_property(TARGET RK45 PROPERTY CXX_STANDARD 14)

add_executable(main_bootstrap_tradestudy main_bootstrap_tradestudy.cpp)
target_link_libraries(main_bootstrap_tradestudy APC)
set_property(TARGET main_bootstrap_tradestudy PROPERTY CXX_STANDARD 14)

add_executable(main_interpolate main_interpolate.cpp)
target_link_libraries(main_interpolate APC)

add_executable(main_interpolate_tradestudy main_interpolate_tradestudy.cpp)
target_link_libraries(main_interpolate_tradestudy APC)

add_executable(CCQ_tests CCQ_tests.cpp)
target_link_libraries(CCQ_tests APC)
set_property(TARGET CCQ_tests PROPERTY CXX_STANDARD 14)