project(APCTests LANGUAGES CXX)

#Include header files
include_directories(${APC_SOURCE_DIR}/include)
include_directories(${APC_SOURCE_DIR}/extern/cspice/include)
#Create binary for the library testing code
add_executable(TestOrbitClass src/OrbitClassTests.cpp)
add_executable(TestEphemerisManager src/EphemerisManagerTests.cpp)
add_executable(TestEphemerisRotation src/EphemerisRotationTests.cpp)
find_package(OpenMP REQUIRED)
add_executable(TestMain src/test.cpp) #build a an executable using the test src file
#Link against APC to be able to test code
target_link_libraries(TestOrbitClass  APC)
target_link_libraries(TestEphemerisManager  APC)
target_link_libraries(TestEphemerisRotation  APC)
#LINK CSPICE AGAINST TEST FUNCTIONS
if (WIN32)
    target_link_libraries(TestEphemerisManager ${APC_SOURCE_DIR}/extern/cspice/lib/cspice.lib) #link cspice library
    target_link_libraries(TestEphemerisRotation ${APC_SOURCE_DIR}/extern/cspice/lib/cspice.lib) #link cspice library
else ()
    target_link_libraries(TestEphemerisManager ${APC_SOURCE_DIR}/extern/cspice/lib/cspice.a) #link cspice library
    target_link_libraries(TestEphemerisRotation ${APC_SOURCE_DIR}/extern/cspice/lib/cspice.a) #link cspice library
endif()
#Add command to be able to run tests
add_test(NAME TestOrbitClass COMMAND TestOrbitClass)
add_test(NAME TestEphemerisManager COMMAND TestEphemerisManager WORKING_DIRECTORY ${APCTests_SOURCE_DIR}/testbin)
add_test(NAME TestEphemerisRotation COMMAND TestEphemerisRotation WORKING_DIRECTORY ${APCTests_SOURCE_DIR}/testbin)
add_test(NAME TestMain COMMAND TestMain WORKING_DIRECTORY ${APCTests_SOURCE_DIR}/testbin)
#Set CXX standard
set_property(TARGET  TestOrbitClass TestEphemerisManager TestEphemerisRotation PROPERTY CXX_STANDARD 14)
# set(CMAKE_CXX_STANDARD 14)





target_link_libraries(TestMain APC) #link against APC library
target_link_libraries(TestMain OpenMP::OpenMP_CXX)   #link openmp library