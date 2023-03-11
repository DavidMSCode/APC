project(APCTests LANGUAGES CXX)

#Include header files
include_directories(${APC_SOURCE_DIR}/include)
include_directories(${APC_SOURCE_DIR}/extern/cspice/include)
#Create binary for the library testing code
add_executable(TestOrbitClass src/OrbitClassTests.cpp)
add_executable(TestEphemerisManager src/EphemerisManagerTests.cpp)
#Link against APC to be able to test code
target_link_libraries(TestOrbitClass  APC)
target_link_libraries(TestEphemerisManager  APC)
#LINK CSPICE AGAINST TEST FUNCTIONS
if (WIN32)
    target_link_libraries(TestEphemerisManager ${APC_SOURCE_DIR}/extern/cspice/lib/cspice.lib) #link cspice library
else ()
    target_link_libraries(TestEphemerisManager ${APC_SOURCE_DIR}/extern/cspice/lib/cspice.a) #link cspice library
endif()
#Add command to be able to run tests
add_test(NAME TestOrbitClass COMMAND TestOrbitClass)
add_test(NAME TestEphemerisManager COMMAND TestEphemerisManager WORKING_DIRECTORY ${APCTests_SOURCE_DIR}/testbin)
#Set CXX standard
set_property(TARGET TestOrbitClass TestEphemerisManager PROPERTY CXX_STANDARD 14)

