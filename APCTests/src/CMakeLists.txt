#Include header files
include_directories(${APC_SOURCE_DIR}/include)
include_directories(${APC_SOURCE_DIR}/extern/cspice/include)


#Add all tests to single test executable
create_test_sourcelist(Tests 
CombinedTests.cpp
OrbitClassTests.cpp
EphemerisManagerTests.cpp
EphemerisRotationTests.cpp
ConstellationTests.cpp
GravTests.cpp
)

add_executable(CombinedTests ${Tests})
# remove test driver from test executable
set(TestsToRun ${Tests})
remove(TestsToRun CombinedTests.cpp)

find_package(OpenMP REQUIRED)
#Link against APC to be able to test code
target_link_libraries(CombinedTests APC) #link against APC library
target_link_libraries(CombinedTests OpenMP::OpenMP_CXX)   #link openmp library
target_link_libraries(CombinedTests Grav) #link gtest library
#LINK CSPICE AGAINST TEST FUNCTIONS
if (WIN32)
    target_link_libraries(CombinedTests ${APC_SOURCE_DIR}/extern/cspice/lib/cspice.lib) #link cspice library
else ()
    target_link_libraries(CombinedTests ${APC_SOURCE_DIR}/extern/cspice/lib/cspice.a) #link cspice library
endif()


# Add all tests to ctest
foreach (test ${TestsToRun})
  get_filename_component (TName ${test} NAME_WE)
  add_test (NAME ${TName} COMMAND CombinedTests ${TName})
endforeach ()


#copy necessary files to build directory
# Set the resource directory
set(RESOURCE_DIR ${CMAKE_SOURCE_DIR}/res)

# Set the build directory
set(BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(copy_assets COMMAND ${CMAKE_COMMAND} -E copy_directory ${RESOURCE_DIR} ${BUILD_DIR})
add_dependencies(CombinedTests copy_assets)
# # Copy all files and subdirectories from the res directory to the build directory on build if changes have been made to the res directory
# file(GLOB RESOURCES ${RESOURCE_DIR}/*)
# foreach(RESOURCE ${RESOURCES})
#     add_custom_command(
#         TARGET CombinedTests POST_BUILD
#         COMMAND ${CMAKE_COMMAND} -E copy_directory
#                 ${RESOURCE}
#                 ${BUILD_DIR}
#         COMMENT "Copying resource ${RESOURCE}"
#     )
# endforeach()
