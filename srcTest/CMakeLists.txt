include_directories(${APC_SOURCE_DIR}/include ${APC_SOURCE_DIR}/extern/cspice/include) #add header file locations for APC and CSPICE
find_package(OpenMP REQUIRED)
add_executable(test test.cpp) #build a an executable using the test src file


target_link_libraries(test APC) #link against APC library
target_link_libraries(test OpenMP::OpenMP_CXX)   #link openmp library