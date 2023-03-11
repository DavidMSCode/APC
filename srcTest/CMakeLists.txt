include_directories(${APC_SOURCE_DIR}/include ${APC_SOURCE_DIR}/extern/cspice/include) #add header file locations for APC and CSPICE
find_package(OpenMP REQUIRED)
add_executable(main test.cpp) #build a an executable using the test src file


target_link_libraries(main APC) #link against APC library
target_link_libraries(main OpenMP::OpenMP_CXX)   #link openmp library