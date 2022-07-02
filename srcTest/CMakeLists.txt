include_directories(${APC_SOURCE_DIR}/include ${APC_SOURCE_DIR}/extern/cspice/include) #add header file locations for APC and CSPICE
SET(CMAKE_OSX_ARCHITECTURES x86_64)
add_executable(test test.cpp) #build a an executable using the test src file

target_link_libraries(test APC) #link against APC library