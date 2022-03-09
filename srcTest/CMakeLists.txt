include_directories(${APC_SOURCE_DIR}/include ${APC_SOURCE_DIR}/extern/cspice/include)
SET(CMAKE_OSX_ARCHITECTURES x86_64)
add_executable(test test.cpp)

target_link_libraries(test APC)