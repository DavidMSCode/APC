#documentation tutorial from https://github.com/smrfeld/cpp_doxygen_sphinx/blob/master/README_doxygen.md  

cmake_minimum_required(VERSION 3.13.0)
# set(CMAKE_OSX_ARCHITECTURES "x86_64" CACHE INTERNAL "" FORCE)
cmake_policy(SET CMP0091 NEW) #allow setting MSVC runtime library instead of using cmake cache
project(APC VERSION 0.1.0)

add_subdirectory(src)     #SOurce code for library
add_subdirectory(srcMB)   #SOurce code for quadrature matrix builder
#Test code
include(CTest)
add_subdirectory(APCTests/src)  #Tests for library
#Compile Flags
set(CMAKE_CXX_STANDARD 14)
set_property(TARGET APC CombinedTests main_bootstrap PROPERTY CXX_STANDARD 14) #sets C++ standard to c++14 (ignored on Windows)

#[[choose runtime library (only Windows). CSPICE is statically linked with the release runtime flags /MT (libcmt.lib) (no debug)
When compiling for debug on Windows LNK4098 warning will occur but it will still be possible to run the resulting binary.]]
set_property(TARGET APC CombinedTests MatrixBuilder PROPERTY
  MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
enable_testing()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -fsanitize=address -fno-omit-frame-pointer")
set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fsanitize=address")

set(ENV{OMP_NUM_THREADS} 4) #set number of threads for OpenMP


set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

# Doxygen (Ignored if doxygen is not installed)
# look for Doxygen package
find_package(Doxygen)

if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs_doxygen/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.out)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")

    # Note: do not put "ALL" - this builds docs together with application EVERY TIME!
    add_custom_target(docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM )
else (DOXYGEN_FOUND)
  message("Doxygen needs to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)

