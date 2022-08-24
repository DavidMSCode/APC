include_directories(${APC_SOURCE_DIR}/include)
add_executable(MatrixBuilder matrix_builder.cpp) #build an executable using the matrix builder src
target_link_libraries(MatrixBuilder APC) #link against APC library