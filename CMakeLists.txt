cmake_minimum_required(VERSION 3.0)
project(BDHW3)

set(CMAKE_CXX_STANDARD 11)

#find_package(OpenMP REQUIRED)

add_executable(BDHW3 main.cpp)
#target_compile_options(BDHW3 PUBLIC ${OpenMP_CXX_FLAGS})