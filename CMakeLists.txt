cmake_minimum_required(VERSION 3.9)
project(BDHW3)

set(CMAKE_CXX_STANDARD 11)

find_package(OpenMP REQUIRED)

add_executable(BDHW3 main.cpp)
target_link_libraries(BDHW3 PUBLIC OpenMP::OpenMP_CXX)