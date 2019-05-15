cmake_minimum_required(VERSION 3.13)
project(parallel_lib)

set(CMAKE_CXX_STANDARD 14)

find_package(Threads REQUIRED)

add_executable(main main.cpp matrix.cpp matrix.hpp)
target_link_libraries(main Threads::Threads)