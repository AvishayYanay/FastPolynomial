cmake_minimum_required(VERSION 3.9)
project(fastpoly)

set(CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall")

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_library(fastpoly recursive.cpp recursive.h utils.cpp utils.h iterative.cpp iterative.h)

add_executable(test_fastpoly main.cpp)

target_link_libraries(test_fastpoly fastpoly)
target_link_libraries(test_fastpoly ntl)
target_link_libraries(test_fastpoly gmp gmpxx)
