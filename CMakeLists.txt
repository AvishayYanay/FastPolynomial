cmake_minimum_required(VERSION 3.9)
project(fastpoly)

set(CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall")

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

include_directories(/path-to/NTL/ntl-10.5.0/include)

add_executable(fastpoly main.cpp recursive.cpp recursive.h utils.cpp utils.h iterative.cpp iterative.h)


target_link_libraries(fastpoly /path-to/NTL/ntl-10.5.0/src/ntl.a)
target_link_libraries(fastpoly gmp gmpxx)


