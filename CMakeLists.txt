cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD_REQUIRED 17)

project(AnasenSim)

set(ASIM_BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(ASIM_LIBRARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)

find_package(ROOT REQUIRED COMPONENTS GenVector)

add_subdirectory(vendor/catima)
add_subdirectory(src)