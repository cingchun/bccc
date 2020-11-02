cmake_minimum_required (VERSION 3.8)

set(CMAKE_C_COMPILER icc)
set(CMAKE_CXX_COMPILER icpc)


set(proj ta)
project (${proj} CXX)

set (${proj}_VERSION_MAJOR 1)
set (${proj}_VERSION_MINOR 0)
set (${proj}_VERSION_PATCH 0)


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

aux_source_directory(${proj} SRC)

# OpenMP
find_package(OpenMP)
if (OPENMP_FOUND OR OpenMP_CXX_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

SET(BUILD_SHARED_LIBS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_BUILD_TYPE "Release")
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O2 -Wall")


add_library(${proj} pub.cpp ${SRC})
set_target_properties(${proj} PROPERTIES VERSION 1.0 SOVERSION 1)

