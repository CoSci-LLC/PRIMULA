cmake_minimum_required (VERSION 3.9)
project (PRIMULA VERSION 1.0.0 LANGUAGES CXX)

# Locate headers
include_directories(src)

#==============================================================================
# PRIMULA
#==============================================================================
include(CMake/KiLib.cmake.in)
include(CMake/CLI11.cmake.in)
include(CMake/CSV.cmake.in)

find_package(OpenMP)

add_executable(
    PRIMULA 
    src/main.cpp 
    src/primula++.cpp
    src/Config.cpp
)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # using Clang
    message(STATUS "Clang detected, adding compile flags")
    target_compile_options(PRIMULA PRIVATE -Wall)
    # Enforce C++17 features
    target_compile_features(PRIMULA PRIVATE cxx_std_17)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # using GCC
    message(STATUS "GCC detected, adding compile flags")
    target_compile_options(PRIMULA PRIVATE -Wall)
    # Enforce C++17 features
    target_compile_features(PRIMULA PRIVATE cxx_std_17)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # using Visual Studio C++
    message(STATUS "msvc detected, adding compile flags")
    target_compile_options(PRIMULA PRIVATE  /std:c++latest /EHsc)
endif()

target_link_libraries(PRIMULA KiLib CLI11 tiff csv OpenMP::OpenMP_CXX --coverage)
