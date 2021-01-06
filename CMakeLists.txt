cmake_minimum_required (VERSION 3.9)
project (PRIMULA)

# Locate headers
include_directories(src)

#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

#==============================================================================
# PRIMULA
#==============================================================================
include(CMake/KiLib.cmake.in)
include(CMake/CLI11.cmake.in)

add_executable(
    PRIMULA 
    src/main.cpp 
    src/primula++.cpp
    src/landslide.cpp
    src/Config.cpp
)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # using Clang
    message(STATUS "Clang detected, adding compile flags")
    target_compile_options(PRIMULA PRIVATE -Wall -Werror)
    # Enforce C++17 features
    target_compile_features(PRIMULA PRIVATE cxx_std_17)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # using GCC
    message(STATUS "GCC detected, adding compile flags")
    target_compile_options(PRIMULA PRIVATE -Wall -Werror)
    # Enforce C++17 features
    target_compile_features(PRIMULA PRIVATE cxx_std_17)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # using Visual Studio C++
    message(STATUS "msvc detected, adding compile flags")
    target_compile_options(PRIMULA PRIVATE  /std:c++latest /EHsc)
endif()

target_link_libraries(PRIMULA KiLib CLI11 tiff)
