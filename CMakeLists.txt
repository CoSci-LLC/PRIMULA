cmake_minimum_required (VERSION 3.9)
project (PRIMULA)

#==============================================================================
# Packages / Libraries
#==============================================================================
find_package(OpenMP)

# Locate headers
include_directories(src)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#==============================================================================
# PRIMULA
#==============================================================================
include(CMake/KiLib.cmake.in)
include(CMake/libtiff.cmake.in)
include(CMake/CLI11.cmake.in)


add_executable(
    PRIMULA 
    src/main.cpp 
    src/primula++.cpp
    src/landslide.cpp
    src/Config.cpp
    src/tiff.cpp
)

target_link_libraries(PRIMULA OpenMP::OpenMP_CXX KiLib CLI11 tiff)
