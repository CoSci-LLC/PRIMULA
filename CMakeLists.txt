cmake_minimum_required (VERSION 3.9)
project (PRIMULA)

#==============================================================================
# Packages / Libraries
#==============================================================================
find_package(OpenMP)

# Locate headers
include_directories(src)

#==============================================================================
# PRIMULA
#==============================================================================
include(CMake/KiLib.cmake.in)
add_executable(
    PRIMULA 
    src/main.cpp 
    src/primula++.cpp
    src/raster.cpp
    src/landslide.cpp
)

target_link_libraries(PRIMULA OpenMP::OpenMP_CXX KiLib)
