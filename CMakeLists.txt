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
add_executable(
    PRIMULA++ 
    src/main.cpp 
    src/primula++.cpp
    src/raster.cpp
    src/raster_functions.cpp
    src/landslide.cpp
)

target_link_libraries(PRIMULA++ OpenMP::OpenMP_CXX)
