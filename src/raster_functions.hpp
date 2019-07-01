// ============================================================================
// Copyright (C) Denis Cohen-Corticchiato. All rights reserved.
//
// Authors: Denis Cohen-Corticchiato (DOHCC)
// Email:   denis.cohen@gmail.com
//
// This file is part of DOHCC C++ Library
// ============================================================================

#ifndef RASTER_FUNCTIONS_HPP
#define RASTER_FUNCTIONS_HPP

#include "raster.hpp"


Raster Slope_from_Raster(const Raster & rDEM);

bool FindRasterCell(const Raster & raster, const double & x, const double & y, std::vector<double> & results);

// New: to replace above
std::vector<double> FindRasterCell(const Raster & raster, const double & x, const double & y);
unsigned int FindRasterCellIndex(const Raster & raster, const double & x, const double & y);
std::vector<double> FindRasterCellCoord(const Raster & raster, const double & x, const double & y);

double BilinearInterpolation(const double & q11, 
                             const double & q12, 
                             const double & q21, 
                             const double & q22, 
                             const double & x1, 
                             const double & x2, 
                             const double & y1, 
                             const double & y2, 
                             const double & x, 
                             const double & y);

#endif
