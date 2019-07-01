// =============================================================================
// Copyright (C) Denis Cohen-Corticchiato. All rights reserved.
//
// Authors: Denis Cohen-Corticchiato (DOHCC)
// Email:   denis.cohen@gmail.com
//
// This file is part of DOHCC C++ Library
// =============================================================================

#include "raster_functions.hpp"

//==============================================================================

bool Read2(Raster & raster, const std::string & file)
{
   return false;
}

//==============================================================================

bool Print2(const Raster & raster, const std::string & file, const std::string & format, const int & precision)
{
   return false;
}

//==============================================================================

Raster Slope_from_Raster(const Raster & rdem)
{
   std::cout << "Raster:Slope_from_Raster \"" << rdem.file_ << "\"" << std::endl;

   // Copy values into slope raster
   Raster rslope(rdem);

   auto boundaryId = rdem.boundaryId_;
   auto itb = boundaryId.begin();

   std::vector<double> slope(rslope.ncells_);
   auto its = slope.begin();

   // Loop over elevation using iterators
   auto elevation = rdem.GetAttribute();
   for (auto ite = elevation.begin(); ite < elevation.end(); ++ite)
   {
      auto fx = 0.0;
      auto fy = 0.0;

      auto ite_2 = ite;                                      // Cell on bottom 
      auto ite_1 = ite;                                      // Cell on bottom-left
      auto ite_3 = ite;                                      // Cell on bottom-right
      auto ite_4 = ite;                                      // Cell on left
      auto ite_6 = ite;                                      // Cell on right
      auto ite_8 = ite;                                      // Cell on top
      auto ite_7 = ite;                                      // Cell on top-left
      auto ite_9 = ite;                                      // Cell on top-right
      if (*itb == 0)
      {
         ite_2 -= rslope.ncols_;
         ite_1 = ite_2; --ite_1;
         ite_3 = ite_2; ++ite_3;
         --ite_4;
         ++ite_6;
         std::advance(ite_8,rslope.ncols_);
         ite_7 = ite_8; --ite_7;
         ite_9 = ite_8; ++ite_9;
      }
      else if (*itb == 1)
      {
         --ite_4;
         ++ite_6;
         std::advance(ite_8,rslope.ncols_);
         ite_7 = ite_8; --ite_7;
         ite_9 = ite_8; ++ite_9;
      }
      else if (*itb == 2)
      {
         ite_2 -= rslope.ncols_;
         ite_1 = ite_2; --ite_1;
         --ite_4;
         std::advance(ite_8,rslope.ncols_);
         ite_7 = ite_8; --ite_7;
      }
      else if (*itb == 3)
      {
         ite_2 -= rslope.ncols_;
         ite_1 = ite_2; --ite_1;
         ite_3 = ite_2; ++ite_3;
         --ite_4;
         ++ite_6;
      }
      else if (*itb == 4)
      {
         ite_2 -= rslope.ncols_;
         ite_3 = ite_2; ++ite_3;
         ++ite_6;
         std::advance(ite_8,rslope.ncols_);
         ite_9 = ite_8; ++ite_9;
      }
      else if (*itb == 5)
      {
         ++ite_6;
         std::advance(ite_8,rslope.ncols_);
         ite_9 = ite_8; ++ite_9;
      }
      else if (*itb == 6)
      {
         --ite_4;
         std::advance(ite_8,rslope.ncols_);
         ite_7 = ite_8; --ite_7;
      }
      else if (*itb == 7)
      {
         ite_2 -= rslope.ncols_;
         ite_1 = ite_2; --ite_1;
         --ite_4;
      }
      else if (*itb == 8)
      {
         ite_2 -= rslope.ncols_;
         ite_3 = ite_2; ++ite_3;
         ++ite_6;
      }
      auto right = *ite_9 + (2 * *ite_6) + *ite_3;
      auto left = *ite_7 + (2 * *ite_4) + *ite_1;
      auto upper = *ite_7 + (2 * *ite_8) + *ite_9;
      auto lower = *ite_1 + (2 * *ite_2) + *ite_3;
      // Compute gradf: fx = df/dx and fy = df/dy using second order finite different
      fx = (right - left) / (8.0 * rslope.cellsizex_);
      fy = (upper - lower) / (8.0 * rslope.cellsizey_);
      // Compute slope (rise/run)
      //const auto slope = sqrt(fx*fx + fy*fy);
      //(*its) = slope;

//    const auto slope = 
      (*its) = sqrt(fx*fx + fy*fy);
      ++itb;
      ++its;
   }

   rslope.SetAttribute(slope);
   std::cout << "Raster:Slope_from_Raster ... Done" << std::endl;
   return rslope;
}

//==============================================================================

bool FindRasterCell(const Raster & raster, const double & x, const double & y, std::vector<double> & results)
{
   auto x0 = raster.xllcorner_ + (raster.cellsizex_/2);
   auto y0 = raster.yllcorner_ + (raster.cellsizey_/2);
   auto xn = raster.xurcorner_ - (raster.cellsizex_/2);
   auto yn = raster.yurcorner_ - (raster.cellsizey_/2);

   if ( x < raster.xllcorner_ || x > raster.xurcorner_ || y < raster.yllcorner_ || y > raster.yurcorner_)
   {
      std::cout << "Point ( " << x << ", " << y << ") outside of raster. Abort" << std::endl;
      return false;
   }
   auto dx = raster.cellsize_;
   auto dy = raster.cellsize_;
   auto ncols = raster.ncols_;
   auto nrows = raster.nrows_;

   auto i = static_cast<unsigned int>(floor((x - x0) / dx));
   auto j = static_cast<unsigned int>(floor((y - y0) / dy));
   if (i < 0) i = 0;
   //if (i >= ncols) i = ncols-1;
   if (j < 0) j = 0;
   //if (j >= nrows) j = nrows-1;

   auto x1 = x0 + static_cast<double>(i) * dx;
   auto x2 = x1 + dx;
   auto y1 = y0 + static_cast<double>(j) * dy;
   auto y2 = y1 + dy;
   auto z11 = raster.GetAttributeIndex(ncols * j + i);
   auto z12 = raster.GetAttributeIndex(ncols * j + i + 1);
   auto z21 = raster.GetAttributeIndex(ncols * (j+1) + i);  
   auto z22 = raster.GetAttributeIndex(ncols * (j+1) + i + 1);

   std::vector<double> tmp{ x1, x2, y1, y2, z11, z12, z21, z22 };
   results = std::move(tmp);

   return true;
}

//==============================================================================

std::vector<double> FindRasterCell(const Raster & raster, const double & x, const double & y)
{
   auto x0 = raster.xllcorner_ + (raster.cellsizex_/2);
   auto y0 = raster.yllcorner_ + (raster.cellsizey_/2);
   auto xn = raster.xurcorner_ - (raster.cellsizex_/2);
   auto yn = raster.yurcorner_ - (raster.cellsizey_/2);

   if ( x < raster.xllcorner_ || x > raster.xurcorner_ || y < raster.yllcorner_ || y > raster.yurcorner_)
   {
      std::cout << "Point ( " << x << ", " << y << ") outside of raster. Abort" << std::endl;
      exit (EXIT_FAILURE);
   }
   auto dx = raster.cellsize_;
   auto dy = raster.cellsize_;
   auto ncols = raster.ncols_;
   auto nrows = raster.nrows_;

   unsigned int i, j;
   if (i < 0) i = 0;
   else i = static_cast<unsigned int>(floor((x - x0) / dx));
   if (j < 0) j = 0;
   else j = static_cast<unsigned int>(floor((y - y0) / dy));

   auto x1 = x0 + static_cast<double>(i) * dx;
   auto x2 = x1 + dx;
   auto y1 = y0 + static_cast<double>(j) * dy;
   auto y2 = y1 + dy;
   auto z11 = raster.GetAttributeIndex(ncols * j + i);
   auto z12 = raster.GetAttributeIndex(ncols * j + i + 1);
   auto z21 = raster.GetAttributeIndex(ncols * (j+1) + i);  
   auto z22 = raster.GetAttributeIndex(ncols * (j+1) + i + 1);

   std::vector<double> result{ x1, x2, y1, y2, z11, z12, z21, z22 };
   return result;
}

//==============================================================================

unsigned int FindRasterCellIndex(const Raster & raster, const double & x, const double & y)
{
   auto x0 = raster.xllcorner_ + (raster.cellsizex_/2);
   auto y0 = raster.yllcorner_ + (raster.cellsizey_/2);
   auto xn = raster.xurcorner_ - (raster.cellsizex_/2);
   auto yn = raster.yurcorner_ - (raster.cellsizey_/2);

   if ( x < raster.xllcorner_ || x > raster.xurcorner_ || y < raster.yllcorner_ || y > raster.yurcorner_)
   {
      std::cout << "Point ( " << x << ", " << y << ") outside of raster. Abort" << std::endl;
      exit (EXIT_FAILURE);
   }
   auto dx = raster.cellsize_;
   auto dy = raster.cellsize_;
   auto ncols = raster.ncols_;
   auto nrows = raster.nrows_;

   unsigned int i, j;
   if (x < x0) i = 0;
   else if (x > xn) i = ncols - 1;
   else i = static_cast<unsigned int>(round((x - x0) / dx));
   if (y < y0) j = 0;
   else if (y > yn) j = nrows - 1;
   else j = static_cast<unsigned int>(round((y - y0) / dy));

   auto i1 = ncols * j + i;
   /*auto i2 = ncols * j + i + 1;
   auto i3 = ncols * (j+1) + i;  
   auto i4 = ncols * (j+1) + i + 1;

   std::vector<unsigned int> result{ i1, i2, i3, i4 };*/
   return i1;
}

//==============================================================================

std::vector<double> FindRasterCellCoord(const Raster & raster, const double & x, const double & y)
{
   auto x0 = raster.xllcorner_ + (raster.cellsizex_/2);
   auto y0 = raster.yllcorner_ + (raster.cellsizey_/2);
   auto xn = raster.xurcorner_ - (raster.cellsizex_/2);
   auto yn = raster.yurcorner_ - (raster.cellsizey_/2);

   if ( x < raster.xllcorner_ || x > raster.xurcorner_ || y < raster.yllcorner_ || y > raster.yurcorner_)
   {
      std::cout << "Point ( " << x << ", " << y << ") outside of raster. Abort" << std::endl;
      exit (EXIT_FAILURE);
   }
   auto dx = raster.cellsizex_;
   auto dy = raster.cellsizey_;
   auto ncols = raster.ncols_;
   auto nrows = raster.nrows_;

   unsigned int i, j;
   if (x < x0) i = 0;
   else if (x > xn) i = ncols - 1;
   else i = static_cast<unsigned int>(round((x - x0) / dx));
   if (y < y0) j = 0;
   else if (y > yn) j = nrows - 1;
   else j = static_cast<unsigned int>(round((y - y0) / dy));

   auto x1 = x0 + static_cast<double>(i) * dx;
   auto x2 = x1 + dx;
   auto y1 = y0 + static_cast<double>(j) * dy;
   auto y2 = y1 + dy;

   std::vector<double> result{ x1, x2, y1, y2};
   return result;
}

//==============================================================================

double BilinearInterpolation(const double & q11, 
                             const double & q12, 
                             const double & q21, 
                             const double & q22, 
                             const double & x1, 
                             const double & x2, 
                             const double & y1, 
                             const double & y2, 
                             const double & x, 
                             const double & y)
{
    double x2x1, y2y1, x2x, y2y, yy1, xx1;
    x2x1 = x2 - x1;
    y2y1 = y2 - y1;
    x2x = x2 - x;
    y2y = y2 - y;
    yy1 = y - y1;
    xx1 = x - x1;
    return 1.0 / (x2x1 * y2y1) * (
        q11 * x2x * y2y +
        q21 * xx1 * y2y +
        q12 * x2x * yy1 +
        q22 * xx1 * yy1
    );
}

