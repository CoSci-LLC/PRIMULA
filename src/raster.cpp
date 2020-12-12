// ============================================================================
// Copyright (C) Denis Cohen-Corticchiato. All rights reserved.
//
// Authors: Denis Cohen-Corticchiato (DOHCC)
// Email:   denis.cohen@gmail.com
//
// This file is part of DOHCC C++ Library
// ============================================================================

#include "raster.hpp"

#include <algorithm>
#include <iomanip>


//=============================================================================
bool Raster::Read(const std::string &file)
{
   file_ = file;
   std::cout << "Raster:Read \"" << file_ << "\"" << std::endl;

   // Open data file
   std::ifstream fin;
   fin.open(file_);
   if (!fin.is_open()) {
      std::cerr << "  File \"" + file_ + "\" open failed" << std::endl;
      return false;
   }

   // Read header
   std::string str;
   fin >> str >> ncols_;
   fin >> str >> nrows_;
   fin >> str >> xllcorner_;
   fin >> str >> yllcorner_;
   fin >> str >> cellsize_;
   fin >> str >> nodata_value_;

   ncells_    = ncols_ * nrows_;
   cellsizex_ = cellsize_;
   cellsizey_ = cellsize_;
   xurcorner_ = xllcorner_ + cellsize_ * static_cast<double>(ncols_);
   yurcorner_ = yllcorner_ + cellsize_ * static_cast<double>(nrows_);

   // Store raster attribute values. Values are stored by rows from bottom so
   // 1st point is lower left, last point is upper right
   attribute_.resize(ncells_);
   points_.resize(ncells_);
   for (unsigned int i = nrows_; i >= 1; i--) {
      for (unsigned int j = 0; j < ncols_; j++) {
         const unsigned int count = (i - 1) * ncols_ + j;
         double val;
         fin >> val;
         attribute_.at(count) = val;
         points_.at(count).x_ = xllcorner_ + (cellsizex_ / 2) + static_cast<double>(j) * cellsizex_;
         points_.at(count).y_ = yllcorner_ + (cellsizey_ / 2) + static_cast<double>(i - 1) * cellsizey_;
      }
   }

   AssignBoundaryId();

   std::cout << "Raster:Read ... Done" << std::endl;
   fin.close();
   return true;
}

//==============================================================================
bool Raster::Print(const std::string &file, const std::string &format, const int &precision)
{
   std::cout << "Raster:Print \"" << file << "\"" << std::endl;

   std::string str = format;
   std::transform(str.begin(), str.end(), str.begin(), ::tolower);

   if (str == "asc") {
      std::ofstream fout;
      fout.open(file);
      if (!fout.is_open()) {
         std::cout << "  File \"" + file + "\" open failed" << std::endl;
         return false;
      }
      fout << std::fixed << std::setprecision(4);
      fout << "ncols " << ncols_ << std::endl;
      fout << "nrows " << nrows_ << std::endl;
      fout << "xllcorner " << xllcorner_ << std::endl;
      fout << "yllcorner " << yllcorner_ << std::endl;
      fout << "cellsize " << cellsize_ << std::endl;
      fout << "NODATA_value -9999" << std::endl;
      fout << std::fixed << std::setprecision(precision);
      for (auto j = nrows_; j > 0; --j) {
         const auto start = (j - 1) * ncols_ + 1;
         const auto end   = j * ncols_;
         for (auto i = start; i <= end; ++i) {
            //          fout << values_.at(i-1) << " ";
            fout << attribute_.at(i - 1) << " ";
         }
         fout << std::endl;
      }
      fout.close();
   } else if (str == "xyz") {
      std::ofstream fout;
      fout.open(file);
      if (!fout.is_open()) {
         std::cout << "  File \"" + file + "\" open failed" << std::endl;
         return false;
      }

      fout << "VARIABLES = \"X\" \"Y\" \"" << attribute_name_ << "\"" << std::endl;
      fout << "ZONE T = \"" << attribute_name_ << "\"" << std::endl;
      fout << std::fixed << std::setprecision(precision);

      for (size_t i = 0; i < points_.size(); ++i) {
         fout << points_.at(i).x_ << " " << points_.at(i).y_ << " " << attribute_.at(i) << std::endl;
      }
   } else {
      std::cout << "Raster write file format not available. Abort." << std::endl;
      return false;
   }

   std::cout << "Raster:Print ... Done" << std::endl;
   return true;
}

//=============================================================================
void Raster::Info()
{
   std::cout << "Raster:Info ... Start" << std::endl;
   std::cout << "Raster:Info File: " << file_ << std::endl;
   std::cout << "File: " << file_ << std::endl;
   std::cout << "Size is " << ncols_ << " x " << nrows_ << std::endl;
   std::cout << "Cell size is " << cellsize_ << std::endl;
   std::cout << std::fixed << std::setprecision(4);
   std::cout << "Upper Left  ( " << xllcorner_ << ", " << yurcorner_ << ")" << std::endl;
   std::cout << "Lower Left  ( " << xllcorner_ << ", " << yllcorner_ << ")" << std::endl;
   std::cout << "Upper Right ( " << xurcorner_ << ", " << yurcorner_ << ")" << std::endl;
   std::cout << "Lower Right ( " << xurcorner_ << ", " << yllcorner_ << ")" << std::endl;
   std::cout << "Raster:Info ... Done" << std::endl;
}

//=============================================================================
void Raster::AssignBoundaryId()
{
   std::cout << "Raster:AssignBoundaryId \"" << file_ << "\"" << std::endl;
   //---------------------------
   // Assignment of boundary id
   //---------------------------
   // Interior node   = 0
   // Bottom boundary = 1
   // Right boundary  = 2
   // Top boundary    = 3
   // Left boundary   = 4
   // Corner cells    = 5 to 8 (counterclockwise from bottom-left)
   //
   //   8--3--3--3--3--3--7
   //   |  |  |  |  |  |  |
   //   4--0--0--0--0--0--2
   //   |  |  |  |  |  |  |
   //   4--0--0--0--0--0--2
   //   |  |  |  |  |  |  |
   //   4--0--0--0--0--0--2
   //   |  |  |  |  |  |  |
   //   5--1--1--1--1--1--6
   //

   boundaryId_.resize(ncells_);

   // Boundary 1: bottom
   for (unsigned int i = 0; i < ncols_; i++) {
      boundaryId_.at(i) = 1;
   }

   // Boundary 2: right
   for (unsigned int i = ncols_ - 1; i <= ncells_; i = i + ncols_) {
      if (boundaryId_.at(i) == 1) {
         boundaryId_.at(i) = 6; // Bottom-right corner cell
      } else {
         boundaryId_.at(i) = 2;
      }
   }

   // Boundary 3: top
   for (unsigned int i = ncells_ - ncols_; i < ncells_; i++) {
      if (boundaryId_.at(i) == 2) {
         boundaryId_.at(i) = 7; // Upper-right corner cell
      } else {
         boundaryId_.at(i) = 3;
      }
   }

   // Boundary 4: left
   for (unsigned int i = 0; i <= ncells_ - ncols_; i = i + ncols_) {
      if (boundaryId_.at(i) == 3) {
         boundaryId_.at(i) = 8; // Upper-left corner cell
      } else if (boundaryId_.at(i) == 1) {
         boundaryId_.at(i) = 5; // Bottom-left corner cell
      } else {
         boundaryId_.at(i) = 4;
      }
   }
   std::cout << "Raster:AssignBoundaryId ... Done" << std::endl;
}
