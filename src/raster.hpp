// ============================================================================
// Copyright (C) Denis Cohen-Corticchiato. All rights reserved.
//
// Authors: Denis Cohen-Corticchiato (DOHCC)
// Email:   denis.cohen@gmail.com
//
// This file is part of DOHCC C++ Library
// ============================================================================

#ifndef RASTER_HPP
#define RASTER_HPP

#include <cstdio>
#include <fstream>
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>

#include "point.hpp"

//=============================================================================
// CLASS Raster
//=============================================================================

class Raster { 

   public:

      // ... Constructors ...
      Raster() { };
      //Raster(const std::string & file);     // TO IMPLEMENT
      Raster(const Raster & raster)  // Copy constructor
      {
            this->ncols_ = raster.ncols_;
            this->nrows_ = raster.nrows_;
            this->ncells_ = raster.ncells_;
            this->cellsize_ = raster.cellsize_;
            this->cellsizex_ = raster.cellsizex_;
            this->cellsizey_ = raster.cellsizey_;
            this->xllcorner_ = raster.xllcorner_;
            this->yllcorner_ = raster.yllcorner_;
            this->nodata_value_ = raster.nodata_value_;
            this->points_ = raster.points_;
            this->attribute_.resize(ncells_);
      };
      //~Raster() { };                        // NEEDED?
	
      // ... Member functions ...
      bool Read(const std::string & file);
      bool Print(const std::string & file, const std::string & format, const int & precision);
      void Info();
      void AssignBoundaryId();

      void SetAttribute(const std::vector<double> & v) { attribute_ = v; }
      std::vector<double> GetAttribute() const { return attribute_; };
      double GetAttributeIndex(const int & i) const { return attribute_[i]; };

      // ... Data ...
      std::string file_;
      unsigned int ncols_{0};
      unsigned int nrows_{0};
      unsigned int ncells_{0};
      double xllcorner_{0.0};
      double yllcorner_{0.0};
      double xurcorner_{0.0};
      double yurcorner_{0.0};
      double cellsize_{0.0};
      double cellsizex_{0.0};
      double cellsizey_{0.0};
      double nodata_value_{9999};
      std::vector<double> attribute_{};    // Attribute value stored in raster
      std::string attribute_name_;
      std::vector<int> boundaryId_{};      // Boundary id
      std::vector<Point> points_{};        // Stores x and y coordinates (and z?)

//-----------------------------------------------------------------------------
   private:

};

//=============================================================================
// NON-MEMBER FUNCTIONS (NOT REALLY USED AT THE MOMENT)
//=============================================================================

bool Read2(Raster & raster, const std::string & file);
bool Print2(const Raster & raster, const std::string & file, const std::string & format, const int & precision);


#endif
