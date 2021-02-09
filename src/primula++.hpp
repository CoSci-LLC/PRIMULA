/**
 *  Copyright (c) 2020-2021 CoSci LLC, USA <software@cosci-llc.com>
 *  
 *  This file is part of PRIMULA.
 *  
 *  PRIMULA is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  PRIMULA is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with PRIMULA.  If not, see <https://www.gnu.org/licenses/>.
*/

#pragma once

#include <vector>

#include <KiLib/KiLib.hpp>
#include <landslide.hpp>
#include <random>

// Structure to conveniently store properties for each soil
struct physProp
{
   double minGamma;
   double maxGamma;

   double minPhi;
   double maxPhi;

   double minCohesion;
   double maxCohesion;

   double minKs;
   double maxKs;
};


class Primula
{
public:
   //==========================================================================
   // ... Constructors, Destructors ...
   //==========================================================================
   Primula(int n, size_t seed) : num_landslides(n)
   {
      this->engine = std::mt19937_64(seed); // this->Engine so our seed is consistent
   }

   //==========================================================================
   // ... Public Member Functions ...
   //==========================================================================
   void GenerateSoilProperties();
   void CalculateSafetyFactor();
   void ReadLandCover(const std::string &landCover);
   void ReadSoilDepth(const std::string &soilDepth);
   void ReadPhysProps(const std::string &physProps);

   // ... Member Data ...
   size_t num_landslides;

   // Model to calculate safety factor
   KiLib::Stability::SafetyFactor::MDSTAB SFModel;
   // Model to wetness
   KiLib::Hydrology::TopModel hydroModel;

   KiLib::Raster slope_;
   KiLib::Raster twi_;
   KiLib::Raster soil_type_;
   KiLib::Raster landuse;

   KiLib::Raster pr_failure_;

   std::vector<Landslide> landslide_;

   double rainfall_ = 0.160; // [m/day]

   // Normal distribution parameters for landslide area
   double area_mu_    = 2.017;          // [m^2]
   double area_sigma_ = sqrt(0.176466); // [m^2]

   // Normal distribution parameters for length-to-width ratio
   double l2w_mu_    = 0.1528;         // [dimensionless]
   double l2w_sigma_ = sqrt(0.037396); // [dimensionless]

   // Looks at all rasters and finds the indices where ALL rasters have valid data (i.e. not nodata_value), stores in a
   // vector called validIndices
   std::vector<size_t> validIndices;
   void                syncValidIndices();


   //==========================================================================
   // ... Private Member Data ...
   //==========================================================================

private:
   std::unordered_map<double, std::pair<double, double>> landcover;
   std::unordered_map<double, std::pair<double, double>> soilDepth;

   std::unordered_map<double, physProp> physProps;

   std::unordered_map<double, std::vector<double>> phi;
   std::unordered_map<double, std::vector<double>> gamma;
   std::unordered_map<double, std::vector<double>> ks;
   std::unordered_map<double, std::vector<double>> cohesion;

   std::mt19937_64 engine; // RNG engine so our results can be consistent

   KiLib::Raster CalcWetness(const KiLib::Raster &ks, const KiLib::Raster &z);
   KiLib::Raster MDSTAB(
      const Landslide &slide, const KiLib::Raster &phi, const KiLib::Raster &m, const KiLib::Raster &gamma_s,
      const KiLib::Raster &z, const KiLib::Raster &Crl, const KiLib::Raster &Crb);
};
