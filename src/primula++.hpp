#pragma once
// ============================================================================
// Copyright (C) PRIMULA. All rights reserved.
//
// Authors: Denis Cohen-Corticchiato (DOHCC)
// Email:   denis.cohen@gmail.com
//
// This file is part of PRIMULA
// ============================================================================

#include <vector>

#include <KiLib/KiLib.hpp>
#include <landslide.hpp>
#include <random>

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
   // ... Public Accessor Functions ...
   //==========================================================================


   //==========================================================================
   // ... Public Member Functions ...
   //==========================================================================
   void GenerateSoilProperties();
   void CalculateSafetyFactor();
   void ReadSoilDataset(const std::string &soil_data, const std::string &root_data);

   // ... Static Member Data ...
   static constexpr double gravity_        = 9.81;         // [m/s^2]
   static constexpr double soil_density_   = 1834.8624;    // [kg/m^3] Schwarz M. R Code
   static constexpr double transmissivity_ = 0.0002644655; // [m/s] Fixed for the moment. From R code
   static constexpr double rain_intensity_ =
      0.001 * 70 / 3600; // [m/s] >> Eventually should live outside of primula class <<


   // ... Member Data ...
   size_t num_landslides;

   // Model to calculate safety factor
   KiLib::Stability::SafetyFactor::MDSTab SFModel;
   // Model to wetness
   KiLib::Hydrology::TopModel hydroModel;

   KiLib::Raster slope_;
   KiLib::Raster probslope_;
   KiLib::Raster twi_;
   KiLib::Raster soil_type_;
   KiLib::Raster soil_depth_;
   KiLib::Raster dusaf_;
   KiLib::Raster pr_failure_;

   std::vector<Landslide>           landslide_;
   std::vector<int>                 soil_id_;
   std::vector<std::vector<double>> z_;

   std::vector<double> Cr_grassland_;
   std::vector<double> Cr_shrubland_;

   double veg_weight_ = 70.0;  // [kg/m2]
   double rainfall_   = 0.100; // [m/day]

   // Normal distribution parameters for landslide area
   double area_mu_    = 2.017;          // [m^2]
   double area_sigma_ = sqrt(0.176466); // [m^2]

   // Normal distribution parameters for length-to-width ratio
   double l2w_mu_    = 0.1528;         // [dimensionless]
   double l2w_sigma_ = sqrt(0.037396); // [dimensionless]


   //==========================================================================
   // ... Private Member Data ...
   //==========================================================================

private:
   std::vector<size_t> iteration_index;

   std::vector<double> Pa400;
   std::vector<double> Pa200;
   std::vector<double> Fs800;
   std::vector<double> Fs200;
   std::vector<double> Cs150;
   std::vector<double> Mf600;
   std::vector<double> Mf300;

   std::vector<double> phi1;
   std::vector<double> phi2;
   std::vector<double> gamma1;
   std::vector<double> ks1;
   std::vector<double> ks2;

   std::mt19937_64 engine; // Engine so our can be consistent

   KiLib::Raster CalcWetness(const KiLib::Raster &ks, const KiLib::Raster &z);
   KiLib::Raster MDSTab_v2(
      const Landslide &slide, const KiLib::Raster &phi, const KiLib::Raster &m, const double &gamma_s,
      const KiLib::Raster &z, const KiLib::Raster &Crl, const KiLib::Raster &Crb);
};
