// ============================================================================
// Copyright (C) PRIMULA. All rights reserved.
//
// Authors: Denis Cohen-Corticchiato (DOHCC)
// Email:   denis.cohen@gmail.com
//
// This file is part of PRIMULA
// ============================================================================

#ifndef PRIMULA_HPP
#define PRIMULA_HPP

#include <vector>

#include "landslide.hpp"
#include <KiLib/KiLib.hpp>

class Primula
{
public:
   //==========================================================================
   // ... Constructors, Destructors ...
   //==========================================================================
   Primula();

   //==========================================================================
   // ... Public Member Functions ...
   //==========================================================================
   void GenerateLandslides(const std::string &file, const unsigned int &num_landslides);
   void ReadCSV(const std::string &file, const unsigned int &num_landslides);


   // ... Member Data ...
   KiLib::Raster slope_;
   KiLib::Raster probslope_;
   KiLib::Raster twi_;
   KiLib::Raster soil_type_;
   KiLib::Raster soil_depth_;
   KiLib::Raster dusaf_;
   std::vector<Landslide> landslide_;
   std::vector<int> soil_id_;
   std::vector<std::vector<double>> z_;
   KiLib::Raster pr_failure_;

   std::vector<double> Crl_Fs200_;
   std::vector<double> Crl_Fs800_;
   std::vector<double> Crl_Pa200_;
   std::vector<double> Crl_Pa400_;
   std::vector<double> Crl_Mf300_;
   std::vector<double> Crl_Mf600_;
   std::vector<double> Crl_Cs150_;
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
   KiLib::Raster TopModel_v3(const KiLib::Raster &ks, const KiLib::Raster &z);
   KiLib::Raster MDSTab_v2(
      const Landslide &slide, const KiLib::Raster &phi, const KiLib::Raster &m, const double &gamma_s,
      const KiLib::Raster &z, const KiLib::Raster &Crl, const KiLib::Raster &Crb);
};

#endif
