// ============================================================================
// Copyright (C) PRIMULA++. All rights reserved.
//
// Authors: Denis Cohen-Corticchiato (DOHCC)
// Email:   denis.cohen@gmail.com
//
// This file is part of PRIMULA++
// ============================================================================

#include <csv.hpp>
#include <KiLib/Utils/Random.hpp>
#include <chrono>
#include <cmath>
#include <primula++.hpp>
#include <random>
#include <spdlog/spdlog.h>
#include <spdlog/stopwatch.h>
#include <stats.hpp>

KiLib::Raster Primula::CalcWetness(const KiLib::Raster &ks, const KiLib::Raster &z)
{
   KiLib::Raster W = KiLib::Raster::nodataLike(this->slope_);

   for (size_t i : this->validIndices)
   {
      W(i) = this->hydroModel.ComputeWetness(rainfall_, ks(i), z(i), slope_(i), twi_(i));
   }

   return W;
}

KiLib::Raster Primula::MDSTab_v2(
   const Landslide &slide, const KiLib::Raster &phi, const KiLib::Raster &m, const double &gamma_s,
   const KiLib::Raster &z, const KiLib::Raster &Crl, const KiLib::Raster &Crb)
{
   KiLib::Raster FS = KiLib::Raster::nodataLike(this->slope_);

   // calculate factor of safety for each raster cell
   for (size_t i : this->validIndices)
   {
      double SF = this->SFModel.ComputeSF(
         phi(i), m(i), z(i), Crl(i), Crb(i), this->slope_(i), this->slope_(i), gamma_s, slide.width_, slide.length_);

      FS(i) = SF;
   }

   return FS;
}

void Primula::ReadLandCover(const std::string &landCover)
{
   csv::CSVReader reader(landCover);

   for (csv::CSVRow &row : reader)
   {
      if (row[0].is_int())
         this->landcover[row[0].get<size_t>()] = {row[4].get<double>(), row[5].get<double>()};
   }
}

void Primula::ReadSoilDepth(const std::string &soilDepth)
{
   csv::CSVReader reader(soilDepth);

   for (csv::CSVRow &row : reader)
   {
      if (row[0].is_int())
         this->soilDepth[row[0].get<size_t>()] = {row[2].get<double>(), row[3].get<double>()};
   }
}

void Primula::ReadPhysProps(const std::string &physProps)
{
   csv::CSVReader reader(physProps);

   for (csv::CSVRow &row : reader)
   {
      if (!row[0].is_int())
         continue;

      this->physProps[row[0].get<size_t>()] = {
         row[2].get<double>(), row[3].get<double>(), // Min Gamma, Max Gamma
         row[4].get<double>(), row[5].get<double>(), // Min Phi, Max Phi
         row[6].get<double>(), row[7].get<double>(), // Min Cohesion, Max Cohesion
         row[8].get<double>(), row[9].get<double>(), // Min Ks, Max Ks
      };
   }
}


/**
 * @brief Uses the data which has been read and generates soil properties based on that
 *
 */
void Primula::GenerateSoilProperties()
{
   spdlog::stopwatch sw;
   
   for (auto const& [soilID, prop] : this->physProps)
   {
      std::vector<double> tPhi      = KiLib::Random::runif(this->num_landslides, 0, 1, this->engine);
      std::vector<double> tGamma    = KiLib::Random::runif(this->num_landslides, 0, 1, this->engine);
      std::vector<double> tKs       = KiLib::Random::runif(this->num_landslides, 0, 1, this->engine);
      std::vector<double> tCohesion = KiLib::Random::runif(this->num_landslides, 0, 1, this->engine);

      tPhi      = stats::qunif(tPhi, prop.minPhi, prop.maxPhi);
      tGamma    = stats::qunif(tGamma, prop.minGamma, prop.maxGamma);
      tKs       = stats::qunif(tKs, prop.minKs, prop.maxKs);
      tCohesion = stats::qunif(tCohesion, prop.minCohesion, prop.maxCohesion);

      std::transform(tGamma.begin(), tGamma.end(), tGamma.begin(), [&](double v) -> double { return v*1000.0; });

      this->phi[soilID]      = tPhi;
      this->gamma[soilID]    = tGamma;
      this->ks[soilID]       = tKs;
      this->cohesion[soilID] = tCohesion;
   }

   for (size_t i = 0; i < this->num_landslides; i++)
   {
      Landslide slide;

      // generate random landslide properties
      slide.area_   = pow(10, stats::qnorm(stats::runif(0, 1, this->engine), this->area_mu_, this->area_sigma_));
      auto l2w      = pow(10, stats::qnorm(stats::runif(0, 1, this->engine), this->l2w_mu_, this->l2w_sigma_));
      slide.width_  = sqrt((slide.area_ * 1.0) / (l2w * 1.0));
      slide.length_ = slide.width_ * l2w;
      this->landslide_.push_back(slide);
   }

   spdlog::info("Soil generation elapsed time: {}", sw);
}


/**
 * @brief Uses the generated soil properties to calculate a safety factor.
 *
 */
void Primula::CalculateSafetyFactor()
{
   spdlog::stopwatch sw;

   this->slope_.nodata_value = -9999;
   this->pr_failure_         = KiLib::Raster::nodataLike(this->slope_);

   for (size_t i : this->validIndices)
   {
      this->pr_failure_(i) = 0.0;
   }

   for (unsigned int i = 0; i < num_landslides; i++)
   {
      spdlog::info("Landslide {}", i);
      KiLib::Raster friction_angle = KiLib::Raster::zerosLike(this->slope_);
      KiLib::Raster permeability   = KiLib::Raster::zerosLike(this->slope_);
      KiLib::Raster depth          = KiLib::Raster::zerosLike(this->slope_);
      KiLib::Raster crl            = KiLib::Raster::zerosLike(this->slope_);
      KiLib::Raster crb            = KiLib::Raster::zerosLike(this->slope_);

      // go through each raster cell
      for (size_t j : this->validIndices)
      {
         // if soil 1 or 2, translate info to rasters
         double phiV      = this->phi.at(this->soil_type_(j))[i];
         double ksV       = this->ks.at(this->soil_type_(j))[i];
         // double cohesionV = this->Cohesion[this->soil_type_(j)][i]; Currently unused!!!!!!!!!!!!!!!
         // double gammaV    = this->gamma[this->soil_type_(j)][i]; Currently unused!!!!!!!!!!!!!!!

         friction_angle(j) = phiV * M_PI / 180.0;
         permeability(j)   = ksV * M_PI / 180.0;

         auto& [minLC, maxLC] = this->landcover.at(this->landuse(j));
         crl(j)               = stats::runif(minLC, maxLC, this->engine);

         auto &[minSD, maxSD] = this->soilDepth.at(this->soil_type_(j));
         depth(j)             = stats::runif(minSD, maxSD, this->engine);

         if (depth(j) >= 0.5)
            crb(j) = 0;
         else
            crb(j) = crl(j);
      }

      auto m = CalcWetness(permeability, depth);

      auto FS = MDSTab_v2(this->landslide_[i], friction_angle, m, this->gamma.at(1)[i], depth, crl, crb);
      for (size_t i : this->validIndices)
      {
         if (FS(i) < 1 && FS(i) > 0)
         {
            this->pr_failure_(i) += FS(i);
         }
      }
   }

   // get average of sum of failure probabilities
   for (size_t i : this->validIndices)
   {
      this->pr_failure_(i) /= num_landslides;
   }

   spdlog::info("Landslide generation elapsed time: {}", sw);
}

void Primula::syncValidIndices()
{
   spdlog::info("Validating raster data");

   for (size_t i = 0; i < this->slope_.nData; i++)
   {
      if (this->slope_(i) == this->slope_.nodata_value)
         continue;
      if (this->twi_(i) == this->twi_.nodata_value)
         continue;
      if (this->soil_type_(i) == this->soil_type_.nodata_value)
         continue;
      if (this->landuse(i) == this->landuse.nodata_value)
         continue;

      this->validIndices.push_back(i);
   }

   spdlog::info(
      "{} / {} Raster indices have valid data. Only computing on valid data.", this->validIndices.size(),
      this->slope_.nData);
}
