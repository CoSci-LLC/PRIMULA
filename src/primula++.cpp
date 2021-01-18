// ============================================================================
// Copyright (C) PRIMULA++. All rights reserved.
//
// Authors: Denis Cohen-Corticchiato (DOHCC)
// Email:   denis.cohen@gmail.com
//
// This file is part of PRIMULA++
// ============================================================================

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
   KiLib::Raster W = KiLib::Raster::nodataLike(this->probslope_);

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
   KiLib::Raster FS = KiLib::Raster::nodataLike(this->probslope_);

   // calculate factor of safety for each raster cell
   for (size_t i : this->validIndices)
   {
      double SF = this->SFModel.ComputeSF(
         phi(i), m(i), z(i), Crl(i), Crb(i), this->probslope_(i), this->slope_(i), gamma_s, slide.width_, slide.length_);

      FS(i) = SF;
   }

   return FS;
}

/**
 * @brief Reads the data from the appropriate CSV files and populates hidden datastructures in the class
 *
 * @param soil_data The relative path to the soil data
 * @param root_data The relative path to the root data
 */
void Primula::ReadSoilDataset(const std::string &soil_data, const std::string &root_data)
{

   spdlog::stopwatch sw;

   size_t last  = 0;
   size_t count = 0;
   int    cod   = -1;
   int    prof  = -1;
   char   state = '\n';

   std::unordered_map<std::string, size_t> col_pos;
   std::stringstream                       ss;
   std::string                             line;
   std::string                             word;

   std::ifstream fin;

   // TODO: Consider replacing this with a CSV parsing library

   // ------------------------------
   // Reading Soil Data
   // ------------------------------
   fin.open(soil_data, std::ios::in);
   if (!fin.is_open())
   {
      spdlog::error("File '{}' failed to open", soil_data);
      exit(EXIT_FAILURE);
   }

   count = 0;
   getline(fin, line);
   ss.str(line);

   while (ss.good() && getline(ss, word, ','))
   {
      col_pos[word] = count++;
   }

   while (getline(fin, line))
   {
      // getline strips the last '\n' which is necessary to check for data in the last column
      line += '\n';
      state = '\n';
      last  = 0;
      count = 0;
      cod   = -1;
      prof  = -1;

      for (size_t it = 0; it < line.size(); it++)
      {

         switch (line[it])
         {
         // Check if it is necessary to consider quataion
         case '"':
         case '\'':
            if (state == line[it])
               state = '\n';
            else if (state == '\n')
               state = line[it];
            break;
         case ',':
            if (state != '\n')
               break;
         case '\n':
            if (count == col_pos["COD_UTS1"])
               cod = std::stoi(line.substr(last, it - last));
            else if (count == col_pos["PROF_UTILE"])
               prof = std::stoi(line.substr(last, it - last));

            state = '\n';
            last  = it + 1;
            count++;
            break;
         }

         if (count > std::max(col_pos["COD_UTS1"], col_pos["PROF_UTILE"]))
            break;
      }

      if (std::find(this->soil_id_.begin(), this->soil_id_.end(), cod) == this->soil_id_.end())
      {
         this->soil_id_.push_back(cod);
         auto max_z = prof / 100.0;

         std::vector<double> rand;

         for (unsigned int i = 0; i < this->num_landslides; i++)
         {
            rand.push_back(
               KiLib::Random::qtri(stats::runif(0, 1, engine), (2.0 / 3.0) * max_z, max_z, (3.0 / 4.0) * max_z));
         }
         this->z_.push_back(rand);
      }
   }

   fin.close();
   ss.clear();
   col_pos.clear();

   // ------------------------------
   // Reading Root Data
   // ------------------------------
   fin.open(root_data, std::ios::in);
   if (!fin.is_open())
   {
      spdlog::error("File '{}' failed to open", root_data);
      exit(EXIT_FAILURE);
   }

   count = 0;
   getline(fin, line);
   ss.str(line);

   while (ss.good() && getline(ss, word, ','))
   {
      // TODO: Find a cleaner way of doing this
      if (
         word.find("Pa400") || word.find("Pa200") || word.find("Fs800") || word.find("Fs200") || word.find("Cs150") ||
         word.find("MF600") || word.find("MF300"))
         col_pos[word] = count;
      count++;
   }

   while (getline(fin, line))
   {
      // getline strips the last '\n' which is necessary to check for data in the last column
      line += '\n';
      state = '\n';
      last  = 0;
      count = 0;

      for (size_t it = 0; it < line.size(); it++)
      {
         switch (line[it])
         {
         // Check if it is necessary to consider quataion
         case '"':
         case '\'':
            if (state == line[it])
               state = '\n';
            else if (state == '\n')
               state = line[it];
            break;
         case ',':
            if (state != '\n')
               break;
         case '\n':
            if (count == col_pos["Pa400"])
               this->Pa400.emplace_back(std::stod(line.substr(last, it - last)) * 1000);
            else if (count == col_pos["Pa200"])
               this->Pa200.emplace_back(std::stod(line.substr(last, it - last)) * 1000);
            else if (count == col_pos["Fs200"])
               this->Fs200.emplace_back(std::stod(line.substr(last, it - last)) * 1000);
            else if (count == col_pos["Fs800"])
               this->Fs800.emplace_back(std::stod(line.substr(last, it - last)) * 1000);
            else if (count == col_pos["MF300"])
               this->Mf300.emplace_back(std::stod(line.substr(last, it - last)) * 1000);
            else if (count == col_pos["MF600"])
               this->Mf600.emplace_back(std::stod(line.substr(last, it - last)) * 1000);
            else if (count == col_pos["Cs150"])
               this->Cs150.emplace_back(std::stod(line.substr(last, it - last)) * 1000);

            last = it + 1;
            count++;
            break;
         }

         // TODO: Find a way of terminating the loop early in a neat manner.
      }
   }

   fin.close();

   spdlog::info("Reading Dataset elapsed time: {}", sw);
}


/**
 * @brief Uses the data which has been read and generates soil properties based on that
 *
 */
void Primula::GenerateSoilProperties()
{
   spdlog::stopwatch sw;

   for (size_t i = 0; i < this->num_landslides; i++)
   {
      Landslide slide;

      // generate random soil properties

      this->phi1.push_back(stats::qunif(stats::runif(0, 1, this->engine), 30, 40));
      this->phi2.push_back(stats::qunif(stats::runif(0, 1, this->engine), 35, 40));
      this->gamma1.push_back(stats::qunif(stats::runif(0, 1, this->engine), 17, 19) * 1000);
      this->ks1.push_back(stats::qunif(stats::runif(0, 1, this->engine), 0.5, 100));
      this->ks2.push_back(stats::qunif(stats::runif(0, 1, this->engine), 0.5, 100));

      // generate random landslide properties
      slide.area_   = pow(10, stats::qnorm(stats::runif(0, 1, this->engine), this->area_mu_, this->area_sigma_));
      auto l2w      = pow(10, stats::qnorm(stats::runif(0, 1, this->engine), this->l2w_mu_, this->l2w_sigma_));
      slide.width_  = sqrt((slide.area_ * 1.0) / (l2w * 1.0));
      slide.length_ = slide.width_ * l2w;
      this->landslide_.push_back(slide);

      // pick random forest density
      // TODO: Replace with vectorized discrete uniform operation
      this->iteration_index.emplace_back(rand() % this->Pa200.size());

      this->Cr_grassland_.push_back(stats::qunif(stats::runif(0, 1, this->engine), 5, 7.5) * 1000);
      this->Cr_shrubland_.push_back(stats::qunif(stats::runif(0, 1, this->engine), 0, 15) * 1000);
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

   this->pr_failure_              = KiLib::Raster::zerosLike(this->probslope_);
   this->pr_failure_.nodata_value = -9999;

   for (unsigned int i = 0; i < num_landslides; i++)
   {
      KiLib::Raster friction_angle = KiLib::Raster::zerosLike(this->probslope_);
      KiLib::Raster permeability   = KiLib::Raster::zerosLike(this->probslope_);
      KiLib::Raster depth          = KiLib::Raster::zerosLike(this->probslope_);
      KiLib::Raster crl            = KiLib::Raster::zerosLike(this->probslope_);
      KiLib::Raster crb            = KiLib::Raster::zerosLike(this->probslope_);

      // go through each raster cell
      for (size_t j : this->validIndices)
      {
         // if soil 1 or 2, translate info to rasters
         if (this->soil_type_(j))
         {
            auto &phiv = (int)this->soil_type_(j) == 1 ? this->phi1 : this->phi2;
            auto &ksv  = (int)this->soil_type_(j) == 1 ? this->ks1 : this->ks2;
            // use the number to determine which element of the vector to access
            friction_angle(j) = phiv[i] * M_PI / 180.0;
            permeability(j)   = ksv[i] * M_PI / 180.0;
         }

         if (this->soil_depth_(j))
         {
            // add the depth of the soil id in the raster to another raster
            for (size_t k = 0; k < this->soil_id_.size(); k++)
            {
               if (this->soil_depth_(j) == this->soil_id_[k])
               {
                  depth(j) = z_[k][i];
                  break;
               }
            }
         }

         // copy dusaf raster, replacing codes with appropriate forest density
         switch ((int)dusaf_(j))
         {
         case 3211:
         case 3212:
         case 3221:
            crl(j) = Cr_grassland_[i];
            break;
         case 332:
         case 333:
            crl(j) = Cr_shrubland_[i];
            break;
         case 3121:
            crl(j) = this->Pa400[this->iteration_index[i]];
            break;
         case 3122:
            crl(j) = this->Pa200[this->iteration_index[i]];
            break;
         case 31111:
            crl(j) = this->Fs800[this->iteration_index[i]];
            break;
         case 31121:
            crl(j) = this->Fs200[this->iteration_index[i]];
            break;
         case 3114:
         case 222:
            crl(j) = this->Cs150[this->iteration_index[i]];
            break;
         case 31311:
            crl(j) = this->Mf600[this->iteration_index[i]];
            break;
         case 31321:
            crl(j) = this->Mf300[this->iteration_index[i]];
            break;
         default:
            crl(j) = 0;
            break;
         }

         if (depth(j) >= 0.5)
            crb(j) = 0;
         else
            crb(j) = crl(j);
      }

      auto m = CalcWetness(permeability, depth);

      auto FS = MDSTab_v2(this->landslide_[i], friction_angle, m, this->gamma1[i], depth, crl, crb);
      for (size_t j = 0; j < this->pr_failure_.nData; j++)
      {
         if (this->probslope_(j) == this->probslope_.nodata_value)
         {
            this->pr_failure_(j) = this->pr_failure_.nodata_value;
         }
         else if (FS(j) < 1 && FS(j) > 0)
         {
            this->pr_failure_(j) += FS(j);
         }
      }
   }

   // get average of sum of failure probabilities
   for (auto &c : this->pr_failure_.data)
   {
      if (c != this->pr_failure_.nodata_value)
         c /= num_landslides;
      else
         c = -9999;
   }

   spdlog::info("Landslide generation elapsed time: {}", sw);
}

void Primula::syncValidIndices()
{
   for (size_t i = 0; i < this->slope_.nData; i++) {
      if (this->slope_(i) == this->slope_.nodata_value)
         continue;
      if (this->twi_(i) == this->twi_.nodata_value)
         continue;
      if (this->soil_type_(i) == this->soil_type_.nodata_value)
         continue;
      if (this->soil_depth_(i) == this->soil_depth_.nodata_value)
         continue;
      if (this->dusaf_(i) == this->dusaf_.nodata_value)
         continue;
      if (this->probslope_(i) == this->probslope_.nodata_value)
         continue;

      this->validIndices.push_back(i);
   }

   spdlog::info("{} / {} Raster indices have valid data. Only computing on valid data.", this->validIndices.size(), this->slope_.nData);
}
