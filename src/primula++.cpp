// ============================================================================
// Copyright (C) PRIMULA++. All rights reserved.
//
// Authors: Denis Cohen-Corticchiato (DOHCC)
// Email:   denis.cohen@gmail.com
//
// This file is part of PRIMULA++
// ============================================================================

#define STATS_USE_OPENMP

#include "primula++.hpp"
#include <chrono>
#include <random>
#include <spdlog/spdlog.h>
#include <spdlog/stopwatch.h>
#include <stats.hpp>
#include <string>
#include <unordered_map>

std::mt19937_64 engine(69420); // Engine so our seed is consistent

/**
 * @brief Returns the quantile of the given value in a triangular distribution
 *
 * @param p The probability
 * @param a The lower bound
 * @param b The upper bound
 * @param c The mode
 * @return double The quantile
 */
static double qtri(const double &p, const double &a, const double &b, const double &c)
{
   if (p < c)
      return a + std::sqrt((b - a) * (c - a) * p);
   else if (p > c)
      return b - std::sqrt((b - a) * (b - c) * (1 - p));

   return c;
}

Primula::~Primula()
{
}

KiLib::Raster Primula::TopModel_v3(const KiLib::Raster &ks, const KiLib::Raster &z)
{
   KiLib::Raster W(ks);

   for (size_t i = 0; i < ks.data.size(); i++) {
      if (slope_.data.at(i) != slope_.nodata_value) {
         W.data.at(i) = ((rainfall_ * 0.9) / (ks.data.at(i) * z.data.at(i) * cos(slope_.data.at(i)))) * twi_.data.at(i);
         if (W.data.at(i) > 1)
            W.data.at(i) = 1.0;
      } else {
         W.data.at(i) = W.nodata_value;
      }
   }

   return W;
}

KiLib::Raster Primula::MDSTab_v2(
   const Landslide &slide, const KiLib::Raster &phi, const KiLib::Raster &m, const double &gamma_s,
   const KiLib::Raster &z, const KiLib::Raster &Crl, const KiLib::Raster &Crb)
{
   auto gamma_w = 9810; // unit weight of water [kN/m3]
   KiLib::Raster FS(m);

   // calculate factor of safety for each raster cell
   for (size_t i = 0; i < FS.data.size(); i++) {
      auto tmp_phi = phi.data.at(i);
      auto tmp_m   = m.data.at(i);
      auto tmp_z   = z.data.at(i);
      auto tmp_Crl = Crl.data.at(i);
      auto tmp_Crb = Crb.data.at(i);
      auto theta   = this->probslope_.data.at(i);
      auto delta   = slope_.data.at(i);

      // run calculation if probslope is non-empty
      if (theta != this->probslope_.nodata_value) {
         auto Fdc = gamma_s * tmp_z * sin(theta) * cos(theta) * slide.width_ * slide.length_;
         if (Fdc) {
            auto K0 = 1.0 - sin(tmp_phi);
            // long equation derived from MDSTAB_v2.m
            auto Frl = 0.5 * K0 * (gamma_s - gamma_w * pow(tmp_m, 2)) * slide.length_ * pow(tmp_z, 2) * cos(theta) *
                          tan(tmp_phi) +
                       tmp_Crl * slide.length_ * tmp_z * cos(theta);

            // Rankine solution for cohesive soils
            // Used in MDSTAB_V2.m
            auto K = 4.0 * pow(cos(theta), 2) * (pow(cos(theta), 2) - pow(cos(tmp_phi), 2)) +
                     (4 * pow(tmp_Crl / (gamma_s * tmp_z), 2) * pow(cos(tmp_phi), 2)) +
                     (8 * (tmp_Crl / (gamma_s * tmp_z)) * pow(cos(theta), 2) * sin(tmp_phi) * cos(tmp_phi));
            if (K < 0)
               K = 0;

            auto Kp =
               (1 / pow(cos(tmp_phi), 2)) *
                  (2 * pow(cos(theta), 2) + 2 * (tmp_Crl / (gamma_s * tmp_z)) * cos(tmp_phi) * sin(tmp_phi) + sqrt(K)) -
               1;
            auto Ka =
               (1 / pow(cos(tmp_phi), 2)) *
                  (2 * pow(cos(theta), 2) + 2 * (tmp_Crl / (gamma_s * tmp_z)) * cos(tmp_phi) * sin(tmp_phi) - sqrt(K)) -
               1;

            // net driving force of the upslope margin
            auto Fdu =
               0.5 * Ka * pow(tmp_z, 2) * (gamma_s - gamma_w * pow(tmp_m, 2)) * slide.width_ * cos(delta - theta);
            auto Fnu =
               0.5 * Ka * pow(tmp_z, 2) * (gamma_s - gamma_w * pow(tmp_m, 2)) * slide.width_ * sin(delta - theta);

            // Passive force on the downslope margin
            auto Frd =
               0.5 * Kp * pow(tmp_z, 2) * (gamma_s - gamma_w * pow(tmp_m, 2)) * slide.width_ * cos(delta - theta);
            // Negligible, so clearly we set it to 0 immediately after calculating it ¯\_(ツ)_/¯
            Frd = 0;
            auto Fnd =
               0.5 * Kp * pow(tmp_z, 2) * (gamma_s - gamma_w * pow(tmp_m, 2)) * slide.width_ * sin(delta - theta);

            // Basal resistance force
            auto Fnc = (gamma_s - gamma_w * tmp_m) * tmp_z * pow(cos(theta), 2) * slide.width_ * slide.length_;
            auto Fnt = Fnc + Fnu - Fnd;
            auto Frb = tmp_Crb * slide.width_ * slide.length_ + Fnt * tan(tmp_phi);

            FS.data.at(i) = (Frb + 2 * Frl + Frd - Fdu) / Fdc;
         }
      }
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

   size_t last  = 0;
   size_t count = 0;
   int cod      = -1;
   int prof     = -1;
   char state   = '\n';

   std::unordered_map<std::string, size_t> col_pos;
   std::stringstream ss;
   std::string line;
   std::string word;

   std::ifstream fin;

   spdlog::info("Primula::ReadSoilDataset Start");

   // TODO: Consider replacing this with a CSV parsing library

   // ------------------------------
   // Reading Soil Data
   // ------------------------------
   fin.open(soil_data, std::ios::in);
   if (!fin.is_open()) {
      spdlog::error("File '{}' failed to open", soil_data);
      exit(EXIT_FAILURE);
   }

   count = 0;
   getline(fin, line);
   ss.str(line);

   while (ss.good() && getline(ss, word, ',')) {
      col_pos[word] = count;
      count++;
   }

   while (getline(fin, line)) {
      // getline strips the last '\n' which is necessary to check for data in the last column
      line += '\n';
      state = '\n';
      last  = 0;
      count = 0;
      cod   = -1;
      prof  = -1;

      for (size_t it = 0; it < line.size(); it++) {

         switch (line[it]) {
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
            last = it + 1;
            count++;
            break;
         }

         if (count > std::max(col_pos["COD_UTS1"], col_pos["PROF_UTILE"]))
            break;
      }

      if (std::find(this->soil_id_.begin(), this->soil_id_.end(), cod) == this->soil_id_.end()) {
         this->soil_id_.push_back(cod);
         auto max_z = prof / 100.0;

         std::vector<double> rand;

         for (unsigned int i = 0; i < this->num_landslides; i++) {
            rand.push_back(qtri(stats::runif(0, 1, engine), (2.0 / 3.0) * max_z, max_z, (3.0 / 4.0) * max_z));
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
   if (!fin.is_open()) {
      spdlog::error("File '{}' failed to open", root_data);
      exit(EXIT_FAILURE);
   }

   count = 0;
   getline(fin, line);
   ss.str(line);

   while (ss.good() && getline(ss, word, ',')) {
      // TODO: Find a cleaner way of doing this
      if (
         word.find("Pa400") || word.find("Pa200") || word.find("Fs800") || word.find("Fs200") || word.find("Cs150") ||
         word.find("MF600") || word.find("MF300"))
         col_pos[word] = count;
      count++;
   }

   while (getline(fin, line)) {
      // getline strips the last '\n' which is necessary to check for data in the last column
      line += '\n';
      state = '\n';
      last  = 0;
      count = 0;

      for (size_t it = 0; it < line.size(); it++) {
         switch (line[it]) {
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

   spdlog::info("Primula::ReadSoilDataset End");
}


/**
 * @brief Uses the data which has been read and generates soil properties based on that
 *
 */
void Primula::GenerateSoilProperties()
{
   spdlog::stopwatch sw;

   for (size_t i = 0; i < this->num_landslides; i++) {
      Landslide slide;

      // generate random soil properties

      this->phi1.push_back(stats::qunif(stats::runif(0, 1, engine), 30, 40));
      this->phi2.push_back(stats::qunif(stats::runif(0, 1, engine), 35, 40));
      this->gamma1.push_back(stats::qunif(stats::runif(0, 1, engine), 17, 19) * 1000);
      this->ks1.push_back(stats::qunif(stats::runif(0, 1, engine), 0.5, 100));
      this->ks2.push_back(stats::qunif(stats::runif(0, 1, engine), 0.5, 100));

      // generate random landslide properties
      slide.area_   = pow(10, stats::qnorm(stats::runif(0, 1, engine), this->area_mu_, this->area_sigma_));
      auto l2w      = pow(10, stats::qnorm(stats::runif(0, 1, engine), this->l2w_mu_, this->l2w_sigma_));
      slide.width_  = sqrt((slide.area_ * 1.0) / (l2w * 1.0));
      slide.length_ = slide.width_ * l2w;
      this->landslide_.push_back(slide);

      // pick random forest density
      // TODO: Replace with vectorized discrete uniform operation
      this->iteration_index.emplace_back(rand() % this->Pa200.size());

      this->Cr_grassland_.push_back(stats::qunif(stats::runif(0, 1, engine), 5, 7.5) * 1000);
      this->Cr_shrubland_.push_back(stats::qunif(stats::runif(0, 1, engine), 0, 15) * 1000);
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

   this->pr_failure_ = KiLib::Raster::zerosLike(this->probslope_);

   for (size_t i = 0; i < this->num_landslides; i++) {
      KiLib::Raster friction_angle(this->soil_type_);
      KiLib::Raster permeability(this->soil_type_);
      KiLib::Raster depth(this->soil_type_);
      KiLib::Raster crl(this->soil_type_);
      KiLib::Raster crb(this->soil_type_);

      // go through each raster cell
      for (size_t j = 0; j < this->soil_type_.data.size(); j++) {
         if (this->probslope_.data.at(j) != this->probslope_.nodata_value) {
            // if soil 1 or 2, translate info to rasters
            if (this->soil_type_.data.at(j)) {
               auto &phiv = ((int)this->soil_type_.data.at(j) - 1) == 0 ? phi1 : phi2;
               auto &ksv  = ((int)this->soil_type_.data.at(j) - 1) == 0 ? ks1 : ks2;
               // use the number to determine which element of the vector to access
               friction_angle.data.at(j) = phiv.at(i) * M_PI / 180.0;
               permeability.data.at(j)   = ksv.at(i) * M_PI / 180.0;
            }

            if (this->soil_depth_.data.at(j)) {
               // add the depth of the soil id in the raster to another raster
               for (size_t k = 0; k < this->soil_id_.size(); k++) {
                  if (this->soil_depth_.data.at(j) == this->soil_id_.at(k)) {
                     depth.data.at(j) = z_.at(k).at(i);
                     break;
                  }
               }
            }

            // copy dusaf raster, replacing codes with appropriate forest density
            switch ((int)dusaf_.data.at(j)) {
            case 3211:
            case 3212:
            case 3221:
               crl.data.at(j) = Cr_grassland_.at(i);
               break;
            case 332:
            case 333:
               crl.data.at(j) = Cr_shrubland_.at(i);
               break;
            case 3121:
               crl.data.at(j) = this->Pa400[this->iteration_index[i]];//Crl_Pa400_.at(i);
               break;
            case 3122:
               crl.data.at(j) = this->Pa200[this->iteration_index[i]];//Crl_Pa200_.at(i);
               break;
            case 31111:
               crl.data.at(j) = this->Fs800[this->iteration_index[i]];//Crl_Fs800_.at(i);
               break;
            case 31121:
               crl.data.at(j) = this->Fs200[this->iteration_index[i]];//Crl_Fs200_.at(i);
               break;
            case 3114:
            case 222:
               crl.data.at(j) = this->Cs150[this->iteration_index[i]];//Crl_Cs150_.at(i);
               break;
            case 31311:
               crl.data.at(j) = this->Mf600[this->iteration_index[i]];//Crl_Mf600_.at(i);
               break;
            case 31321:
               crl.data.at(j) = this->Mf300[this->iteration_index[i]];//Crl_Mf300_.at(i);
               break;
            default:
               crl.data.at(j) = 0;
               break;
            }

            if (depth.data.at(j) >= 0.5)
               crb.data.at(j) = 0;
            else
               crb.data.at(j) = crl.data.at(j);
         }
      }

      auto m = TopModel_v3(permeability, depth);

      auto FS = MDSTab_v2(this->landslide_.at(i), friction_angle, m, this->gamma1.at(i), depth, crl, crb);
      for (size_t j = 0; j < this->pr_failure_.data.size(); j++) {
         if (FS.data.at(j) < 1 && FS.data.at(j) > 0)
            this->pr_failure_.data.at(j) += FS.data.at(j);
      }
   }

   // get average of sum of failure probabilities
   for (size_t i = 0; i < this->pr_failure_.data.size(); i++)
      if (this->pr_failure_.data[i] != this->pr_failure_.nodata_value)
         this->pr_failure_.data[i] /= this->num_landslides;

   spdlog::info("Landslide generation elapsed time: {}", sw);
}
