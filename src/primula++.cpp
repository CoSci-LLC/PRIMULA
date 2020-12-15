// ============================================================================
// Copyright (C) PRIMULA++. All rights reserved.
//
// Authors: Denis Cohen-Corticchiato (DOHCC)
// Email:   denis.cohen@gmail.com
//
// This file is part of PRIMULA++
// ============================================================================

#define STATS_USE_OPENMP

#include <chrono>
#include <random>
#include <stats.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/stopwatch.h>
#include "primula++.hpp"

/**
 * @brief Returns the quantile of the given value in a triangular distribution
 * 
 * @param p The probability
 * @param a The lower bound
 * @param b The upper bound
 * @param c The mode
 * @return double The quantile
 */
static double qtri(const double& p, const double& a, const double& b, const double& c)
{
   if (p < c)
      return a + std::sqrt((b-a)*(c-a)*p);
   else if (p > c)
      return b - std::sqrt((b-a)*(b-c)*(1-p));
   
   return c;
}

Primula::Primula()
{
}

Primula::~Primula()
{
}

void Primula::ReadCSV(const std::string &file, const unsigned int &num_landslides)
{
   spdlog::info("Primula::ReadCSV '{}' Start", file);

   // Open data file
   std::ifstream fin;
   fin.open(file);
   if (!fin.is_open()) {
      spdlog::error("File '{}' failed to open", file);
      exit(EXIT_FAILURE);
   }

   std::string line;
   std::string word;
   int cod_pos, prof_pos, cod, prof;

   // find columns with COD_UTS1 and PROF_UTILE
   int count = 0;
   getline(fin, line);
   std::stringstream ss(line);
   while (ss.good()) {
      getline(ss, word, ',');
      word.erase(
         std::remove_if(word.begin(), word.end(), [](auto const &c) -> bool { return std::iscntrl(c); }), word.end());
      if (word == "COD_UTS1")
         cod_pos = count;
      else if (word == "PROF_UTILE")
         prof_pos = count;
      count++;
   }

   // extract COD_UTS1 and PROF_UTILE info in each line
   while (getline(fin, line)) {
      count      = 0;
      bool quote = false;
      std::vector<char> cstr(line.c_str(), line.c_str() + line.size() + 1);
      std::ostringstream out;
      for (char c : cstr) {
         if (c == '"') {
            if (quote) {
               std::string s(out.str());
               if (s != "") {
                  if (count == cod_pos)
                     std::stringstream(s) >> cod;
                  else if (count == prof_pos)
                     std::stringstream(s) >> prof;
                  count++;
               }
               quote = false;
               out.str("");
               out.clear();
            } else
               quote = true;
         } else if ((c == ',' && !quote) || c == '\n') {
            std::string s(out.str());
            if (s != "") {
               if (count == cod_pos)
                  std::stringstream(s) >> cod;
               else if (count == prof_pos)
                  std::stringstream(s) >> prof;
               count++;
            }
            out.str("");
            out.clear();
         } else
            out << c;
      }

      // store if info not already found
      if (std::find(soil_id_.begin(), soil_id_.end(), cod) == soil_id_.end()) {
         soil_id_.push_back(cod);
         auto max_z = prof / 100.0;

         std::vector<double> rand;

         for (unsigned int i = 0; i < num_landslides; i++) {
            rand.push_back(qtri(stats::runif(0, 1, engine), (2.0 / 3.0) * max_z, max_z, (3.0 / 4.0) * max_z));
         }
         z_.push_back(rand);
      }
   }

   fin.close();
   spdlog::info("Primula::ReadCSV '{}' End", file);
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
      auto theta   = probslope_.data.at(i);
      auto delta   = slope_.data.at(i);

      // run calculation if probslope is non-empty
      if (theta != probslope_.nodata_value) {
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

void Primula::GenerateLandslides(const std::string &file, const unsigned int &num_landslides)
{
   // ----------------------------------------------
   // ... Generate soil properties ...
   // ... Generate random data ...
   // ----------------------------------------------
   std::mt19937_64 engine(69420); // Engine so our seed is consistent
   spdlog::stopwatch sw;

   // Open data file
   std::ifstream fin;
   fin.open(file);
   if (!fin.is_open()) {
      spdlog::error("File '{}' failed to open", file);
      exit(EXIT_FAILURE);
   }

   unsigned int count = 0;
   unsigned int Fs200_pos, Fs800_pos, Pa200_pos, Pa400_pos, Mf300_pos, Mf600_pos, Cs150_pos;
   std::vector<double> Fs200, Fs800, Pa200, Pa400, Mf600, Mf300, Cs150;
   std::string line, word;

   // get column numbers for each data set
   getline(fin, line);
   std::stringstream ss(line);
   while (ss.good()) {
      getline(ss, word, ',');
      word.erase(
         std::remove_if(word.begin(), word.end(), [](auto const &c) -> bool { return std::iscntrl(c); }), word.end());

      if (word == "Pa400")
         Pa400_pos = count;
      else if (word == "Pa200")
         Pa200_pos = count;
      else if (word == "Fs800")
         Fs800_pos = count;
      else if (word == "Fs200")
         Fs200_pos = count;
      else if (word == "Cs150")
         Cs150_pos = count;
      else if (word == "MF600")
         Mf600_pos = count;
      else if (word == "MF300")
         Mf300_pos = count;
      count++;
   }

   // store data in each line
   while (getline(fin, line)) {
      count = 0;
      std::stringstream ss(line);
      while (ss.good()) {
         getline(ss, word, ',');
         if (count == Pa400_pos)
            Pa400.push_back(std::stod(word) * 1000);
         else if (count == Pa200_pos)
            Pa200.push_back(std::stod(word) * 1000);
         else if (count == Fs800_pos)
            Fs800.push_back(std::stod(word) * 1000);
         else if (count == Fs200_pos)
            Fs200.push_back(std::stod(word) * 1000);
         else if (count == Cs150_pos)
            Cs150.push_back(std::stod(word) * 1000);
         else if (count == Mf600_pos)
            Mf600.push_back(std::stod(word) * 1000);
         else if (count == Mf300_pos)
            Mf300.push_back(std::stod(word) * 1000);
         count++;
      }
   }

   std::vector<std::vector<double>> phi;   // soil friction angle (rad)
   std::vector<std::vector<double>> gamma; // specific weight falues [N/m^3]
   std::vector<std::vector<double>> ks;    // soil permeability [m/day]

   std::vector<double> phi1;   // phi values for soil 1
   std::vector<double> phi2;   // phi values for soil 2
   std::vector<double> gamma1; // gamma values (for soil 1?)
   std::vector<double> ks1;    // ks values for soil 1
   std::vector<double> ks2;    // ks values for soil 2
   for (unsigned int i = 0; i < num_landslides; i++) {
      Landslide slide;

      // generate random soil properties

      phi1.push_back(stats::qunif(stats::runif(0, 1, engine), 30, 40));
      phi2.push_back(stats::qunif(stats::runif(0, 1, engine), 35, 40));
      gamma1.push_back(stats::qunif(stats::runif(0, 1, engine), 17, 19) * 1000);
      ks1.push_back(stats::qunif(stats::runif(0, 1, engine), 0.5, 100));
      ks2.push_back(stats::qunif(stats::runif(0, 1, engine), 0.5, 100));

      // generate random landslide properties
      slide.area_   = pow(10, stats::qnorm(stats::runif(0, 1, engine), area_mu_, area_sigma_));
      auto l2w      = pow(10, stats::qnorm(stats::runif(0, 1, engine), l2w_mu_, l2w_sigma_));
      slide.width_  = sqrt((slide.area_ * 1.0) / (l2w * 1.0));
      slide.length_ = slide.width_ * l2w;
      landslide_.push_back(slide);

      // pick random forest density
      auto n = rand() % Pa200.size();
      Crl_Fs200_.push_back(Fs200.at(n));
      Crl_Fs800_.push_back(Fs800.at(n));
      Crl_Pa200_.push_back(Pa200.at(n));
      Crl_Pa400_.push_back(Pa400.at(n));
      Crl_Mf300_.push_back(Mf300.at(n));
      Crl_Mf600_.push_back(Mf600.at(n));
      Crl_Cs150_.push_back(Cs150.at(n));

      Cr_grassland_.push_back(stats::qunif(stats::runif(0, 1, engine), 5, 7.5) * 1000);
      Cr_shrubland_.push_back(stats::qunif(stats::runif(0, 1, engine), 0, 15) * 1000);
   }
   // add to vector for easier access
   phi.push_back(phi1);
   phi.push_back(phi2);
   gamma.push_back(gamma1);
   ks.push_back(ks1);
   ks.push_back(ks2);

   spdlog::info("Soil generation elapsed time: {}", sw);
   sw.reset();

   // ----------------------------------------------
   // ... Landslide generation ...
   // ----------------------------------------------
   KiLib::Raster Pr_failure(probslope_);
   Pr_failure.nodata_value = -9999;

   for (unsigned int i = 0; i < num_landslides; i++) {
      KiLib::Raster friction_angle(soil_type_);
      KiLib::Raster permeability(soil_type_);
      KiLib::Raster depth(soil_type_);
      KiLib::Raster crl(soil_type_);
      KiLib::Raster crb(soil_type_);

      // go through each raster cell
      for (size_t j = 0; j < soil_type_.data.size(); j++) {
         if (probslope_.data.at(j) != probslope_.nodata_value) {
            // if soil 1 or 2, translate info to rasters
            if (soil_type_.data.at(j)) {
               // use the number to determine which element of the vector to access
               friction_angle.data.at(j) = phi.at((int)soil_type_.data.at(j) - 1).at(i) * M_PI / 180.0;
               permeability.data.at(j)   = ks.at((int)soil_type_.data.at(j) - 1).at(i) * M_PI / 180.0;
            }

            if (soil_depth_.data.at(j)) {
               // add the depth of the soil id in the raster to another raster
               for (size_t k = 0; k < soil_id_.size(); k++) {
                  if (soil_depth_.data.at(j) == soil_id_.at(k)) {
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
               crl.data.at(j) = Crl_Pa400_.at(i);
               break;
            case 3122:
               crl.data.at(j) = Crl_Pa200_.at(i);
               break;
            case 31111:
               crl.data.at(j) = Crl_Fs800_.at(i);
               break;
            case 31121:
               crl.data.at(j) = Crl_Fs200_.at(i);
               break;
            case 3114:
            case 222:
               crl.data.at(j) = Crl_Cs150_.at(i);
               break;
            case 31311:
               crl.data.at(j) = Crl_Mf600_.at(i);
               break;
            case 31321:
               crl.data.at(j) = Crl_Mf300_.at(i);
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

      auto FS = MDSTab_v2(landslide_.at(i), friction_angle, m, gamma.at(0).at(i), depth, crl, crb);
      for (size_t j = 0; j < Pr_failure.data.size(); j++) {
         if (probslope_.data.at(j) == probslope_.nodata_value)
            Pr_failure.data.at(j) = Pr_failure.nodata_value;
         else if (FS.data.at(j) < 1 && FS.data.at(j) > 0) {
            Pr_failure.data.at(j) += FS.data.at(j);
         }
      }
   }

   // get average of sum of failure probabilities
   for (auto &c : Pr_failure.data) {
      if (c != Pr_failure.nodata_value)
         c /= num_landslides;
      else
         c = -9999;
   }

   pr_failure_              = Pr_failure;
   pr_failure_.nodata_value = -9999;

   spdlog::info("Landslide generation elapsed time: {}", sw);
}
