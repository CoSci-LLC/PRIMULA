// ============================================================================
// Copyright (C) PRIMULA++. All rights reserved.
//
// Authors: Denis Cohen-Corticchiato (DOHCC)
// Email:   denis.cohen@gmail.com
//
// This file is part of PRIMULA++  
// ============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <string>
#include <sstream>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/lognormal_distribution.hpp>
#include <boost/math/distributions/inverse_gamma.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/distributions/triangular.hpp>
#include <boost/math/distributions/uniform.hpp>
#include <boost/math/distributions/normal.hpp>

#include "primula++.hpp"
//#include "raster_functions.hpp"
#include "random_number_generator.hpp"
//#include "tree.hpp"

Primula::Primula()
{
}

Primula::~Primula()
{
}

bool Primula::ReadCSV(const std::string & file, const unsigned int & num_landslides)
{
// ------------------------------------------------------
// ... Uniform random number generator for landslides ...
// ------------------------------------------------------
   boost::mt19937 rng;                // Always same sequence for the moment
 //boost::mt19937 rng(std::time(0));  // Randomise generator
   static boost::uniform_01<boost::mt19937> rng_uniform_01(rng);

   std::cout << "Primula:ReadCSV \"" << file << "\"" << std::endl;

   // Open data file
   std::ifstream fin;
   fin.open(file);
   if (!fin.is_open())
   {
      std::cerr << "  File \"" + file + "\" open failed" << std::endl;
      return false;
   }

   std::string line;
   std::string word;
   int cod_pos, prof_pos, cod, prof;

   // find columns with COD_UTS1 and PROF_UTILE
   int count = 0;
   getline(fin, line);
   std::stringstream ss(line);
   while (ss.good())
   {
      getline(ss,word,',');
      word.erase(std::remove_if(word.begin(),word.end(),
         [](auto const & c) -> bool {return std::iscntrl(c);}),word.end());
      if (word == "COD_UTS1") cod_pos = count;
      else if (word == "PROF_UTILE") prof_pos = count;
      count++;
   }

   // extract COD_UTS1 and PROF_UTILE info in each line
   while (getline(fin, line))
   {
      count = 0;
      bool quote = false;
      std::vector<char> cstr(line.c_str(), line.c_str() + line.size() + 1);
      std::ostringstream out;
      for (char c: cstr){
         if (c == '"')
         {
            if (quote)
            {
               std::string s(out.str());
               if (s != "")
               {
                  if (count == cod_pos) std::stringstream(s) >> cod;
                  else if (count == prof_pos) std::stringstream(s) >> prof;
                  count++;
               }
               quote = false;
               out.str("");
               out.clear();
            } else quote = true;
         } else if ((c == ',' && !quote) || c == '\n')
         {
            std::string s(out.str());
            if (s != "")
            {
               if (count == cod_pos) std::stringstream(s) >> cod;
               else if (count == prof_pos) std::stringstream(s) >> prof;
               count++;
            }
            out.str("");
            out.clear();
         } else out << c;
      }

      // store if info not already found
      if (std::find(soil_id_.begin(), soil_id_.end(), cod)==soil_id_.end())
      {
         soil_id_.push_back(cod);
         auto max_z = prof/100.0;

         std::vector<double> rand;
         boost::math::triangular_distribution<> tri((2.0/3.0)*max_z,(3.0/4.0)*max_z,max_z);

         for (unsigned int i = 0; i < num_landslides; i++)
         {
            rand.push_back(quantile(tri, rng_uniform_01()));
         }
         z_.push_back(rand);
      }
   }

   fin.close();
   std::cout << "Primula:ReadCSV ... Done" << std::endl;
   return true;
}


Raster Primula::TopModel_v3(const Raster & ks, const Raster & z)
{
   Raster W(ks);

   for (size_t i = 0; i < ks.attribute_.size(); i++)
   {
      if (slope_.attribute_.at(i) != slope_.nodata_value_)
      {
         W.attribute_.at(i) = ((rainfall_ * 0.9) / (ks.attribute_.at(i) * z.attribute_.at(i) * cos(slope_.attribute_.at(i)))) * twi_.attribute_.at(i);
         if (W.attribute_.at(i) > 1) W.attribute_.at(i) = 1.0;
      } else
      {
         W.attribute_.at(i) = W.nodata_value_;
      }
   }

   return W;
}

Raster Primula::MDSTab_v2(const Landslide & slide, const Raster & phi, const Raster & m, const double & gamma_s, const Raster & z, const Raster & Crl, const Raster & Crb)
{
   auto gamma_w = 9810; // unit weight of water [kN/m3]
   Raster FS(m);

   // calculate factor of safety for each raster cell
   for (size_t i = 0; i < FS.attribute_.size(); i++)
   {
      auto tmp_phi = phi.attribute_.at(i);
      auto tmp_m = m.attribute_.at(i);
      auto tmp_z = z.attribute_.at(i);
      auto tmp_Crl = Crl.attribute_.at(i);
      auto tmp_Crb = Crb.attribute_.at(i);
      auto theta = probslope_.attribute_.at(i);
      auto delta = slope_.attribute_.at(i);

      // run calculation if probslope is non-empty
      if (theta != probslope_.nodata_value_)
      {
         auto Fdc = gamma_s * tmp_z * sin(theta) * cos(theta) * slide.width_ * slide.length_;
         if (Fdc)
         {
            auto K0 = 1.0 - sin(tmp_phi);
            // long equation derived from MDSTAB_v2.m
            auto Frl = 0.5 * K0 * (gamma_s - gamma_w * pow(tmp_m,2)) * slide.length_ * pow(tmp_z,2) * cos(theta) * tan(tmp_phi) + tmp_Crl * slide.length_ * tmp_z * cos(theta);
            
            // Rankine solution for cohesive soils
            // Used in MDSTAB_V2.m
            auto K = 4.0 * pow(cos(theta),2) * (pow(cos(theta),2) - pow(cos(tmp_phi),2)) + (4 * pow(tmp_Crl/(gamma_s * tmp_z),2) * pow(cos(tmp_phi),2)) + (8 * (tmp_Crl/(gamma_s * tmp_z)) * pow(cos(theta),2) * sin(tmp_phi) * cos(tmp_phi));
            if (K < 0) K = 0;

            auto Kp = (1 / pow(cos(tmp_phi),2)) * (2 * pow(cos(theta),2) + 2 * (tmp_Crl/(gamma_s * tmp_z)) * cos(tmp_phi) * sin(tmp_phi) + sqrt(K)) - 1;
            auto Ka = (1 / pow(cos(tmp_phi),2)) * (2 * pow(cos(theta),2) + 2 * (tmp_Crl/(gamma_s * tmp_z)) * cos(tmp_phi) * sin(tmp_phi) - sqrt(K)) - 1;

            // net driving force of the upslope margin
            auto Fdu = 0.5 * Ka * pow(tmp_z,2) * (gamma_s - gamma_w * pow(tmp_m,2)) * slide.width_ * cos(delta - theta);
            auto Fnu = 0.5 * Ka * pow(tmp_z,2) * (gamma_s - gamma_w * pow(tmp_m,2)) * slide.width_ * sin(delta - theta);
            //std::cout << Fdu << " " << Fnu << "\n";
            //std::cout << 0.5 * Ka * pow(tmp_z,2) << "\n";
            //std::cout << pow(tmp_m,2) << "\n";

            // Passive force on the downslope margin
            auto Frd = 0.5 * Kp * pow(tmp_z,2) * (gamma_s - gamma_w * pow(tmp_m,2)) * slide.width_ * cos(delta - theta);
            // Negligible, so clearly we set it to 0 immediately after calculating it ¯\_(ツ)_/¯
            Frd = 0;
            auto Fnd = 0.5 * Kp * pow(tmp_z,2) * (gamma_s - gamma_w * pow(tmp_m,2)) * slide.width_ * sin(delta - theta);

            // Basal resistance force
            auto Fnc = (gamma_s - gamma_w * tmp_m) * tmp_z * pow(cos(theta),2) * slide.width_ * slide.length_;
            auto Fnt = Fnc + Fnu - Fnd;
            auto Frb = tmp_Crb * slide.width_ * slide.length_ + Fnt * tan(tmp_phi);

            FS.attribute_.at(i) = (Frb + 2 * Frl + Frd - Fdu) / Fdc;
         }
      }
   }

   return FS;
}

bool Primula::GenerateLandslides(const std::string & file, const unsigned int & num_landslides)
{
   //landslide_.resize(num_landslides);

   //GenerateRandomUniform01(landlide_.x_, dem_.xllcorner_, dem_.yllcorner_, dem_xurcorner_, dem_.yurcorner_);

// ------------------------------------------------------
// ... Uniform random number generator for landslides ...
// ------------------------------------------------------
   boost::mt19937 rng;                // Always same sequence for the moment
 //boost::mt19937 rng(std::time(0));  // Randomise generator
   static boost::uniform_01<boost::mt19937> rng_uniform_01(rng);
   
// ------------------------------------------------------
// ... uniform distribution functions ...
// ------------------------------------------------------
   boost::math::uniform_distribution<> rng_phi1(30,40);
   boost::math::uniform_distribution<> rng_phi2(35,40);
   boost::math::uniform_distribution<> rng_gamma(17,19);
   boost::math::uniform_distribution<> rng_ks(0.5,100);
   boost::math::uniform_distribution<> rng_cr_grass(5,7.5);
   boost::math::uniform_distribution<> rng_cr_shrub(0,15);

// --------------------------------------------------------
// ... normal distribution for landslide area ...
// --------------------------------------------------------
   boost::math::normal_distribution<> rng_area(area_mu_, area_sigma_);

// --------------------------------------------------------
// ... normal distribution for length-to-width ratio ...
// --------------------------------------------------------
   boost::math::normal_distribution<> rng_l2w(l2w_mu_, l2w_sigma_);

// ----------------------------------------------
// ... Generate soil properties ...
// ... Generate random data ...
// ----------------------------------------------
   auto start_rng = std::chrono::high_resolution_clock::now(); // Time loop
   
   // Open data file
   std::ifstream fin;
   fin.open(file);
   if (!fin.is_open())
   {
      std::cerr << "  File \"" + file + "\" open failed" << std::endl;
      return false;
   }

   unsigned int count = 0;
   unsigned int Fs200_pos, Fs800_pos, Pa200_pos, Pa400_pos, Mf300_pos, Mf600_pos, Cs150_pos;
   std::vector<double> Fs200, Fs800, Pa200, Pa400, Mf600, Mf300, Cs150;
   std::string line, word;

   // get column numbers for each data set
   getline(fin, line);
   std::stringstream ss(line);
   while (ss.good())
   {
      getline(ss,word,',');
      word.erase(std::remove_if(word.begin(),word.end(),
         [](auto const & c) -> bool {return std::iscntrl(c);}),word.end());
      
      if (word == "Pa400") Pa400_pos = count;
      else if (word == "Pa200") Pa200_pos = count;
      else if (word == "Fs800") Fs800_pos = count;
      else if (word == "Fs200") Fs200_pos = count;
      else if (word == "Cs150") Cs150_pos = count;
      else if (word == "MF600") Mf600_pos = count;
      else if (word == "MF300") Mf300_pos = count;
      count++;
   }

   // store data in each line
   while (getline(fin, line))
   {
      count = 0;
      std::stringstream ss(line);
      while (ss.good())
      {
         getline(ss,word,',');
         if (count == Pa400_pos) Pa400.push_back(std::stod(word)*1000);
         else if (count == Pa200_pos) Pa200.push_back(std::stod(word)*1000);
         else if (count == Fs800_pos) Fs800.push_back(std::stod(word)*1000);
         else if (count == Fs200_pos) Fs200.push_back(std::stod(word)*1000);
         else if (count == Cs150_pos) Cs150.push_back(std::stod(word)*1000);
         else if (count == Mf600_pos) Mf600.push_back(std::stod(word)*1000);
         else if (count == Mf300_pos) Mf300.push_back(std::stod(word)*1000);
         count++;
      }
   }

   std::vector<std::vector<double>> phi;  // soil friction angle (rad)
   std::vector<std::vector<double>> gamma;  // specific weight falues [N/m^3]
   std::vector<std::vector<double>> ks;  // soil permeability [m/day]

   std::vector<double> phi1;  // phi values for soil 1
   std::vector<double> phi2;  // phi values for soil 2
   std::vector<double> gamma1;  // gamma values (for soil 1?)
   std::vector<double> ks1;  // ks values for soil 1
   std::vector<double> ks2;  // ks values for soil 2
   for (unsigned int i = 0; i < num_landslides; i++)
   {
      Landslide slide;

      // generate random soil properties
      phi1.push_back(quantile(rng_phi1, rng_uniform_01()));
      phi2.push_back(quantile(rng_phi2, rng_uniform_01()));
      gamma1.push_back(quantile(rng_gamma, rng_uniform_01())*1000);
      ks1.push_back(quantile(rng_ks, rng_uniform_01()));
      ks2.push_back(quantile(rng_ks, rng_uniform_01()));

      // generate random landslide properties
      slide.area_ = pow(10,quantile(rng_area, rng_uniform_01()));
      auto l2w = pow(10,quantile(rng_l2w, rng_uniform_01()));
      slide.width_ = sqrt((slide.area_ * 1.0) / (l2w * 1.0));
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

      Cr_grassland_.push_back(quantile(rng_cr_grass,rng_uniform_01())*1000);
      Cr_shrubland_.push_back(quantile(rng_cr_shrub,rng_uniform_01())*1000);
   }
   // add to vector for easier access
   phi.push_back(phi1);
   phi.push_back(phi2);
   gamma.push_back(gamma1);
   ks.push_back(ks1);
   ks.push_back(ks2);

   auto finish_rng = std::chrono::high_resolution_clock::now(); // End time bi-linear interpolation
   std::chrono::duration<double> elapsed_rng = finish_rng - start_rng;
   std::cout << "Soil generation elapsed time: " << elapsed_rng.count() << " s\n";

// ----------------------------------------------
// ... Landslide generation ...
// ----------------------------------------------
   auto start_sli = std::chrono::high_resolution_clock::now(); // Time bi-linear interpolation

   Raster Pr_failure(probslope_);
   Pr_failure.nodata_value_ = -9999;

   for (unsigned int i = 0; i < num_landslides; i++)
   {
      Raster friction_angle(soil_type_);
      Raster permeability(soil_type_);
      Raster depth(soil_type_);
      Raster crl(soil_type_);
      Raster crb(soil_type_);

      // go through each raster cell
      for (size_t j = 0; j < soil_type_.attribute_.size(); j++)
      {
         if (probslope_.attribute_.at(j) != probslope_.nodata_value_) {
            // if soil 1 or 2, translate info to rasters
            if (soil_type_.attribute_.at(j))
            {
               // use the number to determine which element of the vector to access
               friction_angle.attribute_.at(j) = phi.at((int)soil_type_.attribute_.at(j)-1).at(i) * M_PI/180.0;
               permeability.attribute_.at(j) = ks.at((int)soil_type_.attribute_.at(j)-1).at(i) * M_PI/180.0;
            }

            if (soil_depth_.attribute_.at(j))
            {
               // add the depth of the soil id in the raster to another raster
               for (size_t k = 0; k < soil_id_.size(); k++)
               {
                  if (soil_depth_.attribute_.at(j) == soil_id_.at(k))
                  {
                     depth.attribute_.at(j) = z_.at(k).at(i);
                     break;
                  }
               }
            }

            // copy dusaf raster, replacing codes with appropriate forest density
            switch ((int)dusaf_.attribute_.at(j))
            {
               case 3211:
               case 3212:
               case 3221:
                  crl.attribute_.at(j) = Cr_grassland_.at(i);
                  break;
               case 332:
               case 333:
                  crl.attribute_.at(j) = Cr_shrubland_.at(i);
                  break;
               case 3121:
                  crl.attribute_.at(j) = Crl_Pa400_.at(i);
                  break;
               case 3122:
                  crl.attribute_.at(j) = Crl_Pa200_.at(i);
                  break;
               case 31111:
                  crl.attribute_.at(j) = Crl_Fs800_.at(i);
                  break;
               case 31121:
                  crl.attribute_.at(j) = Crl_Fs200_.at(i);
                  break;
               case 3114:
               case 222:
                  crl.attribute_.at(j) = Crl_Cs150_.at(i);
                  break;
               case 31311:
                  crl.attribute_.at(j) = Crl_Mf600_.at(i);
                  break;
               case 31321:
                  crl.attribute_.at(j) = Crl_Mf300_.at(i);
                  break;
               default:
                  crl.attribute_.at(j) = 0;
                  break;
            }

            if (depth.attribute_.at(j) >= 0.5) crb.attribute_.at(j) = 0;
            else crb.attribute_.at(j) = crl.attribute_.at(j);
         }
      }

      auto m = TopModel_v3(permeability,depth);

      auto FS = MDSTab_v2(landslide_.at(i), friction_angle, m, gamma.at(0).at(i), depth, crl, crb);
      for (size_t j = 0; j < Pr_failure.attribute_.size(); j++)
      {
         if (probslope_.attribute_.at(j) == probslope_.nodata_value_) Pr_failure.attribute_.at(j) = Pr_failure.nodata_value_;
         else if (FS.attribute_.at(j) < 1 && FS.attribute_.at(j) > 0) {
            Pr_failure.attribute_.at(j) += FS.attribute_.at(j);
         }
      }
   }

   // get average of sum of failure probabilities
   for (auto & c : Pr_failure.attribute_)
   {
      if (c != Pr_failure.nodata_value_) c /= num_landslides;
      else c = -9999;
   }

   pr_failure_ = Pr_failure;
   pr_failure_.nodata_value_ = -9999;

   auto finish_sli = std::chrono::high_resolution_clock::now(); // End time bi-linear interpolation
   std::chrono::duration<double> elapsed_sli = finish_sli - start_sli;
   std::cout << "Landslide generation elapsed time: " << elapsed_sli.count() << " s\n";

   return true;
}

