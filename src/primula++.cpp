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

   std::vector<double> max_z_;
   std::string line;
   std::string word;
   int cod_pos, prof_pos, cod, prof;

   int count = 0;
   getline(fin, line);
   std::stringstream ss(line);
   while (ss.good())
   {
      getline(ss,word,',');
      if (word == "COD_UTS1") cod_pos = count;
      else if (word == "PROF_UTILE") prof_pos = count;
      count++;
   }

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

      if (std::find(soil_id_.begin(), soil_id_.end(), cod)==soil_id_.end())
      {
         soil_id_.push_back(cod);
         double max_z = prof/100.0;
         max_z_.push_back(max_z);

         std::vector<double> rand;
         boost::math::triangular_distribution<> tri(2/3*max_z,3/4*max_z,max_z);

         for (auto i = 0; i < num_landslides; i++)
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

   for (auto i = 0; i < ks.attribute_.size(); i++)
   {
      if (slope_.attribute_.at(i) != slope_.nodata_value_)
      {
         W.attribute_.at(i) = rainfall_ * 0.9 / (ks.attribute_.at(i) * z.attribute_.at(i) * cos(slope_.attribute_.at(i))) * twi_.attribute_.at(i);
         if (W.attribute_.at(i) > 1) W.attribute_.at(i) = 1;
      } else
      {
         W.attribute_.at(i) = slope_.nodata_value_;
      }
   }

   return W;
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
   boost::math::uniform_distribution<> rng_phi1(30,40);
   boost::math::uniform_distribution<> rng_phi2(35,40);
   boost::math::uniform_distribution<> rng_gamma(17,19);
   boost::math::uniform_distribution<> rng_ks(0.5,100);
   boost::math::uniform_distribution<> rng_cr_grass(5,7.5);
   boost::math::uniform_distribution<> rng_cr_shrub(0,15);

/*
// -----------------------------------------------------
// ... normal random number generator for soil depth ...
// -----------------------------------------------------
   boost::variate_generator<boost::mt19937, boost::normal_distribution<> >
      rng_soil_depth(rng, boost::normal_distribution<>(soil_depth_mean_, soil_depth_sd_));
   boost::math::normal snormal(1.35*33.15*M_PI/180.0,0.75*5.5*M_PI/180.0);

// --------------------------------------------------------
// ... normal random number generator for soil cohesion ...
// --------------------------------------------------------
   boost::variate_generator<boost::mt19937, boost::normal_distribution<> >
      rng_soil_cohesion(rng, boost::normal_distribution<>(soil_cohesion_mean_, soil_cohesion_sd_));

// --------------------------------------------------------
// ... normal random number generator for soil friction_angle ...
// --------------------------------------------------------
   boost::variate_generator<boost::mt19937, boost::normal_distribution<> >
      rng_soil_friction_angle(rng, boost::normal_distribution<>(soil_friction_angle_mean_, soil_friction_angle_sd_));
*/

// --------------------------------------------------------
// ... normal random number generator for landslide area ...
// --------------------------------------------------------
   boost::math::normal_distribution<> rng_area(area_mu_, area_sigma_);

// --------------------------------------------------------
// ... normal random number generator for length-to-width ratio ...
// --------------------------------------------------------
   boost::math::normal_distribution<> rng_l2w(l2w_mu_, l2w_sigma_);

// -----------------------------------------------------
// ... Inverse gamma distribution for landslide area ...
// -----------------------------------------------------
   //boost::math::inverse_gamma_distribution<> inv_gamma_dist(inv_gamma_rho_, inv_gamma_a_);

   // Plot density distribution as a check. Here using Tecplot 360 output format
   /*std::ofstream fout;
   fout.open("inv_gamma_dist.dat");
   fout << "VARIABLES = \"Area [km<sup>2</sup>]\" \"p [km<sup>-2</sup>]\"" << std::endl;
   fout << "DATASETAUXDATA rho = \"" << inv_gamma_rho_ << "\"" << std::endl;
   fout << "DATASETAUXDATA a = \"" << inv_gamma_a_ << "\"" << std::endl;
   fout << "DATASETAUXDATA s = \"" << inv_gamma_s_ << "\"" << std::endl;
   fout << "ZONE T = \"Inverse Gamma Probability Distribution\"" << std::endl;
   const auto area_max = 0.005; // [km2] 
   const auto num_intervals = 10000;
   for (auto i = 0; i < num_intervals; i++)
   {
      const auto area = static_cast<double>(i)/static_cast<double>(num_intervals)*area_max;
      fout << area << " " << boost::math::pdf(inv_gamma_dist, area + inv_gamma_s_) << std::endl;
   }
   fout.close();*/

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

   auto count = 0;
   unsigned int Fs200_pos, Fs800_pos, Pa200_pos, Pa400_pos, Mf300_pos, Mf600_pos, Cs150_pos;
   std::vector<double> Fs200, Fs800, Pa200, Pa400, Mf600, Mf300, Cs150;
   std::string line, word;

   getline(fin, line);
   std::stringstream ss(line);
   while (ss.good())
   {
      getline(ss,word,',');
      word.erase(std::remove_if(word.begin(),word.end(),
         [](auto const & c) -> bool {return !std::isalnum(c);}),word.end());
      
      if (word == "Pa400") Pa400_pos = count;
      else if (word == "Pa200") Pa200_pos = count;
      else if (word == "Fs800") Fs800_pos = count;
      else if (word == "Fs200") Fs200_pos = count;
      else if (word == "Cs150") Cs150_pos = count;
      else if (word == "MF600") Mf600_pos = count;
      else if (word == "MF300") Mf300_pos = count;
      count++;
   }

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

   std::vector<std::vector<double>> phi;
   std::vector<std::vector<double>> gamma;
   std::vector<std::vector<double>> ks;

   std::vector<double> phi1;
   std::vector<double> phi2;
   std::vector<double> gamma1;
   std::vector<double> ks1;
   std::vector<double> ks2;
   for (int i = 0; i < num_landslides; i++)
   {
      Landslide slide;

      phi1.push_back(quantile(rng_phi1, rng_uniform_01()));
      phi2.push_back(quantile(rng_phi2, rng_uniform_01()));
      gamma1.push_back(quantile(rng_gamma, rng_uniform_01()));
      ks1.push_back(quantile(rng_ks, rng_uniform_01()));
      ks2.push_back(quantile(rng_ks, rng_uniform_01()));

      slide.area_ = pow(10,quantile(rng_area, rng_uniform_01()));
      auto l2w = pow(10,quantile(rng_l2w, rng_uniform_01()));
      slide.width_ = sqrt(slide.area_ / l2w);
      slide.length_ = slide.width_ * l2w;
      landslide_.push_back(slide);

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

   for (auto i = 0; i < num_landslides; i++)
   {
      Raster friction_angle(soil_type_);
      Raster permeability(soil_type_);
      Raster depth(soil_depth_);
      Raster crl(dusaf_);
      Raster crb(soil_depth_);

      for (auto j = 0; j < soil_type_.attribute_.size(); j++)
      {
         if (soil_type_.attribute_.at(j))
         {
            friction_angle.attribute_.at(j) = phi.at((int)soil_type_.attribute_.at(j)-1).at(i);
            permeability.attribute_.at(j) = ks.at((int)soil_type_.attribute_.at(j)-1).at(i);
         }

         if (soil_depth_.attribute_.at(j))
         {
            for (auto k = 0; k < soil_id_.size(); k++)
            {
               if (soil_depth_.attribute_.at(j) == soil_id_.at(k))
               {
                  depth.attribute_.at(j) = z_.at(k).at(i);
                  break;
               }
            }
         }

         switch ((int)dusaf_.attribute_.at(i))
         {
            case 3211:
            case 3212:
            case 3221:
               crl.attribute_.at(i) = Cr_grassland_.at(i);
               break;
            case 332:
            case 333:
               crl.attribute_.at(i) = Cr_shrubland_.at(i);
               break;
            case 3121:
               crl.attribute_.at(i) = Crl_Pa400_.at(i);
               break;
            case 3122:
               crl.attribute_.at(i) = Crl_Pa200_.at(i);
               break;
            case 31111:
               crl.attribute_.at(i) = Crl_Fs800_.at(i);
               break;
            case 31121:
               crl.attribute_.at(i) = Crl_Fs200_.at(i);
               break;
            case 3114:
            case 222:
               crl.attribute_.at(i) = Crl_Cs150_.at(i);
               break;
            case 31311:
               crl.attribute_.at(i) = Crl_Mf600_.at(i);
               break;
            case 31321:
               crl.attribute_.at(i) = Crl_Mf300_.at(i);
               break;
            default:
               crl.attribute_.at(i) = dusaf_.attribute_.at(i);
               break;
         }

         if (soil_depth_.attribute_.at(i) >= 0.5) crb.attribute_.at(i) = 0;
         else crb.attribute_.at(i) = soil_depth_.attribute_.at(i);
      }

      auto m = TopModel_v3(permeability,depth);
   }

   auto finish_sli = std::chrono::high_resolution_clock::now(); // End time bi-linear interpolation
   std::chrono::duration<double> elapsed_sli = finish_sli - start_sli;
   std::cout << "Landslide generation elapsed time: " << elapsed_sli.count() << " s\n";

/*
// ----------------------------------------------
// ... Calculate slopes, twi, and related attributes ...
// ----------------------------------------------
   auto start_slope_ls = std::chrono::high_resolution_clock::now(); // Time loop

   for (auto & it : landslide_)
   {
      Point p = Point(it.x_,it.y_);
      double minx = std::max(it.x_ - it.width_, dem_.xllcorner_);
      double maxx = std::min(it.x_ + it.width_, dem_.xurcorner_);
      double miny = std::max(it.y_ - it.width_, dem_.yllcorner_);
      double maxy = std::min(it.y_ + it.width_, dem_.yurcorner_);
      int counter = 0;

      auto llIndex = FindRasterCellIndex(dem_, minx, miny);
      auto lrIndex = FindRasterCellIndex(dem_, maxx, miny);
      auto ulIndex = FindRasterCellIndex(dem_, minx, maxy);
      for (int j = 0; j <= (ulIndex - llIndex)/dem_.ncols_; j++)
      {
         for (int i = llIndex; i <= lrIndex; i++)
         {
            if (p.dist(dem_.points_.at(i + j*dem_.ncols_)) <= it.width_)
            {
               counter++;
               it.slope_ += slope_.attribute_.at(i + j*dem_.ncols_);
               it.twi_ += twi_.attribute_.at(i + j*dem_.ncols_);
            }
         }
      }

      if (counter == 0)
      {
         auto cellIndex = FindRasterCellIndex(dem_, it.x_, it.y_);
         it.slope_ = atan(slope_.attribute_.at(cellIndex)) * 180 / M_PI;
         it.twi_ = twi_.attribute_.at(cellIndex);
      } else
      {
         it.slope_ = atan(it.slope_ / counter) * 180 / M_PI;
         it.twi_ = it.twi_ / counter;
      }

      it.volume_ = it.area_ * it.soildepth_ * cos(it.GetSlopeRad());
      it.mass_ = it.volume_ * soil_density_;

      it.m_ = (rain_intensity_/transmissivity_) * it.twi_; // incorrect but in original R code
      //it.m_ = (rain_intensity_/transmissivity_) * exp(it.twi_) / cos(it.GetSlopeRad()); // correct
      if (it.m_ > 1.0) it.m_ = 1.0;

      it.pw_ = soil_density_ * gravity_ * it.soildepth_ * it.m_; // [Pa]
      if (it.m_ > 1) it.m_ = 1;
   }

   auto finish_slope_ls = std::chrono::high_resolution_clock::now(); // End time extraction
   std::chrono::duration<double> elapsed_slope_ls = finish_slope_ls - start_slope_ls;
   std::cout << "Extraction elapsed time: " << elapsed_slope_ls.count() << " s\n";
*/

   return true;

}

/*void Primula::AddTrees(const std::std::string & file)
{
// -----------------------------------------------------
// ... Inverse gamma distribution for root basal ...
// -----------------------------------------------------
   boost::math::inverse_gamma_distribution<> inv_gamma_dist_reinf(0.01220221, 0.004682089);

   auto start_trees = std::chrono::high_resolution_clock::now(); // Time loop

   unsigned int tree_count = 0;
   double sum_dbh = 0.0;

   TreeRead(file);

   for (auto & t : trees_)
   {
      tree_count++;
      sum_dbh += t.dbh_;

      Point p(t.x_, t.y_);
      for (auto & it : landslide_)
      {
         Point q(it.x_, it.y_);
         if (p.dist(q) < t.root_len_)
         {
            it.root_lateral_ += (276601 * t.dbh_ * boost::math::pdf(inv_gamma_dist_reinf,p.dist(q)/(t.root_len_)))/1000;
         }
      }
   }

   sum_dbh /= tree_count;
   for (auto & t : trees_)
      t.SetWeight(sum_dbh);
   std::sort(trees_.begin(), trees_.end(), [](const Tree & a, const Tree & b) -> bool
   {
      return a.weight_ < b.weight_;
   });
   if (trees_.size() % 2)
      veg_weight_ = trees_.at(trees_.size()/2).weight_;
   else veg_weight_ = (trees_.at((trees_.size()-1)/2).weight_ + trees_.at(trees_.size()/2).weight_)/2.0;
   veg_weight_ *= 0.1;

   auto finish_trees = std::chrono::high_resolution_clock::now(); // End time extraction
   std::chrono::duration<double> elapsed_trees = finish_trees - start_trees;
   std::cout << "Tree adding elapsed time: " << elapsed_trees.count() << " s\n";
}

void Primula::FindFOS()
{
// -----------------------------------------------------
// ... Inverse gamma distribution for root basal ...
// -----------------------------------------------------
   boost::math::inverse_gamma_distribution<> inv_gamma_dist_basal(1.3, 3.7);

   for (auto & it : landslide_)
   {
      // Root reinforcement values. Fixed and constant for the moment. To be computed later
//    it.root_lateral_ = 0.0;    // Should be user input so passed to object from main perhaps
      it.root_basal_ = boost::math::pdf(inv_gamma_dist_basal,it.soildepth_)*it.root_lateral_; // Fixed tree type for now
      
      auto Fd_parallel = gravity_ * sin(it.GetSlopeRad()) * (it.mass_ + veg_weight_ * it.area_);
      auto Fd_perpendicular = gravity_ * cos(it.GetSlopeRad()) * (it.mass_ + veg_weight_ * it.area_);
      auto Fr_basal = (it.area_ * (it.cohesion_ + it.root_basal_)) + ((Fd_perpendicular-(it.area_*it.pw_)) * tan(it.friction_angle_));
      auto Fr_lateral = it.circum_ / 2 * it.root_lateral_;
      auto Fr = Fr_basal + Fr_lateral;
      it.fos_ = Fr / Fd_parallel;
   }
}*/
