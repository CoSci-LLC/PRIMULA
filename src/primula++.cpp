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
   

   std::cout << "Primula:ReadCSV ... Done" << std::endl;
   fin.close();
   return true;
}

void Primula::GenerateLandslides(const unsigned int & num_landslides)
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
// ... Generate LS coordinates using rng ...
// ... Generate random data ...
// ----------------------------------------------
   auto start_rng = std::chrono::high_resolution_clock::now(); // Time loop
   
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
   }
   phi.push_back(phi1);
   phi.push_back(phi2);
   gamma.push_back(gamma1);
   ks.push_back(ks1);
   ks.push_back(ks2);

   auto finish_rng = std::chrono::high_resolution_clock::now(); // End time bi-linear interpolation
   std::chrono::duration<double> elapsed_rng = finish_rng - start_rng;
   std::cout << "Landslide generation elapsed time: " << elapsed_rng.count() << " s\n";

/*
// ----------------------------------------------
// ... Bi-linear interpolation ...
// ... <TODO>
// ... SHOULD USE OMP
// ... </TODO>
// ----------------------------------------------
   auto start_bli = std::chrono::high_resolution_clock::now(); // Time bi-linear interpolation
   std::vector<double> raster_cell;
   for (auto & it : landslide_)
   {
//    auto rtn = FindRasterCell(dem_, it.x_, it.y_, raster_cell); // Returns {x1, x2, y1, y2, z11, z12, z21, z22}
      auto raster_cell = FindRasterCell(dem_, it.x_, it.y_);      // Returns {x1, x2, y1, y2, z11, z12, z21, z22}
//    auto cellCoord = FindRasterCellCoord(dem_, it.x_, it.y_);   // Returns {x1, x2, y1, y2}
//    auto cellIndex = FindRasterCellIndex(dem_, it.x_, it.y_);   // Returns {i, j, k, l}
//    if (!rtn) return EXIT_SUCCESS;
      it.z_ = BilinearInterpolation(raster_cell.at(4),
                                    raster_cell.at(5),
                                    raster_cell.at(6),
                                    raster_cell.at(7),
                                    raster_cell.at(0),
                                    raster_cell.at(1),
                                    raster_cell.at(2),
                                    raster_cell.at(3),
                                    it.x_,
                                    it.y_);
   }
   auto finish_bli = std::chrono::high_resolution_clock::now(); // End time bi-linear interpolation
   std::chrono::duration<double> elapsed_bli = finish_bli - start_bli;
   std::cout << "Bi-linear interpolation elapsed time: " << elapsed_bli.count() << " s\n";

   // Print to file
   std::std::string filename = "ls.xyz.dat";
// std::ofstream fout;
   fout.open(filename);
   fout << std::fixed << std::setprecision(4);
   for (auto & it : landslide_)
   {
      fout << it.x_ << " " << it.y_ << " " << it.z_ << std::endl;
   }
   fout.close();

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
