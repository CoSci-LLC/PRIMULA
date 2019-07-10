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

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <memory>
#include <map>
#include <unordered_map>
#include <stdio.h>
#include <math.h>

#include "raster.hpp"
#include "landslide.hpp"
//#include "tree.hpp"

class Primula
{
   public:

   //==========================================================================
   // ... Constructors, Destructors ...
   //==========================================================================
      Primula();

      ~Primula();

   //==========================================================================
   // ... Public Accessor Functions ...
   //==========================================================================


   //==========================================================================
   // ... Public Member Functions ...
   //==========================================================================
      bool GenerateLandslides(const std::string & file, const unsigned int & num_landslides);
      bool ReadCSV(const std::string & file, const unsigned int & num_landslides);
      //void FindFOS();

      // ... Static Member Data ...
      static constexpr double gravity_ = 9.81;                  // [m/s^2]
      static constexpr double soil_density_ = 1834.8624;        // [kg/m^3] Schwarz M. R Code
      static constexpr double transmissivity_ = 0.0002644655;   // [m/s] Fixed for the moment. From R code
      static constexpr double rain_intensity_ = 0.001*70/3600;  // [m/s] >> Eventually should live outside of slideformap class <<

 
      // ... Member Data ...
      //Raster dem_;
      //Raster twi_;
      //Raster slope_;
      Raster soil_type_;
      Raster soil_depth_;
      std::vector<Landslide> landslide_;
      //std::vector<Tree> trees_;
      std::vector<int> soil_id_;
      std::vector<std::vector<double>> z_;

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
      double rainfall_ = 0.100;   // [m/day]
      
      // Inverse Gamma distribution parameters
      //double inv_gamma_rho_ = 1.6;          // BAFU landslide databank (2018)
      //double inv_gamma_a_ = 9.99E-05;       // BAFU landslide databank (2018)
      //double inv_gamma_s_ = 5.26163E-08;    // BAFU landslide databank (2018)
//    double inv_gamma_rho_ = 1.4;          // Malamud et al (2004)
//    double inv_gamma_a_ = 1.28E-03;       // [km2] Malamud et al (2004)
//    double inv_gamma_s_ = 1.32E-04;       // [km2] Malamud et al (2004)

      // Normal distribution parameters for soil depth
      //double soil_depth_mean_ = 1.14;       // [m]
      //double soil_depth_sd_ = 0.4;          // [m]

      // Normal distribution parameters for soil cohesion
      //double soil_cohesion_mean_ = 1500.0;  // [Pa]
      //double soil_cohesion_sd_ = 100.0;     // [Pa]

      // Normal distribution parameters for soil friction angle
      //double soil_friction_angle_mean_ = 38.0;    // [degree]
      //double soil_friction_angle_sd_ = 1.0;       // [degree]

      // Normal distribution parameters for landslide area
      double area_mu_ = 2.017;                // [m^2]
      double area_sigma_ = sqrt(0.176466);    // [m^2]

      // Normal distribution parameters for length-to-width ratio
      double l2w_mu_ = 0.1528;                // [dimensionless]
      double l2w_sigma_ = sqrt(0.037396);     // [dimensionless]


   //==========================================================================
   // ... Private Member Data ...
   //==========================================================================

   private:
      //bool TreeRead(const std::string & file);

};

#endif


