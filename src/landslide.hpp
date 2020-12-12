// ============================================================================
// Copyright (C) SlideforMAP++. All rights reserved.
//
// Authors: Denis Cohen-Corticchiato (DOHCC)
// Email:   denis.cohen@gmail.com
//
// This file is part of SlideforMAP++
// ============================================================================

#ifndef LANDSLIDE_HPP
#define LANDSLIDE_HPP

// #include <fstream>
// #include <iostream>
// #include <map>
// #include <math.h>
// #include <memory>
// #include <stdio.h>
// #include <string>
// #include <unordered_map>
// #include <vector>

#include <boost/math/distributions/normal.hpp>

class Landslide
{
public:
   //==========================================================================
   // ... Constructors, Destructors ...
   //==========================================================================
   Landslide();

   ~Landslide();

   //==========================================================================
   // ... Public Accessor Functions ...
   //==========================================================================

   double GetSlopeDeg();
   void SetFrictionAngle(const unsigned int &angle);

   //==========================================================================
   // ... Public Member Functions ...
   //==========================================================================

   // ... Static Member Data ...

   static constexpr double gravity_ = 9.81; // [m/s^2]
   // static constexpr double rhosoil_ = 1834.86;  // [kg/m^3] Fixed for the moment
   static constexpr double transmissivity_ = 0.0002644655; // [m/s] Fixed for the moment. From R code

   // ... Member Data ...

   double x_     = 0.0; // [m]                            Random, uniform
   double y_     = 0.0; // [m]                            Random, uniform
   double z_     = 0.0; // [m]                            Interpolated
   double slope_ = 0.0; // [m/m]                          Interpolated

   // Geometry
   double area_           = 0.0; // [m2]                           Random, inverse gamma
   double length_         = 0.0; // [m]                            Calculated
   double width_          = 0.0; // [m]                            Calculated
   double volume_         = 0.0; // [m3]                           Calculated
   double circum_         = 0.0; // [m]                            Calculated
   double mass_           = 0.0; // [kg]                           Calculated
   double soildepth_      = 0.0; // [m]                            Random, normal
   double cohesion_       = 0.0; // [Pa]                           Random, normal
   double friction_angle_ = 0.0; // [rad]                          Random, normal
   double twi_            = 0.0; // Topographic Wetness Index [-]  Interpolated
   double m_              = 0.0; //                                Calculated
   double pw_             = 0.0; // Pore-water pressure [Pa]       Calculated
   double root_basal_     = 0.0; // [Pa]                         Calculated
   double root_lateral_   = 0.0; // [Pa]                         User input
   double fos_            = 0.0; // Factor of Safety [-]           Calculated

   //==========================================================================
   // ... Private Member Data ...
   //==========================================================================

private:
};

#endif
