#pragma once
// ============================================================================
// Copyright (C) PRIMULA. All rights reserved.
//
// Authors: Denis Cohen-Corticchiato (DOHCC)
// Email:   denis.cohen@gmail.com
//
// This file is part of PRIMULA
// ============================================================================

class Landslide
{
public:
   Landslide(){};

   // Geometry
   double area_           = 0.0; // [m2]    Random, inverse gamma
   double length_         = 0.0; // [m]     Calculated
   double width_          = 0.0; // [m]     Calculated
};
