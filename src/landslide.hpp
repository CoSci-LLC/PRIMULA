/**
 *  Copyright (c) 2020-2021 CoSci LLC, USA <software@cosci-llc.com>
 *
 *  This file is part of PRIMULA.
 *
 *  PRIMULA is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  PRIMULA is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PRIMULA.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

class Landslide
{
public:
   Landslide(){};

   // Geometry
   double area_   = 0.0; // [m2]    Random, inverse gamma
   double length_ = 0.0; // [m]     Calculated
   double width_  = 0.0; // [m]     Calculated
};
