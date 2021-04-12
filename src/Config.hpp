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

#include <CLI11.hpp>
#include <primula++.hpp>
#include <spdlog/spdlog.h>

class Config
{
public:
   Config(int argc, char **argv);

   Primula configModel();

   std::string outputPath;

   const std::string &extension()
   {
      return this->defaultExtension;
   };

private:
   std::string slopePath;
   std::string twiPath;
   std::string soilTypePath;
   std::string landUsePath;

   std::string landCoverPath;
   std::string soilDepthPath;
   std::string physPropPath;

   std::string defaultExtension = ".tif";

   std::string rastOutPath;

   double rainfall_ = 0.160; // [m/day]

   size_t num_landslides = 100;
   size_t seed           = 333;

   bool calcFrd = true; //Default value

   void dumpConfigFile(std::string path);

   CLI::App parser{"PRIMULA++: A C++ version of PRobabilistIc MUltidimensional shallow Landslide Analysis"};
};
