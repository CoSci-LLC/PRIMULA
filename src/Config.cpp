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

#include <CLI11.hpp>
#include <Config.hpp>
#include <filesystem>
#include <spdlog/fmt/ostr.h>

namespace fs = std::filesystem;

// LandCover.csv          : Soil probability data [Codice, Uso del suolo, Descrizione, Description, Min, Max]
// LANDUSE.tif            : Maps into LandCover.csv
// PhysicalProperties.csv : gamma, phi, cohesion, ks
// slope_rad.tif          : slope in radians
// SoilDepth.csv          : min/max soil depth for soil ID
// SOILS.tif              : Maps into soil depth
// TWI.tif                : Topographical wetness index

Config::Config(int argc, char **argv)
{
   parser.set_config("--config");
   parser.config_formatter(std::make_shared<CLI::ConfigTOML>());

   // clang-format off
   parser.add_option("--twiPath",         this->twiPath,          "Path to TWI Raster")->check(CLI::ExistingFile)->required();
   parser.add_option("--slopePath",       this->slopePath,        "Path to slope Raster")->check(CLI::ExistingFile)->required();
   parser.add_option("--soilTypePath",    this->soilTypePath,     "Path to soil type Raster")->check(CLI::ExistingFile)->required();
   parser.add_option("--landUsePath",     this->landUsePath,      "Path to land use Raster")->check(CLI::ExistingFile)->required();

   parser.add_option("--soilDepthPath",   this->soilDepthPath,    "Path to soil depth CSV")->check(CLI::ExistingFile)->required();
   parser.add_option("--landCoverPath",   this->landCoverPath,    "Path to the landcover CSV")->check(CLI::ExistingFile)->required();
   parser.add_option("--physPropPath",    this->physPropPath,     "Path to the Physical Properties CSV")->check(CLI::ExistingFile)->required();

   parser.add_option("--outputPath",      this->outputPath,       "Path to output directory")->required();

   parser.add_option("--rainfall",        this->rainfall_,        "Rainfall in meters/day", true);

   parser.add_option("--outputExtension", this->defaultExtension, "File extension for output rasters", true);
   parser.add_option("--seed",            this->seed,             "Seed for RNG", true);
   parser.add_option("--numLandslides",   this->num_landslides,   "The number of landslides to simulate", true);
   parser.add_option("--calcFrd",         this->calcFrd,          "Include slope-parallel passive force on downslope margin", true);
   // clang-format on

   try
   {
      parser.parse(argc, argv);
   }
   catch (const CLI::ParseError &e)
   {
      parser.exit(e);
      std::cout << "\nWriting out a config file, called 'default_config.toml', for you to fill in." << std::endl;
      this->dumpConfigFile("default_config.toml");
      exit(EXIT_FAILURE);
   }

   // Output raster name
   this->rastOutPath = (fs::path(this->outputPath) / fs::path("prob_failure")).string() + this->extension();

   // Create output dir
   fs::create_directories(this->outputPath);

   // Output configuration
   std::ofstream outFile = std::ofstream(fs::path(this->outputPath) / fs::path("config.toml"));
   outFile << parser.config_to_str(true, true);
   outFile.close();
}

Primula Config::configModel()
{
   spdlog::info("Configuring model");
   Primula model(this->num_landslides, this->seed);
   model.rainfall_ = this->rainfall_;

   spdlog::info("Loading Rasters");
   model.twi_       = KiLib::Raster(this->twiPath);
   model.slope_     = KiLib::Raster(this->slopePath);
   model.soil_type_ = KiLib::Raster(this->soilTypePath);
   model.landuse    = KiLib::Raster(this->landUsePath);

   // Make sure raster dimension agree
   for (const auto rast : {&model.slope_, &model.twi_, &model.soil_type_, &model.landuse})
   {
      if (rast->nRows != model.slope_.nRows)
      {
         spdlog::error("Raster row sizes dont agree!");
         exit(EXIT_FAILURE);
      }

      if (rast->nCols != model.slope_.nCols)
      {
         spdlog::error("Raster col sizes dont agree!");
         exit(EXIT_FAILURE);
      }
   }

   spdlog::info("Reading CSVs");
   model.ReadLandCover(this->landCoverPath);
   model.ReadSoilDepth(this->soilDepthPath);
   model.ReadPhysProps(this->physPropPath);

   model.validateData();

   spdlog::info("Done configuring model");

   return model;
}

// IMPORTANT:
// This wont work with more advanced CLI11 usage, keep that in mind when adding options.
void Config::dumpConfigFile(std::string path)
{
   std::ofstream outFile = std::ofstream(path);
   for (auto opt : this->parser.get_options())
   {
      // Dont print help or config
      if (opt->get_lnames()[0] == "help" || opt->get_lnames()[0] == "config")
      {
         continue;
      }

      // Print description
      fmt::print(outFile, "# {} | Type: {}\n", opt->get_description(), opt->get_type_name());

      // Print name and default value, wrap strings in quotes
      if (opt->get_type_name().find("TEXT") != std::string::npos)
         fmt::print(outFile, "{} = \"{}\"\n\n", opt->get_lnames()[0], opt->get_default_str());
      else
         fmt::print(outFile, "{} = {}\n\n", opt->get_lnames()[0], opt->get_default_str());
   }
   outFile.close();
}
