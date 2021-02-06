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

   size_t num_landslides = 100;
   size_t seed           = 333;

   void dumpConfigFile(std::string path);

   CLI::App parser{"PRIMULA: Code that does stuff."};
};