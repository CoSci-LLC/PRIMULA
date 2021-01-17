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

   std::string landCoverPath;

private:
   std::string slopePath;
   std::string twiPath;
   std::string soilTypePath;
   std::string soilDepthPath;
   std::string dusafPath;
   std::string probslopePath;

   std::string defaultExtension = ".tif";

   size_t num_landslides = 100;
   size_t seed           = 69420;

   void dumpConfigFile(std::string path);

   CLI::App parser{"PRIMULA: Code that does stuff."};
};