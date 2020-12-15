#pragma once

#include <CLI11.hpp>
#include <primula++.hpp>
#include <spdlog/spdlog.h>

class Config
{
public:
   Config(int argc, char **argv);

   Primula configModel();

private:
   std::string slopePath;
   std::string twiPath;
   std::string soilTypePath;
   std::string soilDepthPath;
   std::string dusafPath;
   std::string probslopePath;

   std::string outputPath;

   void dumpConfigFile(std::string path);

   CLI::App parser{"PRIMULA: Code that does stuff."};
};