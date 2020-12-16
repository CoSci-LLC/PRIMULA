#include <CLI11.hpp>
#include <Config.hpp>
#include <filesystem>
#include <spdlog/fmt/ostr.h>

namespace fs = std::filesystem;

Config::Config(int argc, char **argv)
{
   parser.set_config("--config");
   parser.config_formatter(std::make_shared<CLI::ConfigTOML>());

   // clang-format off
   parser.add_option("--slopePath",     this->slopePath,     "Path to slope data")->check(CLI::ExistingFile);
   parser.add_option("--twiPath",       this->twiPath,       "Path to TWI data")->check(CLI::ExistingFile);
   parser.add_option("--soilTypePath",  this->soilTypePath,  "Path to soil type data")->check(CLI::ExistingFile);
   parser.add_option("--soilDepthPath", this->soilDepthPath, "Path to soil depth data")->check(CLI::ExistingFile);
   parser.add_option("--dusafPath",     this->dusafPath,     "Path to dusaf data")->check(CLI::ExistingFile);
   parser.add_option("--probslopePath", this->probslopePath, "Path to probslope data")->check(CLI::ExistingFile);
   parser.add_option("--seed",          this->seed,          "Seed for RNG", true);

   parser.add_option("--outputPath",    this->outputPath,    "Path to output directory");
   // clang-format on

   try {
      parser.parse(argc, argv);
   } catch (const CLI::ParseError &e) {
      parser.exit(e);
      std::cout << "\nDumping config to default_config.toml" << std::endl;
      this->dumpConfigFile("default_config.toml");
      exit(EXIT_FAILURE);
   }

   fs::create_directories(this->outputPath);
   std::ofstream outFile = std::ofstream(fs::path(this->outputPath) / fs::path("config.toml"));
   outFile << parser.config_to_str(true, true);
   outFile.close();
}

Primula Config::configModel()
{
   Primula model(this->seed);

   model.slope_      = KiLib::Raster(this->slopePath);
   model.twi_        = KiLib::Raster(this->twiPath);
   model.soil_type_  = KiLib::Raster(this->soilTypePath);
   model.soil_depth_ = KiLib::Raster(this->soilDepthPath);
   model.dusaf_      = KiLib::Raster(this->dusafPath);
   model.probslope_  = KiLib::Raster(this->probslopePath);

   return model;
}

// IMPORTANT:
// This wont work with more advanced CLI11 usage, keep that in mind when adding options.
void Config::dumpConfigFile(std::string path)
{
   std::ofstream outFile = std::ofstream(path);
   for (auto opt : this->parser.get_options()) {
      // Dont print help or config
      if (opt->get_lnames()[0] == "help" || opt->get_lnames()[0] == "config") {
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