#include <CLI11.hpp>
#include <Config.hpp>
#include <filesystem>
#include <spdlog/fmt/ostr.h>

namespace fs = std::filesystem;

// Checks if a path exists and exits if it doesnt
void pathExists(std::string argName, std::string path);

Config::Config(int argc, char **argv)
{
   parser.set_config("--config");
   parser.config_formatter(std::make_shared<CLI::ConfigTOML>());

   // clang-format off
   parser.add_option("--slopePath",     this->slopePath,     "Path to elevation data")->check(CLI::ExistingFile);
   parser.add_option("--twiPath",       this->twiPath,       "Path to contributing area data")->check(CLI::ExistingFile);
   parser.add_option("--soilTypePath",  this->soilTypePath,  "Path to tree file (empty for no trees)")->check(CLI::ExistingFile);
   parser.add_option("--soilDepthPath", this->soilDepthPath, "Path to output directory")->check(CLI::ExistingFile);
   parser.add_option("--dusafPath",     this->dusafPath,     "Path to output directory")->check(CLI::ExistingFile);
   parser.add_option("--probslopePath", this->probslopePath, "Path to output directory")->check(CLI::ExistingFile);
   parser.add_option("--outputPath",    this->outputPath,    "Path to output directory")->check(CLI::NonexistentPath);
   // clang-format on

   std::ofstream outFile = std::ofstream(fs::path(this->outputPath) / fs::path("config.toml"));
   outFile << parser.config_to_str(true, true);
   outFile.close();
}

Primula Config::configModel()
{
   Primula model;

   model.slope_ = KiLib::Raster("../tests/malonno/slope_MALONNO.asc");
   model.twi_ = KiLib::Raster("../tests/malonno/twi_MALONNO.asc");
   model.soil_type_ = KiLib::Raster("../tests/malonno/soils_MALONNO_v2.asc");
   model.soil_depth_ = KiLib::Raster("../tests/malonno/soils_MALONNO.asc");
   model.dusaf_ = KiLib::Raster("../tests/malonno/dusaf_MALONNO.asc");
   model.probslope_ = KiLib::Raster("../tests/malonno/PROBSLOPE_MALONNO.asc");

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

void pathExists(std::string argName, std::string path)
{
   if (fs::exists(path) == false) {
      std::cerr << fmt::format("{}: the file `{}` does not exist!\n", argName, path);
      exit(EXIT_FAILURE);
   }
}