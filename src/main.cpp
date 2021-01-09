// ============================================================================
// Copyright (C) model++. All rights reserved.
//
// Authors: Denis Cohen-Corticchiato (DOHCC)
// Email:   denis.cohen@gmail.com
//
// This file is part of PRIMULA++
// ============================================================================

#include <Config.hpp>
#include <filesystem>

namespace fs = std::filesystem;

int main(int argc, char **argv)
{
   Config config(argc, argv);

   Primula model = config.configModel();

   for (size_t i = 0; i < model.probslope_.nData; i++)
   {
      if (model.slope_(i) == model.slope_.nodata_value)
         model.probslope_(i) = model.probslope_.nodata_value;
      else if (model.probslope_(i) < 0)
         model.probslope_(i) = 1.1028656e-06;
   }

   model.ReadSoilDataset("../tests/malonno/Pedologia_25k_MALONNO.csv", "../tests/malonno/RootReinforcement.csv");
   model.GenerateSoilProperties();
   model.CalculateSafetyFactor();

   model.pr_failure_.writeToFile((fs::path(config.outputPath) / fs::path("prob_failure.asc")).string());

   return EXIT_SUCCESS;
}
